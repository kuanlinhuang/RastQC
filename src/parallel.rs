use crate::config::FastQCConfig;
use crate::io::{Sequence, SequenceReader};
use crate::modules::{merge_module_sets, ModuleFactory, QCModule};
use anyhow::Result;
use rayon::prelude::*;
use std::path::Path;

const CHUNK_SIZE: usize = 100_000;
const MIN_SEQS_FOR_PARALLEL: usize = 200_000;

/// Process a file using intra-file parallelism.
///
/// Reads all sequences into memory, splits into chunks, processes each
/// chunk with independent module instances in parallel, then merges
/// results by combining module accumulator states.
///
/// This gives real speedup for large files by distributing the
/// CPU-intensive per-sequence processing across multiple threads.
pub fn process_file_parallel(
    path: &Path,
    config: &FastQCConfig,
) -> Result<(Vec<Box<dyn QCModule>>, u64)> {
    let mut reader = SequenceReader::open(path)?;

    // Phase 1: Read all sequences into memory
    let mut all_sequences: Vec<Sequence> = Vec::with_capacity(CHUNK_SIZE);
    while let Some(seq) = reader.next_sequence()? {
        all_sequences.push(seq);
    }

    let total_count = all_sequences.len() as u64;

    if total_count == 0 {
        let modules = ModuleFactory::create_modules(config);
        return Ok((modules, 0));
    }

    // If the file is small, process sequentially (no overhead from chunking)
    if all_sequences.len() < MIN_SEQS_FOR_PARALLEL {
        let mut modules = ModuleFactory::create_modules(config);
        for seq in &all_sequences {
            for module in modules.iter_mut() {
                module.process_sequence(seq);
            }
        }
        for module in modules.iter_mut() {
            module.calculate_results(config);
        }
        return Ok((modules, total_count));
    }

    // Phase 2: Process chunks in parallel with independent module sets
    let num_threads = rayon::current_num_threads();
    let chunk_size = (all_sequences.len() + num_threads - 1) / num_threads;

    let mut chunk_modules: Vec<Vec<Box<dyn QCModule>>> = all_sequences
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut modules = ModuleFactory::create_modules(config);
            for seq in chunk {
                for module in modules.iter_mut() {
                    module.process_sequence(seq);
                }
            }
            modules
        })
        .collect();

    // Phase 3: Merge all chunk results into the first chunk
    if chunk_modules.len() > 1 {
        // Take the first chunk as the accumulator
        let (first, rest) = chunk_modules.split_at_mut(1);
        let target = &mut first[0];

        for source in rest.iter_mut() {
            merge_module_sets(target, source);
        }
    }

    let mut final_modules = chunk_modules.into_iter().next().unwrap();

    // Phase 4: Calculate results on the merged modules
    for module in final_modules.iter_mut() {
        module.calculate_results(config);
    }

    Ok((final_modules, total_count))
}

/// Check if a file is large enough to benefit from parallel processing.
pub fn should_use_parallel(path: &Path) -> bool {
    path.metadata()
        .map(|m| m.len() > 50 * 1024 * 1024) // > 50 MB
        .unwrap_or(false)
}
