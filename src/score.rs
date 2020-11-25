use crate::align;
use crate::reference_library;
use crate::utils;

use reference_library::ReferenceMetadata;
use std::io::Error;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

/* Takes a list of sequences and optionally reverse sequences, a reference library index, reference library metadata,
 * and an aligner configuration object, and returns a vector of scores and relative match percentages generated from an alignment
 * of the sequences to the reference library. */
pub fn score<I>(sequences: I, reverse_sequences: Option<I>,
  reference_index: Pseudoaligner<debruijn_mapping::config::KmerType>,
  reference_metadata: &ReferenceMetadata, align_config: align::AlignFilterConfig) -> Vec<(String, i32, f32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>
  {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, &reference_metadata, &align_config);

  // Remove scores below the score threshold
  let reference_scores: Vec<(String, i32)> = reference_scores.into_iter().filter(|(_, val)| val > &align_config.score_filter).collect();

  // Append match percentages for each reference result
  let num_reads: i32 = reference_scores.iter().map(|(_, val)| val).sum();
  let results = utils::append_match_percent(reference_scores, num_reads as usize);

  utils::sort_score_vector(results)
}