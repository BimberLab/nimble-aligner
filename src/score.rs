use crate::align;
use crate::reference_library;
use crate::utils;

use reference_library::ReferenceMetadata;
use std::io::Error;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

pub fn score<I>(sequences: I, reverse_sequences: Option<I>,
  reference_index: Pseudoaligner<debruijn_mapping::config::KmerType>,
  reference_metadata: &ReferenceMetadata, align_config: align::AlignFilterConfig) -> Vec<(String, i32, f32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>
  {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, &reference_metadata, &align_config);


  // Begin report formatting pipeline
  let reference_scores: Vec<(String, i32)> = reference_scores.into_iter().filter(|(_, val)| val > &align_config.score_filter).collect();
  //let results: Vec<(String, f32)> = reference_scores.into_iter().map(|(name, val)| (name, val as f32)).collect(); // Debug line for looking at raw scores
  let num_reads: i32 = reference_scores.iter().map(|(_, val)| val).sum();
  let results = utils::convert_scores_to_percentage(reference_scores, num_reads as usize);

  utils::sort_score_vector(results)
}