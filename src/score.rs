use crate::filter;
use crate::align;
use crate::reference_library;
use crate::utils;

use reference_library::ReferenceMetadata;
use std::io::Error;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

pub fn score<I>(sequences: I, reverse_sequences: Option<I>,
  reference_index: Pseudoaligner<debruijn_mapping::config::KmerType>,
  reference_metadata: ReferenceMetadata, align_config: align::AlignFilterConfig) -> Vec<(String, f32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>
  {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, reference_metadata, &align_config);
  let reference_scores_len = reference_scores.len();

  // Begin report formatting pipeline
  let results = utils::convert_scores_to_percentage(reference_scores, reference_scores_len);
  let results = filter::report::threshold_percentage(results, align_config.percent_threshold);
  utils::sort_score_vector(results)
}