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
    I: Iterator<Item = Result<DnaString, Error>>,
  {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, align_config, reference_metadata);

  // TODO REMOVE
  let results = reference_scores.into_iter().map(|v| (v.0, v.1 as f32)).collect();

  //let results = utils::convert_scores_to_percentage(results, READS_SIZE);
  filter::report::threshold_percentage(results, 0.0)
}