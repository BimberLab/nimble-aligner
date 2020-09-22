use crate::filter;
use crate::align;

use std::io::{Error, Read};
use csv::StringRecordsIntoIter;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

pub fn score<I, R>(sequences: I, reverse_sequences: Option<I>,
  reference_index: Pseudoaligner<debruijn_mapping::config::KmerType>,
  library: StringRecordsIntoIter<R>, align_config: align::AlignFilterConfig, group_column: usize) -> Vec<(String, f32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>,
    R: Read
  {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, align_config);

  println!("Filtering results by lineage");

  // Post-alignment filtration pipeline
  let results = filter::report::collapse_results_by_lineage(library, reference_scores.iter(), group_column);

  let results = results.into_iter().map(|v| (v.0, v.1 as f32)).collect();

  //let results = utils::convert_scores_to_percentage(results, READS_SIZE);
  filter::report::threshold_percentage(results, 0.0)
}