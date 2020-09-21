use crate::filter;
use crate::align;
use crate::utils;

use std::path;
use std::io::Error;
use std::fs::File;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;

pub fn score(sequences: impl Iterator<Item = Result<DnaString, Error>>, reverse_sequences: Option<impl Iterator<Item = Result<DnaString, Error>>>,
  reference_index: Pseudoaligner<debruijn_mapping::config::KmerType>, library: &str, align_config: align::AlignFilterConfig) -> Vec<(String, f32)> {

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, align_config);

  println!("Filtering results by lineage");

  // Create reference library iterator to filter the scores by lineage
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library for filtration"
  )).into_records();

  // Post-alignment filtration pipeline
  let results = filter::report::collapse_results_by_lineage(reference_library, reference_scores.iter());
  let results = results.into_iter().map(|v| (v.0, v.1 as f32)).collect();

  //let results = utils::convert_scores_to_percentage(results, READS_SIZE);
  filter::report::threshold_percentage(results, 10.0)
}