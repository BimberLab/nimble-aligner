mod filter;
mod utils;
mod align;

use std::path;
use std::fs::File;
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;

fn main() {
  println!("Loading and preprocessing data");

  // Parse arguments based on the yaml schema
  let yaml = load_yaml!("cli.yml");
  let matches = App::from_yaml(yaml).get_matches();

  let library = matches.value_of("library").unwrap();
  let library_fasta = matches.value_of("library_fasta").unwrap();
  let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
  let num_cores = matches.value_of("num_cores").unwrap_or("1").parse::<usize>().expect("Error -- please provide an integer value for the number of cores");

  // Parameters for alignment, alignment filtering, and report filtering
  const SCORE_THRESHOLD: usize = 60;
  const REFERENCE_GENOME_SIZE: usize = 1209;
  const READS_SIZE: usize = 2786342;
  const PERCENT_THRESHOLD: f32 = 0.0;
  const NUM_MISMATCHES: usize = 150;
  const DISCARD_DIFFERING_READ_PAIRS: bool = false;
  const DISCARD_NONZERO_MISMATCH: bool = false;
  const DISCARD_MULTIPLE_MATCHES: bool = true;

  // Get iterators to records in the reference genome file and library
  let reference_genome  = fasta::Reader::from_file(library_fasta).expect(
    "Error -- could not read reference genome"
  ).records();

  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library for alignment"
  )).into_records();

  // Parse out reference/reference name pairs that cannot be read for whatever reason
  let (reference_seqs, reference_names) = utils::get_valid_reference_pairs(reference_genome, reference_library);

  // Create debruijn-mapped index of the reference library
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    num_cores
  ).expect("Error -- could not create pseudoaligner index of the reference library");

  println!("Reading sequences");

  /* Get error-checked iterators to the sequences that will be aligned to the reference from the
   * sequence genome file(s) */
  let sequences = utils::get_error_checked_fastq_reader(input_files[0]);

  // Only get reverse sequences if a file is provided
  let reverse_sequences = if input_files.len() > 1 {
    println!("Reading reverse sequences");
    Some(utils::get_error_checked_fastq_reader(input_files[1]))
  } else {
    None
  };

  println!("Pseudo-aligning reads to reference index");

  // Configure aligner
  let align_config = align::AlignFilterConfig {
    reference_genome_size: REFERENCE_GENOME_SIZE,
    score_threshold: SCORE_THRESHOLD,
    num_mismatches: NUM_MISMATCHES,
    discard_differing_read_pairs: DISCARD_DIFFERING_READ_PAIRS,
    discard_nonzero_mismatch: DISCARD_NONZERO_MISMATCH,
    discard_multiple_matches: DISCARD_MULTIPLE_MATCHES
  };

  // Perform filtered pseudoalignment 
  let reference_scores = align::score(sequences, reverse_sequences, reference_index, align_config);

  println!("Filtering results by lineage");

  // Create reference library iterator to filter the scores by lineage
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library for filtration"
  )).into_records();

  // Post-alignment filtration pipeline
  let results = filter::report::collapse_results_by_lineage(reference_library, reference_scores.iter());
  let results = utils::convert_scores_to_percentage(results, READS_SIZE);
  let results = filter::report::threshold_percentage(results, PERCENT_THRESHOLD);

  println!("Writing results to file");

  utils::write_to_tsv(results);

  print!("Output results written to ./results.tsv")
}