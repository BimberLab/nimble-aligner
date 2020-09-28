extern crate immuno_genotyper;

use immuno_genotyper::score;
use immuno_genotyper::utils;
use immuno_genotyper::reference_library;

use std::path::Path;
use std::fs::File;
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;

fn main() {
  println!("Loading and preprocessing data");

  // Parse arguments based on the yaml schema
  let yaml = load_yaml!("cli.yml");
  let matches = App::from_yaml(yaml).get_matches();

  let library_path = matches.value_of("library").unwrap();
  let library_fasta = matches.value_of("library_fasta").unwrap();
  let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
  let num_cores = matches.value_of("num_cores").unwrap_or("1").parse::<usize>().expect("Error -- please provide an integer value for the number of cores");

  // Get iterators to records in the reference genome file and library
  let reference_genome  = fasta::Reader::from_file(library_fasta).expect(
    "Error -- could not read reference genome"
  ).records();

  let library = utils::get_tsv_reader(File::open(Path::new(library_path)).expect(
    "Error -- could not read reference library for alignment"
  )).into_records();

  // Parse out reference/reference name pairs that cannot be read for whatever reason
  let (reference_seqs, reference_names) = utils::get_valid_reference_pairs(reference_genome, library);

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
  let align_config = reference_library::get_reference_library(Path::new(library_path));

  let results = score::score(sequences, reverse_sequences, reference_index, reference_library, align_config, GROUP_COLUMN);

  println!("Writing results to file");

  utils::write_to_tsv(results);

  print!("Output results written to ./results.tsv");
}