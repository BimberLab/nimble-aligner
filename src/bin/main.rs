extern crate nimble;

use nimble::score;
use nimble::utils;
use nimble::reference_library;

use std::path::Path;
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;

fn main() {
  // Parse command line arguments based on the yaml schema
  let yaml = load_yaml!("cli.yml");
  let matches = App::from_yaml(yaml).get_matches();

  let reference_json_path = matches.value_of("reference").unwrap();
  let reference_fasta = matches.value_of("reference_fasta").unwrap();
  let output_path = matches.value_of("output").unwrap();
  let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
  let num_cores = matches.value_of("num_cores").unwrap_or("1").parse::<usize>().expect("Error -- please provide an integer value for the number of cores");

  println!("Loading and preprocessing reference data");

  // Read library alignment config info and reference library metadata from library json
  let (align_config, reference_metadata) = reference_library::get_reference_library(Path::new(reference_json_path));

  // Get iterator to records in the reference genome .fasta
  let reference_genome  = fasta::Reader::from_file(reference_fasta).expect(
    "Error -- could not read reference genome"
  ).records();

  // Generate error-checked vectors of seqs and names for the debrujin index
  let (reference_seqs, reference_names) = utils::validate_reference_pairs(reference_genome, reference_metadata.columns[reference_metadata.nt_sequence_idx].iter());

  // Create debruijn index of the reference library
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    num_cores
  ).expect("Error -- could not create pseudoaligner index of the reference library");

  println!("Loading read sequences");

  /* Get error-checked iterators to the sequences that will be aligned to the reference from the
   * sequence genome file(s) */
  let sequences = utils::get_error_checked_fastq_reader(input_files[0]);

  // Only get reverse sequences if a reverse sequence file is provided
  let reverse_sequences = if input_files.len() > 1 {
    println!("Reading reverse sequences");
    Some(utils::get_error_checked_fastq_reader(input_files[1]))
  } else {
    None
  };

  println!("Pseudo-aligning reads to reference index");

  // Perform alignment and filtration using the score package
  let results = score::score(sequences, reverse_sequences, reference_index, &reference_metadata, align_config);

  println!("Writing results to file");

  utils::write_to_tsv(results, output_path);

  print!("Output results written to output path");
}