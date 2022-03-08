extern crate nimble;

use nimble::process::{bam, fastq};
use nimble::reference_library;
use nimble::utils;

use clap::{load_yaml, App};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

fn main() {
    // Parse command line arguments based on the yaml schema
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    let reference_json_path = matches.value_of("reference").unwrap();
    let output_path = matches.value_of("output").unwrap();
    let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
    let num_cores = matches
        .value_of("num_cores")
        .unwrap_or("1")
        .parse::<usize>()
        .expect("Error -- please provide an integer value for the number of cores");
    let debug_file = matches.value_of("log").unwrap_or("").to_owned();
    let debug_file = if debug_file == "" {
        None
    } else {
        Some(debug_file)
    };
    let alignment_file = matches.value_of("alignment").unwrap_or("").to_owned();
    let alignment_file = if alignment_file == "" {
        None
    } else {
        Some(alignment_file)
    };

    println!("Loading and preprocessing reference data");

    // Read library alignment config info, reference library metadata, and sequences from library json
    let (align_config, reference_metadata) =
        reference_library::get_reference_library(Path::new(reference_json_path));

    // Generate error-checked vectors of seqs and names for the debrujin index
    let (reference_seqs, reference_seqs_rev, reference_names) =
        utils::validate_reference_pairs(&reference_metadata);

    // Create debruijn index of the reference library
    let reference_index_forward =
        debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
            &reference_seqs,
            &reference_names,
            &HashMap::new(),
            num_cores,
        )
        .expect("Error -- could not create pseudoaligner index of the reference library");

    // Create debruijn index of the reverse-comp of the reference library
    let reference_index_reverse =
        debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
            &reference_seqs_rev,
            &reference_names,
            &HashMap::new(),
            num_cores,
        )
        .expect("Error -- could not create reverse pseudoaligner index of the reference library");

    let reference_index = (reference_index_forward, reference_index_reverse);

    println!("Loading read sequences");

    if utils::is_fastq(input_files[0]) {
        fastq::process(
            input_files,
            &reference_index,
            &reference_metadata,
            &align_config,
            output_path,
            debug_file,
            alignment_file
        );
    } else {
        bam::process(
            input_files,
            &reference_index,
            &reference_metadata,
            &align_config,
            output_path,
            debug_file,
            alignment_file
        );
    };

    // Ensure we've created the output_path file
    OpenOptions::new()
        .create(true)
        .append(true)
        .open(output_path)
        .expect("Error -- could not create results file");
}
