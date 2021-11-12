extern crate nimble;

use std::path::Path;
use std::fs::OpenOptions;
use std::collections::HashMap;

use clap::{load_yaml, App};

use nimble::process::{bam, fastq};
use nimble::reference_library;
use nimble::utils;

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

    println!("Loading and preprocessing reference data");

    // Read the reference library's aligner configuration and data from the given .json file
    let (aligner_config, reference_metadata) =
        reference_library::get_reference_library(Path::new(reference_json_path));

    // Create error-checked vectors of reference sequences and names, which will be passed to the debrujin index
    // Since we create a debruijn index for both forward and reverse versions of the library, we
    // get two reference vectors
    let (reference_seqs, reference_seqs_rev, reference_names) =
        utils::get_valid_reference_sequence_lists(&reference_metadata);

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

    // Pack the indices into one variable and pass it into one of two processing pipelines based on the
    // file extension of the first input file
    let reference_index = (reference_index_forward, reference_index_reverse);

    println!("Loading read sequences");

    if utils::is_fastq(input_files[0]) {
        let results = fastq::process(
            input_files,
            &reference_index,
            &reference_metadata,
            &aligner_config,
            debug_file
        );

        utils::write_to_tsv(utils::filter_scores(results, &aligner_config.score_filter), None, true, output_path);
    } else {
        let (results, cell_barcodes) = bam::process(
            input_files,
            &reference_index,
            &reference_metadata,
            &aligner_config,
            debug_file
        );

        utils::write_to_tsv(results, Some(cell_barcodes), false, output_path);
    };

    // Ensure we've created the output_path file
    OpenOptions::new()
        .create(true)
        .append(true)
        .open(output_path)
        .expect("Error -- could not create results file");
}
