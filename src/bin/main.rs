extern crate nimble;

use nimble::align::LibraryChemistry;
use nimble::process::{bam, fastq};
use nimble::reference_library;
use nimble::utils;

use clap::{load_yaml, App};
use std::collections::HashMap;
use std::path::Path;

fn main() {
    let yaml = load_yaml!("cli.yml");

    // Parse command line arguments based on the yaml schema
    let matches = App::from_yaml(yaml).get_matches();

    // Set of reference input paths
    let reference_json_paths: Vec<String> = matches
        .values_of("reference")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

    // Set of alignment output paths
    let output_paths: Vec<String> = matches
        .values_of("output")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

    // Number of cores to use for sequence alignment
    let num_cores = matches
        .value_of("num_cores")
        .unwrap_or("1")
        .parse::<usize>()
        .expect("Error -- please provide an integer value for the number of cores");

    // Flag defining library chemistry, which impacts orientation filters
    let strand_filter = matches.value_of("strand_filter").unwrap_or("unstranded");
    let strand_filter = match strand_filter {
        "unstranded" => LibraryChemistry::Unstranded,
        "fiveprime" => LibraryChemistry::FivePrime,
        "threeprime" => LibraryChemistry::ThreePrime,
        "none" => LibraryChemistry::None,
        _ => panic!("Could not parse strand_filter option."),
    };

    // Set of input files (one or two .fastqs, or one .bam)
    let input_files: Vec<String> = matches
        .values_of("input")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

    let force_bam_paired = matches.is_present("force_bam_paired");

    // The first input file, used to determine which pipeline to run
    let first_input_file = &input_files[0];
    let file_path = Path::new(first_input_file);
    let file_extension_sequence = file_path
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .unwrap_or("")
        .to_lowercase();
    let is_fastq_gz = file_path
        .file_name()
        .and_then(std::ffi::OsStr::to_str)
        .map(|name| name.ends_with(".fastq.gz"))
        .unwrap_or(false);

    let mut reference_indices = Vec::new();
    let mut references = Vec::new();
    let mut aligner_configs = Vec::new();

    // Parse the trim option if provided
    let trim_option = matches.value_of("trim");
    let mut trim_pairs = Vec::new();
    if let Some(trim_option) = trim_option {
        trim_pairs = trim_option.split(',')
            .map(|s| {
                let mut parts = s.split(':');
                let length = parts.next().unwrap().parse::<usize>().expect("Invalid length");
                let strictness = parts.next().unwrap().parse::<f64>().expect("Invalid strictness");
                (length, strictness)
            })
            .collect::<Vec<_>>();

        if trim_pairs.len() != reference_json_paths.len() {
            panic!("The number of trim options does not match the number of reference libraries");
        }
    }

    // Produce two debruijn graphs per reference library
    for (i, reference_json_path) in reference_json_paths.iter().enumerate() {
        println!(
            "Loading and preprocessing reference data for {}",
            reference_json_path
        );

        // Load the reference library into memory
        let (mut aligner_config, reference) =
            reference_library::get_reference_library(
                Path::new(&reference_json_path),
                strand_filter,
            );

        // Update trim options if provided
        if let Some((length, strictness)) = trim_pairs.get(i) {
            aligner_config.trim_target_length = *length;
            aligner_config.trim_strictness = *strictness;
            println!("Manually setting trim settings for library {} | target length: {}, strictness: {}",
                reference_json_path, aligner_config.trim_target_length, aligner_config.trim_strictness);
        }

        // Get sequences and feature names from the reference library for producing the library index
        let (reference_sequences, reference_feature_names) =
            utils::get_reference_sequence_data(&reference);

        // Produce reference index from the library sequence data
        let reference_index =
            debruijn_mapping::build_index::build_index::<debruijn::kmer::Kmer30>(
                &reference_sequences,
                &reference_feature_names,
                &HashMap::new(),
                num_cores,
            )
            .expect("Error -- could not create pseudoaligner index of the reference library");

        reference_indices.push(reference_index);
        references.push(reference);
        aligner_configs.push(aligner_config);
    }

    println!("Loading read sequences and aligning");

    if is_fastq_gz || file_extension_sequence == "fastq" {
        println!("Processing as FASTQ file");
        fastq::process(
            input_files,
            reference_indices,
            references,
            aligner_configs,
            output_paths
        );
    } else if file_extension_sequence == "bam" {
        println!("Processing as BAM file");
        bam::process(
            input_files,
            reference_indices,
            references,
            aligner_configs,
            output_paths,
            num_cores,
            force_bam_paired
        );
    } else {
        panic!("Unsupported file format: {}", file_extension_sequence);
    }

    println!("Alignment successful, terminating.")
}
