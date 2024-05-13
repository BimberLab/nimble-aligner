extern crate nimble;

use nimble::align::StrandFilter;
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
        "unstranded" => StrandFilter::Unstranded,
        "fiveprime" => StrandFilter::FivePrime,
        "threeprime" => StrandFilter::ThreePrime,
        "none" => StrandFilter::None,
        _ => panic!("Could not parse strand_filter option."),
    };

    // Set of input files (one or two .fastqs, or one .bam)
    let input_files: Vec<String> = matches
        .values_of("input")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

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
    let mut reference_metadata = Vec::new();
    let mut align_config = Vec::new();
    let mut reference_seqs = Vec::<Vec<String>>::new();
    let mut reference_names = Vec::<Vec<String>>::new();

    // Produce two debruijn graphs per reference library
    for reference_json_path in reference_json_paths {
        println!(
            "Loading and preprocessing reference data for {}",
            reference_json_path
        );

        // Load the reference library into memory
        let (align_config_thread, reference_metadata_thread) =
            reference_library::get_reference_library(
                Path::new(&reference_json_path),
                strand_filter,
            );


        // TODO refactor + test all below
        let (reference_seqs_thread, reference_seqs_rev_thread, reference_names_thread) =
            utils::validate_reference_pairs(&reference_metadata_thread);

        let reference_index_forward =
            debruijn_mapping::build_index::build_index::<debruijn::kmer::Kmer30>(
                &reference_seqs_thread,
                &reference_names_thread,
                &HashMap::new(),
                num_cores,
            )
            .expect("Error -- could not create pseudoaligner index of the reference library");

        let reference_index_reverse = debruijn_mapping::build_index::build_index::<
            debruijn::kmer::Kmer30,
        >(
            &reference_seqs_rev_thread,
            &reference_names_thread,
            &HashMap::new(),
            num_cores,
        )
        .expect("Error -- could not create reverse pseudoaligner index of the reference library");

        reference_indices.push((reference_index_forward, reference_index_reverse));
        reference_metadata.push(reference_metadata_thread);
        align_config.push(align_config_thread);
        reference_seqs.push(reference_seqs_thread.into_iter().map(|s| s.to_string()).collect());
        reference_names.push(reference_names_thread);
    }

    println!("Loading read sequences and aligning");

    if is_fastq_gz || file_extension_sequence == "fastq" {
        println!("Processing as FASTQ file");
        fastq::process(
            input_files,
            reference_indices,
            reference_metadata,
            align_config,
            output_paths,
            num_cores,
        );
    } else if file_extension_sequence == "bam" {
        println!("Processing as BAM file");
        bam::process(
            input_files,
            reference_indices,
            reference_metadata,
            align_config,
            output_paths,
            num_cores,
        );
    } else {
        panic!("Unsupported file format: {}", file_extension_sequence);
    }

    println!("Alignment successful, terminating.")
}
