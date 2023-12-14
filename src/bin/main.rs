extern crate nimble;

use nimble::align::StrandFilter;
use nimble::process::{bam, fastq};
use nimble::reference_library;
use nimble::utils;

use clap::{load_yaml, App};
use std::collections::HashMap;
use std::path::Path;

use nimble::ALLOCATOR;

fn main() {

    let yaml = load_yaml!("cli.yml");

    // Parse command line arguments based on the yaml schema
    let matches = App::from_yaml(yaml).get_matches();

    // Modified to handle multiple reference libraries
    let reference_json_paths: Vec<String> = matches
        .values_of("reference")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

    let output_paths: Vec<String> = matches
        .values_of("output")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

    let _debug_files: Vec<Option<String>> =
        matches.values_of("log").map_or_else(Vec::new, |values| {
            values
                .map(|v| if v == "" { None } else { Some(v.to_string()) })
                .collect()
        });

    let _alignment_files: Vec<Option<String>> =
        matches
            .values_of("alignment")
            .map_or_else(Vec::new, |values| {
                values
                    .map(|v| if v == "" { None } else { Some(v.to_string()) })
                    .collect()
            });

    let num_cores = matches
        .value_of("num_cores")
        .unwrap_or("1")
        .parse::<usize>()
        .expect("Error -- please provide an integer value for the number of cores");

    let hard_memory_limit = matches
        .value_of("hard_memory_limit")
        .unwrap_or("0")
        .parse::<usize>()
        .expect("Error -- please provide a positive integer for the hard memory limit");

    let strand_filter = matches.value_of("strand_filter").unwrap_or("unstranded");
    let strand_filter = match strand_filter {
        "unstranded" => StrandFilter::Unstranded,
        "fiveprime" => StrandFilter::FivePrime,
        "threeprime" => StrandFilter::ThreePrime,
        "none" => StrandFilter::None,
        _ => panic!("Could not parse strand_filter option."),
    };

    let input_files: Vec<String> = matches
        .values_of("input")
        .unwrap()
        .map(|v| v.to_string())
        .collect();

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

    // Set memory allocator hard limit if it was passed
    if hard_memory_limit > 0 {
        ALLOCATOR.set_limit(hard_memory_limit * 1024 * 1024).unwrap();
        println!("Set memory allocator hard limit to {} bytes", hard_memory_limit * 1024 * 1024);
    }

    let mut reference_indices = Vec::new();
    let mut reference_metadata = Vec::new();
    let mut align_config = Vec::new();
    let mut reference_seqs = Vec::<Vec<String>>::new();
    let mut reference_names = Vec::<Vec<String>>::new();

    for reference_json_path in reference_json_paths {
        println!(
            "Loading and preprocessing reference data for {}",
            reference_json_path
        );

        let (align_config_thread, reference_metadata_thread) =
            reference_library::get_reference_library(
                Path::new(&reference_json_path),
                strand_filter,
            );

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
