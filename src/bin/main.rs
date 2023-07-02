extern crate nimble;

use nimble::align::StrandFilter;
use nimble::process::bam;
use nimble::reference_library;
use nimble::utils;

use clap::{load_yaml, App};
use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

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

    // Thread-safe collections to store setup data
    let reference_indices = Arc::new(Mutex::new(Vec::new()));
    let reference_metadata = Arc::new(Mutex::new(Vec::new()));
    let align_config = Arc::new(Mutex::new(Vec::new()));
    let reference_seqs = Arc::new(Mutex::new(Vec::<Vec<String>>::new()));
    let reference_names = Arc::new(Mutex::new(Vec::<Vec<String>>::new()));

    let mut handles = vec![];

    for reference_json_path in reference_json_paths {
        // Clone Arcs to move into the thread
        let reference_indices = Arc::clone(&reference_indices);
        let reference_metadata = Arc::clone(&reference_metadata);
        let align_config = Arc::clone(&align_config);
        let reference_seqs = Arc::clone(&reference_seqs);
        let reference_names = Arc::clone(&reference_names);

        let strand_filter = strand_filter.clone();

        let handle = thread::spawn(move || {
            println!(
                "Loading and preprocessing reference data for {}",
                reference_json_path
            );

            // Read library alignment config info, reference library metadata, and sequences from library json
            let (align_config_thread, reference_metadata_thread) =
                reference_library::get_reference_library(
                    Path::new(&reference_json_path),
                    strand_filter,
                );

            // Generate error-checked vectors of seqs and names for the debruijn index
            let (reference_seqs_thread, reference_seqs_rev_thread, reference_names_thread) =
                utils::validate_reference_pairs(&reference_metadata_thread);

            // Create debruijn index of the reference library
            let reference_index_forward =
                debruijn_mapping::build_index::build_index::<debruijn::kmer::Kmer30>(
                    &reference_seqs_thread,
                    &reference_names_thread,
                    &HashMap::new(),
                    num_cores,
                )
                .expect("Error -- could not create pseudoaligner index of the reference library");

            // Create debruijn index of the reverse-comp of the reference library
            let reference_index_reverse = debruijn_mapping::build_index::build_index::<
                debruijn::kmer::Kmer30,
            >(
                &reference_seqs_rev_thread,
                &reference_names_thread,
                &HashMap::new(),
                num_cores,
            )
            .expect(
                "Error -- could not create reverse pseudoaligner index of the reference library",
            );

            // Update the shared data
            let mut reference_indices_lock = reference_indices.lock().unwrap();
            reference_indices_lock.push((reference_index_forward, reference_index_reverse));
            drop(reference_indices_lock);

            let mut reference_metadata_lock = reference_metadata.lock().unwrap();
            reference_metadata_lock.push(reference_metadata_thread);
            drop(reference_metadata_lock);

            let mut align_config_lock = align_config.lock().unwrap();
            align_config_lock.push(align_config_thread);
            drop(align_config_lock);

            let mut reference_seqs_lock = reference_seqs.lock().unwrap();
            reference_seqs_lock.push(
                reference_seqs_thread
                    .into_iter()
                    .map(|s| s.to_string())
                    .collect(),
            );
            drop(reference_seqs_lock);

            let mut reference_names_lock = reference_names.lock().unwrap();
            reference_names_lock.push(reference_names_thread);
            drop(reference_names_lock);
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    println!("Loading read sequences and aligning");

    let reference_indices = Arc::try_unwrap(reference_indices)
        .unwrap()
        .into_inner()
        .unwrap();
    let reference_metadata = Arc::try_unwrap(reference_metadata)
        .unwrap()
        .into_inner()
        .unwrap();
    let align_config = Arc::try_unwrap(align_config).unwrap().into_inner().unwrap();

    bam::process(
        input_files,
        reference_indices,
        reference_metadata,
        align_config,
        output_paths,
        num_cores,
    );

    println!("Alignment successful, terminating.")
}
