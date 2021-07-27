#![feature(allocator_api)]
extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

use nimble::align;
use nimble::parse;
use nimble::reference_library;
use nimble::utils;
use std::collections::HashMap;
use std::path::PathBuf;

use debruijn::dna_string::DnaString;
use std::io::Error;

// Shared function for generating basic single strand test data
fn get_basic_single_strand_data(
    reverse_comp_ref: bool,
) -> (
    (
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
    ),
    (align::PseudoAligner, align::PseudoAligner),
    reference_library::ReferenceMetadata,
    align::AlignFilterConfig,
) {
    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut library = data_path.clone();

    if !reverse_comp_ref {
        library.push("libraries/basic-rev.json");
    } else {
        library.push("libraries/basic.json");
    }

    let mut sequences = data_path.clone();
    sequences.push("reads/basic.fastq");

    let (align_config, reference_metadata) =
        reference_library::get_reference_library(library.as_path());

    let (reference_seqs, reference_seqs_rev, reference_names) =
        utils::validate_reference_pairs(&reference_metadata);

    let reference_index_forward = debruijn_mapping::build_index::build_index::<
        debruijn_mapping::config::KmerType,
    >(&reference_seqs, &reference_names, &HashMap::new(), 1)
    .expect("Error -- could not create pseudoaligner index of the unit test reference library");

    let reference_index_reverse = debruijn_mapping::build_index::build_index::<
        debruijn_mapping::config::KmerType,
    >(&reference_seqs_rev, &reference_names, &HashMap::new(), 1)
    .expect(
        "Error -- could not create reverse pseudoaligner index of the unit test reference library",
    );

    let reference_index = (reference_index_forward, reference_index_reverse);

    let sequences = parse::fastq::get_error_checked_fastq_readers(
        &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice."),
    );

    (sequences, reference_index, reference_metadata, align_config)
}

fn get_group_by_data(
    reverse_comp_ref: bool,
) -> (
    (
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
    ),
    (align::PseudoAligner, align::PseudoAligner),
    reference_library::ReferenceMetadata,
    align::AlignFilterConfig,
) {
    let (sequences, reference_index, mut reference_metadata, align_config) =
        get_basic_single_strand_data(reverse_comp_ref);

    reference_metadata.group_on = 4;
    reference_metadata.headers.push("test_group_on".to_string());
    reference_metadata.columns.push(
        vec!["g1", "g2", "g2", "g1", "g1"]
            .into_iter()
            .map(|column| column.to_string())
            .collect(),
    );

    (sequences, reference_index, reference_metadata, align_config)
}

fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_forward() {
    let (sequences, reference_index, reference_metadata, align_config) =
        get_basic_single_strand_data(false);

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with one mismatch
fn basic_single_strand_one_mismatch_forward() {
    let (sequences, reference_index, reference_metadata, mut align_config) =
        get_basic_single_strand_data(false);

    align_config.num_mismatches = 1;

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with two mismatches
fn basic_single_strand_two_mismatch_forward() {
    let (sequences, reference_index, reference_metadata, mut align_config) =
        get_basic_single_strand_data(false);

    align_config.num_mismatches = 2;

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_reverse() {
    let (sequences, reference_index, reference_metadata, align_config) =
        get_basic_single_strand_data(true);

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with one mismatch
fn basic_single_strand_one_mismatch_reverse() {
    let (sequences, reference_index, reference_metadata, mut align_config) =
        get_basic_single_strand_data(true);

    align_config.num_mismatches = 1;

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with two mismatches
fn basic_single_strand_two_mismatch_reverse() {
    let (sequences, reference_index, reference_metadata, mut align_config) =
        get_basic_single_strand_data(true);

    align_config.num_mismatches = 2;

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            1,
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], 1),
        (vec![String::from("A02-1")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with group_by instead of basic forward allele-level reporting
fn group_by_forward() {
    let (sequences, reference_index, reference_metadata, align_config) = get_group_by_data(false);

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (vec![String::from("g1")], 1),
        (vec![String::from("g1"), String::from("g2")], 1),
        (vec![String::from("g2")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with group_by instead of basic reverse allele-level reporting
fn group_by_rev() {
    let (sequences, reference_index, reference_metadata, align_config) = get_group_by_data(true);

    let results = nimble::align::score(
        sequences,
        None,
        reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = sort_score_vector(results);

    let expected_results = vec![
        (vec![String::from("g1")], 1),
        (vec![String::from("g1"), String::from("g2")], 1),
        (vec![String::from("g2")], 2),
    ];
    let expected_results = sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}
