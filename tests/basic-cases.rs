extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

use nimble::align;
use nimble::reference_library;

use debruijn::dna_string::DnaString;
use std::io::Error;

#[path = "./utils.rs"]
mod utils;

fn get_group_by_data(
    seq_filename: &str,
    lib_filename: &str
) -> (
    (
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
    ),
    align::PseudoAligner,
    reference_library::Reference,
    align::AlignFilterConfig,
) {
    let (sequences, reference_index, mut reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

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


#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_forward() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("A02-0"),
                String::from("A02-1"),
                String::from("A02-2"),
                String::from("A02-LC"),
            ],
            (1, Vec::new(), Vec::new()),
        ),
        (vec![String::from("A02-0"), String::from("A02-LC")], (1, Vec::new(), Vec::new())),
        (vec![String::from("A02-1")], (2, Vec::new(), Vec::new())),
    ];
    //let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}
/*
#[test]
// Case with one mismatch
fn basic_single_strand_one_mismatch_forward() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, mut align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    align_config.num_mismatches = 1;

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

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
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with two mismatches
fn basic_single_strand_two_mismatch_forward() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, mut align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    align_config.num_mismatches = 2;

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

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
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_reverse() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic-rev.json";
    let (sequences, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

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
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with one mismatch
fn basic_single_strand_one_mismatch_reverse() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic-rev.json";
    let (sequences, reference_index, reference_metadata, mut align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    align_config.num_mismatches = 1;

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

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
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with two mismatches
fn basic_single_strand_two_mismatch_reverse() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic-rev.json";
    let (sequences, reference_index, reference_metadata, mut align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::LibraryChemistry::None);

    align_config.num_mismatches = 2;

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

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
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with group_by instead of basic forward allele-level reporting
fn group_by_forward() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, align_config) = get_group_by_data(seq_filename, lib_filename);

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![
        (vec![String::from("g1")], 1),
        (vec![String::from("g1"), String::from("g2")], 1),
        (vec![String::from("g2")], 2),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
// Case with group_by instead of basic reverse allele-level reporting
fn group_by_rev() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, align_config) = get_group_by_data(seq_filename, lib_filename);

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![
        (vec![String::from("g1")], 1),
        (vec![String::from("g1"), String::from("g2")], 1),
        (vec![String::from("g2")], 2),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}
*/