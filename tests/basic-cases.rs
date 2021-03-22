#![feature(min_type_alias_impl_trait)]
extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

use bio::alphabets::rna::revcomp;
use debruijn::dna_string::DnaString;
use nimble::align;
use nimble::align::IntersectLevel;
use nimble::reference_library;
use nimble::utils::{validate_reference_pairs, DNAStringIter};
use std::collections::HashMap;
use std::io::Error;

fn get_raw_data(reverse_comp_ref: bool) -> (Vec<String>, Vec<String>) {
    let reference_sequences = vec![
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCcTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGgAAGGAGACGCTgCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCgTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CACTCCCCCACTGAGTGGTCGGCACCCAGCAACCCCCTGGTGATCATGGTCACAGGTCTATATGAGAAACCTTCTCTCTCAGCCCAGCCGGGCCCCACGGTTCCCACAGGAGAGAACATGACCTTGTCCTGCAGTTCCCGGCGCTCCTTTGACATGTACCATCTATCCAGGGAGGGGGAG"];

    // Test sequences
    let sequences = vec![
    "TACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCT",                  // Test-Data-1: exact match to A02-0
    "TACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCT",                  // Test-Data-2: exact match to A02-1
    "TACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTatgatgatagatag",    // Test-Data-3: exact match to A02-1, except has extraneous bases at end
    "CAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGC"                  // Test-Data-4: exact match to 5' end of A02-0 through A02-2
  ];

    let sequences = sequences.into_iter().map(String::from).collect();

    let reference_sequences = if reverse_comp_ref {
        reference_sequences
            .into_iter()
            .map(String::from)
            .map(|seq| seq.into_bytes())
            .map(revcomp)
            .map(|seq| String::from_utf8(seq).unwrap())
            .collect()
    } else {
        reference_sequences.into_iter().map(String::from).collect()
    };

    (reference_sequences, sequences)
}

// Shared function for generating basic single strand test data
fn get_basic_single_strand_data(
    reverse_comp_ref: bool,
) -> (
    (DNAStringIter, DNAStringIter),
    (align::PseudoAligner, align::PseudoAligner),
    reference_library::ReferenceMetadata,
) {
    let (reference_sequences, sequences) = get_raw_data(reverse_comp_ref);

    /* 'A02' is a portion of a the macaque MHC sequence Mamu-A1*002. A02-1 and A02-2 are 1bp and 2bp changes form A02 (see lower-case bases).
    A02-LC is the same sequence as A02, just with some upper -> lower case changes to ensure that our results are case-insensitive. */
    let mut columns: Vec<Vec<String>> = vec![
        vec!["test", "test", "test"],
        vec!["A02-0", "A02-1", "A02-2", "A02-LC", "KIR2DL-4"],
        vec!["180", "180", "180"],
    ]
    .into_iter()
    .map(|column| column.into_iter().map(|val| val.to_string()).collect())
    .collect();
    columns.push(reference_sequences);

    let sequences: Vec<Result<DnaString, Error>> = sequences
        .into_iter()
        .map(|seq| Ok(DnaString::from_dna_string(&seq)))
        .collect();

    let reference_metadata = reference_library::ReferenceMetadata {
        group_on: 1,
        headers: vec![
            "reference_genome",
            "sequence_name_idx",
            "nt_length",
            "sequence_idx",
        ]
        .into_iter()
        .map(|header| header.to_string())
        .collect(),
        columns,
        sequence_name_idx: 1,
        sequence_idx: 3,
    };

    let (reference_seqs, reference_seqs_rev, reference_names) =
        validate_reference_pairs(&reference_metadata);

    let reference_index_forward = debruijn_mapping::build_index::build_index::<
        debruijn_mapping::config::KmerType,
    >(&reference_seqs, &reference_names, &HashMap::new(), 1)
    .expect("Error -- could not create pseudoaligner index of the reference library");

    let reference_index_reverse = debruijn_mapping::build_index::build_index::<
        debruijn_mapping::config::KmerType,
    >(&reference_seqs_rev, &reference_names, &HashMap::new(), 1)
    .expect("Error -- could not create pseudoaligner index of the reference library");

    let reference_index = (reference_index_forward, reference_index_reverse);

    (cast_to_iter_pair(sequences), reference_index, reference_metadata)
}

fn cast_to_iter_pair(sequences: Vec<Result<DnaString, Error>>) -> (DNAStringIter, DNAStringIter) {
  (sequences.into_iter(), sequences.into_iter())
}

fn get_group_by_data() -> (
    (DNAStringIter, DNAStringIter),
    (align::PseudoAligner, align::PseudoAligner),
    reference_library::ReferenceMetadata,
) {
    let (sequences, reference_index, mut reference_metadata) = get_basic_single_strand_data(false);

    reference_metadata.group_on = 4;
    reference_metadata.headers.push("test_group_on".to_string());
    reference_metadata.columns.push(
        vec!["g1", "g2", "g2", "g1", "g1"]
            .into_iter()
            .map(|column| column.to_string())
            .collect(),
    );

    (sequences, reference_index, reference_metadata)
}

fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_forward() {
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(false);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 0,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(false);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 1,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(false);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 2,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(true);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 0,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(true);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 1,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
    let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data(true);

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 2,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
// Case with group_by instead of basic allele-level reporting
fn group_by() {
    let (sequences, reference_index, reference_metadata) = get_group_by_data();

    // Configure aligner
    let align_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 0,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 0,
        require_valid_pair: false,
        intersect_level: IntersectLevel::NoIntersect,
        discard_multi_hits: 0,
    };

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
