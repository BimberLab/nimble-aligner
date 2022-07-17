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

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch_forward() {
    let seq_filename = "basic.fastq";
    let lib_filename = "basic.json";
    let (sequences, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename);

    let (results, _) = nimble::align::score(
        sequences,
        None,
        &reference_index,
        &reference_metadata,
        &align_config,
        None
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