extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

#[path = "./utils.rs"]
mod utils;

#[test]
// Case with zero mismatches
fn zero_mismatch() {
    let seq_filename = "mismatch.fastq";
    let lib_filename = "mismatch.json";
    let (sequences, reference_index, reference_metadata, align_config) = utils::get_data(
        seq_filename,
        lib_filename,
        nimble::align::LibraryChemistry::None,
    );

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![(vec![String::from("gene")], (1, Vec::new(), Vec::new()))];

    assert_eq!(results, expected_results);
}
#[test]
// Case with one mismatches
fn one_mismatch() {
    let seq_filename = "mismatch.fastq";
    let lib_filename = "mismatch.json";
    let (sequences, reference_index, reference_metadata, mut align_config) = utils::get_data(
        seq_filename,
        lib_filename,
        nimble::align::LibraryChemistry::None,
    );

    align_config.num_mismatches = 1;

    let (results, _, _) = nimble::align::get_calls(
        sequences,
        None,
        &Vec::new(),
        &reference_index,
        &reference_metadata,
        &align_config,
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![(vec![String::from("gene")], (2, Vec::new(), Vec::new()))];

    assert_eq!(results, expected_results);
}