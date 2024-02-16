/*extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

#[path = "./utils.rs"]
mod utils;

#[test]
// Case with zero mismatches
fn mismatch() {
    let seq_filename = "mismatch.fastq";
    let lib_filename = "mismatch.json";
    let (sequences, reference_index, reference_metadata, align_config) = utils::get_data(
        seq_filename,
        lib_filename,
        nimble::align::StrandFilter::None,
    );

    let (results, _) = nimble::align::score(
        sequences,
        None,
        &reference_index,
        &reference_metadata,
        &align_config,
        None,
    );
    let results = utils::sort_score_vector(results);

    let expected_results = vec![(vec![String::from("NKG2E_NM_001104593")], 2)];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}*/
