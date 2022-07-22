use std::path::PathBuf;

extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

#[path = "./utils.rs"]
mod utils;

#[test]
fn strandedness() {
    let seq_filename = "test_unstranded_filter.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::Unstranded);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");

    let results = nimble::process::bam::process(
        vec![sequences],
        &reference_index,
        &reference_metadata,
        &align_config,
        "/dev/null",
        None,
        None
    );

    let results = utils::sort_score_vector(results);

    let expected_results = vec![
        (
            vec![
                String::from("mismatch_test"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}