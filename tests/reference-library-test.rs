#![feature(allocator_api)]
extern crate nimble;

use std::path::PathBuf;

#[test]
fn reference_library_test() {
    let lib_filename = "reference-library-test.json";

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut library = data_path.clone();

    library.push("libraries/");
    library.push(lib_filename);

    let (result_config, result_data) = nimble::reference_library::get_reference_library(
        library.as_path()
    );

    let expected_config = nimble::align::AlignFilterConfig {
        reference_genome_size: 5,
        score_threshold: 60,
        num_mismatches: 0,
        discard_nonzero_mismatch: false,
        discard_multiple_matches: false,
        score_filter: 25,
        intersect_level: nimble::align::IntersectLevel::NoIntersect,
        require_valid_pair: false,
        discard_multi_hits: 0,
    };

    let expected_data = nimble::reference_library::ReferenceData {
        group_on: 1,
        headers: vec!["reference_genome", "sequence_name", "nt_length", "sequence"].into_iter().map(String::from).collect(),
        columns: vec![
            vec!["basic", "basic", "basic", "basic", "basic"].into_iter().map(String::from).collect(),
            vec!["A02-0", "A02-1", "A02-2", "A02-LC", "KIR2DL-4"].into_iter().map(String::from).collect(),
            vec!["180", "180", "180", "180", "180"].into_iter().map(String::from).collect(),
            vec!["A", "A", "A", "A", "A"].into_iter().map(String::from).collect()
        ],
        sequence_name_idx: 1,
        sequence_idx: 3,
        data_type: String::from("DNA")
    };

    assert_eq!(result_config, expected_config);
    assert_eq!(result_data, expected_data);
}

#[test]
#[should_panic(expected = "Error -- could not parse discard_multiple_mismatches as boolean")]
fn get_reference_library_fail_discard_multiple_matches() {
    let lib_filename = "discard_multiple_matches_error.json";

    let mut library = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    library.push("tests/test-sequences/libraries/parse-error-libs");
    library.push(lib_filename);

    nimble::reference_library::get_reference_library(library.as_path());
}

#[test]
#[should_panic(expected = "Error -- could not parse group_on as string")]
fn get_reference_library_fail_group_on() {
    let lib_filename = "group_on_error.json";

    let mut library = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    library.push("tests/test-sequences/libraries/parse-error-libs");
    library.push(lib_filename);

    nimble::reference_library::get_reference_library(library.as_path());
}

#[test]
#[should_panic(expected = "Error -- could not parse score_threshold as int64")]
fn get_reference_library_fail_score_threshold() {
    let lib_filename = "score_threshold_error.json";

    let mut library = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    library.push("tests/test-sequences/libraries/parse-error-libs");
    library.push(lib_filename);

    nimble::reference_library::get_reference_library(library.as_path());
}