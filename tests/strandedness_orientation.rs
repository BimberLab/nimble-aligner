use std::path::PathBuf;

extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

#[path = "./utils.rs"]
mod utils;

#[test]
fn unstranded_filter_fr() {
    // SETUP
    let seq_filename = "strandedness_FR.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            2,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_filter_rf() {
    // SETUP
    let seq_filename = "strandedness_RF.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            2,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_filter_fu() {
    // SETUP
    let seq_filename = "strandedness_FU.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_filter_uf() {
    // SETUP
    let seq_filename = "strandedness_UF.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_filter_ru() {
    // SETUP
    let seq_filename = "strandedness_RU.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_filter_ur() {
    // SETUP
    let seq_filename = "strandedness_UR.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn fiveprime_filter_rf_failure() {
    // SETUP
    let seq_filename = "strandedness_RF.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::FivePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}


#[test]
fn fiveprime_filter_ru_failure() {
    // SETUP
    let seq_filename = "strandedness_RU.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::FivePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}


#[test]
fn fiveprime_filter_uf_failure() {
    // SETUP
    let seq_filename = "strandedness_UF.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::FivePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn fiveprime_filter_fr_success() {
    // SETUP
    let seq_filename = "strandedness_FR.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::FivePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            2,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn threeprime_filter_fr_failure() {
    // SETUP
    let seq_filename = "strandedness_FR.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::ThreePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}


#[test]
fn threeprime_filter_fu_failure() {
    // SETUP
    let seq_filename = "strandedness_FU.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::ThreePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}


#[test]
fn threeprime_filter_ur_failure() {
    // SETUP
    let seq_filename = "strandedness_UR.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::ThreePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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

    let expected_results = vec![];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn fiveprime_filter_rf_success() {
    // SETUP
    let seq_filename = "strandedness_RF.sam";
    let lib_filename = "strandedness.json";
    let (_, reference_index, reference_metadata, align_config) = utils::get_data(seq_filename, lib_filename, nimble::align::StrandFilter::ThreePrime);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

    let sequences = &sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice.");


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            2,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}

#[test]
fn unstranded_revcomp_ur_to_uf() {
    // SETUP
    let seq_filename = "UR_reversecomp.sam";
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


    // TEST
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
                String::from("strandedness_test_first"),
            ],
            1,
        ),
    ];
    let expected_results = utils::sort_score_vector(expected_results);

    assert_eq!(results, expected_results);
}