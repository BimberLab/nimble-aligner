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

pub fn get_data(
    seq_filename: &str,
    lib_filename: &str,
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

    library.push("libraries/");
    library.push(lib_filename);

    let mut sequences = data_path.clone();
    sequences.push("reads/");
    sequences.push(seq_filename);

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

pub fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}
