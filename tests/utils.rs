extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

use nimble::align;
use nimble::align::LibraryChemistry;
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
    strand_filter: LibraryChemistry 
) -> (
    (
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
        Box<dyn Iterator<Item = Result<DnaString, Error>>>,
    ),
    align::PseudoAligner,
    reference_library::Reference,
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
        reference_library::get_reference_library(library.as_path(), strand_filter);

    let (reference_seqs, reference_names) =
        utils::get_reference_sequence_data(&reference_metadata);

    let reference_index = debruijn_mapping::build_index::build_index::<
        debruijn::kmer::Kmer30,
    >(&reference_seqs, &reference_names, &HashMap::new(), 1)
    .expect("Error -- could not create pseudoaligner index of the unit test reference library");

    let sequences = parse::fastq::get_error_checked_fastq_readers(
        sequences
            .into_os_string()
            .into_string()
            .expect("Could not convert unit test sequence to OsStr slice."),
    );

    (sequences, reference_index, reference_metadata, align_config)
}

pub fn sort_score_vector(mut scores: Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>) -> Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}