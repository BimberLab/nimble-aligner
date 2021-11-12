extern crate csv;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate nimble;

use std::collections::HashMap;
use std::collections::HashSet; 
use std::path::PathBuf;

use nimble::reference_library;
use nimble::utils;

#[test]
fn bam_pipeline() {
    let seq_filename = "real-data.bam";
    let lib_filename = "real-data.json";

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
        utils::get_valid_reference_sequence_lists(&reference_metadata);

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

    let (results, _) = nimble::process::bam::process(
        vec![sequences.to_str().unwrap()],
        &reference_index,
        &reference_metadata,
        &align_config,
        None
    );
    let results = utils::sort_score_vector(results);
    let mut lengths = HashSet::new();

    for (result, _) in results {
        lengths.insert(result.len());
    };

    let mut expected_lengths = HashSet::new(); 
    expected_lengths.insert(1usize);
    expected_lengths.insert(2usize);
    expected_lengths.insert(3usize);
    expected_lengths.insert(4usize);
    expected_lengths.insert(5usize);
    expected_lengths.insert(6usize);
    expected_lengths.insert(7usize);
    expected_lengths.insert(8usize);
    expected_lengths.insert(9usize);

    assert_eq!(lengths, expected_lengths);
}
