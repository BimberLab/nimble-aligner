extern crate immuno_genotyper;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate csv;

use std::io::Error;
use immuno_genotyper::align::IntersectLevel;
use immuno_genotyper::reference_library;
use std::collections::HashMap;
use debruijn::dna_string::DnaString;
use debruijn_mapping::pseudoaligner::Pseudoaligner;


// Shared function for generating basic single strand test data
fn get_basic_single_strand_data() -> (Vec<Result<DnaString, Error>>, Pseudoaligner<debruijn_mapping::config::KmerType>, reference_library::ReferenceMetadata) {
  /* 'A02' is a portion of a the macaque MHC sequence Mamu-A1*002. A02-1 and A02-2 are 1bp and 2bp changes form A02 (see lower-case bases).  
  A02-LC is the same sequence as A02, just with some upper -> lower case changes to ensure that our results are case-insensitive. */
  let columns: Vec<Vec<String>> = vec![vec!["test", "test", "test"], vec!["A02-0", "A02-1", "A02-2", "A02-LC", "KIR2DL-4"], vec!["180", "180", "180"]].into_iter().map(|column| column.into_iter().map(|val| val.to_string()).collect()).collect();

  // Data associated with the aforementioned names
  let reference_sequences = vec![
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCcTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CGCAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGgAAGGAGACGCTgCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCgTCTCTGACCAAGAGGCCACCCTGAGGTGCTGG",
    "CACTCCCCCACTGAGTGGTCGGCACCCAGCAACCCCCTGGTGATCATGGTCACAGGTCTATATGAGAAACCTTCTCTCTCAGCCCAGCCGGGCCCCACGGTTCCCACAGGAGAGAACATGACCTTGTCCTGCAGTTCCCGGCGCTCCTTTGACATGTACCATCTATCCAGGGAGGGGGAG"];
  let reference_sequences: Vec<DnaString> = reference_sequences.into_iter().map(|seq| DnaString::from_dna_string(seq)).collect();

  // Test sequences
  let sequences = vec![
    "TACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCT",                  // Test-Data-1: exact match to A02-0
    "TACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCT",                  // Test-Data-2: exact match to A02-1
    "TACCTGGAGAACGGGAAGGAGACGCTcCAGCGCGCGGACCCCCCCAAGACACATGTGACCCACCACCCCGTCTCTGACCAAGAGGCCACCCTGAGGTGCTatgatgatagatag",    // Test-Data-3: exact match to A02-1, except has extraneous bases at end
    "CAAGTGGGAGGCGGCGGGTGAGGCGGAGCAGCACAGAACCTACCTGGAGGGCGAGTGCCTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGC"                  // Test-Data-4: exact match to 5' end of A02-0 through A02-2
  ];
  let sequences: Vec<Result<DnaString, Error>> = sequences.into_iter().map(|seq| Ok(DnaString::from_dna_string(seq))).collect();

  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_sequences.as_slice(),
    &columns[1],
    &HashMap::new(),
    1
  ).expect("Error -- could not create pseudoaligner index of the reference library");

  let reference_metadata = reference_library::ReferenceMetadata {
    group_on: 1,
    headers: vec!["reference_genome", "nt_sequence", "nt_length"].into_iter().map(|header| header.to_string()).collect(),
    columns,
    nt_sequence_idx: 1
  };

  (sequences, reference_index, reference_metadata)
}


fn get_group_by_data() -> (Vec<Result<DnaString, Error>>, Pseudoaligner<debruijn_mapping::config::KmerType>, reference_library::ReferenceMetadata) {
  let (sequences, reference_index, mut reference_metadata) = get_basic_single_strand_data();

  reference_metadata.group_on = 3;
  reference_metadata.headers.push("test_group_on".to_string());
  reference_metadata.columns.push(vec!["g1", "g2", "g2", "g1", "g1"].into_iter().map(|column| column.to_string()).collect());

  (sequences, reference_index, reference_metadata)
}

fn sort_score_vector(mut scores: Vec<(String, i32)>) -> Vec<(String, i32)> {
  scores.sort_by(|a, b| a.0.cmp(&b.0));
  scores
}

#[test]
// Case with zero mismatches
fn basic_single_strand_no_mismatch() {
  let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data();

  // Configure aligner
  let align_config = immuno_genotyper::align::AlignFilterConfig {
    reference_genome_size: 5,
    score_threshold: 60,
    num_mismatches: 0,
    discard_nonzero_mismatch: false,
    discard_multiple_matches: false,
    score_filter: 0,
    intersect_level: IntersectLevel::NoIntersect,
    debug_reference: String::new()
  };

  let results = immuno_genotyper::align::score(sequences.into_iter(), None, reference_index, &reference_metadata, &align_config);
  let results = sort_score_vector(results);

  let expected_results = vec![
    (String::from("A02-0"), 2),
    (String::from("A02-1"), 3),
    (String::from("A02-2"), 1),
    (String::from("A02-LC"), 2)];
  let expected_results = sort_score_vector(expected_results);

  assert_eq!(results, expected_results);
}


#[test]
// Case with one mismatch
fn basic_single_strand_one_mismatch() {
  let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data();

  // Configure aligner
  let align_config = immuno_genotyper::align::AlignFilterConfig {
    reference_genome_size: 5,
    score_threshold: 60,
    num_mismatches: 1,
    discard_nonzero_mismatch: false,
    discard_multiple_matches: false,
    score_filter: 0,
    intersect_level: IntersectLevel::NoIntersect,
    debug_reference: String::new()
  };

  let results = immuno_genotyper::align::score(sequences.into_iter(), None, reference_index, &reference_metadata, &align_config);
  let results = sort_score_vector(results);

  let expected_results = vec![
    (String::from("A02-0"), 2),
    (String::from("A02-1"), 3),
    (String::from("A02-2"), 1),
    (String::from("A02-LC"), 2)];
  let expected_results = sort_score_vector(expected_results);

  assert_eq!(results, expected_results);
}


#[test]
// Case with two mismatches
fn basic_single_strand_two_mismatch() {
  let (sequences, reference_index, reference_metadata) = get_basic_single_strand_data();

  // Configure aligner
  let align_config = immuno_genotyper::align::AlignFilterConfig {
    reference_genome_size: 5,
    score_threshold: 60,
    num_mismatches: 2,
    discard_nonzero_mismatch: false,
    discard_multiple_matches: false,
    score_filter: 0,
    intersect_level: IntersectLevel::NoIntersect,
    debug_reference: String::new()
  };

  let results = immuno_genotyper::align::score(sequences.into_iter(), None, reference_index, &reference_metadata, &align_config);
  let results = sort_score_vector(results);

  let expected_results = vec![
    (String::from("A02-0"), 2),
    (String::from("A02-1"), 3),
    (String::from("A02-2"), 1),
    (String::from("A02-LC"), 2)];
  let expected_results = sort_score_vector(expected_results);

  assert_eq!(results, expected_results);
}


#[test]
// Case with group_by instead of basic allele-level reporting
fn group_by() {
  let (sequences, reference_index, reference_metadata) = get_group_by_data();

  // Configure aligner
  let align_config = immuno_genotyper::align::AlignFilterConfig {
    reference_genome_size: 5,
    score_threshold: 60,
    num_mismatches: 0,
    discard_nonzero_mismatch: false,
    discard_multiple_matches: false,
    score_filter: 0,
    intersect_level: IntersectLevel::NoIntersect,
    debug_reference: String::new()
  };

  let results = immuno_genotyper::align::score(sequences.into_iter(), None, reference_index, &reference_metadata, &align_config);
  let results = sort_score_vector(results);
  
  let expected_results = vec![
    (String::from("g1"), 2),
    (String::from("g2"), 3)];
  let expected_results = sort_score_vector(expected_results);

  assert_eq!(results, expected_results);
} 