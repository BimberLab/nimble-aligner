use std::env;
use std::path;
use std::fs::File;
use std::io::prelude::*;
use rayon::prelude::*;
use rayon::iter::ParallelBridge;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::sparse::hash_kmers;
use bio::alignment::pairwise::{Scoring,MatchFunc};
use bio::alignment::pairwise::banded::*;

fn main() {
  const K: usize = 6; // kmer match length, used to configure aligners and kmer hashing function
  const MATCH_THRESHOLD: i32 = 60; // Threshold to determine whether an alignment is a match

  // TODO: Make this safely parse arguments
  let args: Vec<String> = env::args().collect();

  // Get iterator to records in the reference genome file 
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Replace unwraps() with error handling where appropriate
  let reference_genome = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap().records();

  // Parallelized map over the reference library iterator, producing an iterator of score vectors for each reference
  let scores = reference_genome.par_bridge().map(|reference| { 
    // Get the reference and its hash
    let reference = reference.unwrap();
    let reference_kmers_hash = hash_kmers(reference.seq(), K);

    // Get iterators to the sequences that will be aligned to the reference from the sequence genome file(s)
    // TODO: Handle single-end sequences, currently only does paired-end
    let sequences = fastq::Reader::from_file(path::Path::new(&args[2])).unwrap().records();
    let reverse_sequences = fastq::Reader::from_file(path::Path::new(&args[3])).unwrap().records();

    // Create banded aligners for the sequence and the reverse sequence
    let match_fn = |a: u8, b: u8| if a == b { 1 } else { -1 };
    let mut sequence_aligner = get_mutable_aligner(K, match_fn);
    let mut reverse_sequence_aligner = get_mutable_aligner(K, match_fn);

    // Align the sequence and reverse sequence against the current reference and return an iterator to the scores
    let sequence_scores = sequences.map(|record| {
      match record {
        Ok(seq) => sequence_aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score,
        Err(_) => 0
      }
    });

    let reverse_sequence_scores = reverse_sequences.map(|record| {
      match record {
        Ok(seq) => reverse_sequence_aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score,
        Err(_) => 0
      }
    });

    // Concatenate the iterators together
    let mut scores = sequence_scores.chain(reverse_sequence_scores);

    // Subset the iterator for development
    // TODO: remove this subsetting
    let mut subset_values = Vec::new();
    for _ in 1..2000 {
      subset_values.push(scores.next().unwrap());
    }

    subset_values
  });

  // Run the iterator and collect the values into memory
  // TODO: Consume lazily so that we're not loading the whole matrix into memory
  let score_matrix: Vec<Vec<i32>> = scores.collect();

  // Fold the score matrix into a first-pass results matrix
  let results: Vec<f64> = score_matrix.iter().map(|score_vec| {
    let score_vec_len = score_vec.len();
    let ratio = score_vec.iter().fold(0.0, |acc, score| if score >= &MATCH_THRESHOLD { acc+1.0 } else { acc }) / score_vec_len as f64;
    ratio
  }).collect();

  print!("{:?}\n", results);
}

// Creates a pre-configured aligner
fn get_mutable_aligner<F: MatchFunc>(k: usize, function: F) -> Aligner<F> {
    // Define aligner scoring parameters
    const W: usize = 20; // Window size for creating the band
    const GAP_OPEN: i32 = -3;
    const GAP_EXTEND: i32 = -1;

    let scoring = Scoring {
      gap_open: GAP_OPEN,
      gap_extend: GAP_EXTEND,
      match_fn: function,
      match_scores: Some((1, -1)),
      xclip_prefix: 0,
      xclip_suffix: 0,
      yclip_prefix: 0,
      yclip_suffix: 0,
    };

  /* TODO: We will want to make the choice of global vs local vs semiglobal alignment depend on the reference library
   * For now, we're running the local banded alignment with hashing.
   */
  Aligner::with_scoring(scoring, k, W)
}

// Convert to TSV and write to disc
// TODO: Remove this, it's only needed for testing the aligner
fn _write_alignment_scores_to_tsv(score_matrix: Vec<Vec<i32>>) {
  let mut str_rep = String::new();
  for score_vec in score_matrix {
    for n in score_vec {
      str_rep += &(n.to_string() + "\t")[..];
    }

    str_rep += "\n";
  }

  let mut file = File::create("aligner_score_matrix.tsv").unwrap();
  file.write(str_rep.as_bytes()).unwrap();
}