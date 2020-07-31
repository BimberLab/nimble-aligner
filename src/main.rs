mod filter;
mod utils;

use std::env;
use std::path;
use std::fs::File;
use std::io::prelude::*;
//use rayon::prelude::*;
//use rayon::iter::ParallelBridge;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::sparse::hash_kmers;
use bio::alignment::pairwise::Scoring;
use bio::alignment::pairwise::banded::*;

fn main() {

  // TODO: Configure these parameters in the reference library
  // Aligner parameters
  const K: usize = 6; // kmer match length, used to configure aligners and kmer hashing function
  const W: usize = 20; // Window size for creating the band
  const GAP_OPEN: i32 = -3;
  const GAP_EXTEND: i32 = -1;

  // Filter parameters
  const MATCH_THRESHOLD: i32 = 60; // Threshold to determine whether an alignment is a match
  const REPORT_THRESHOLD: f32 = 0.004; // Number of matches below this value is discarded

  // TODO: Get these from reference library and/or data files
  const REFERENCE_GENOME_SIZE: usize = 1209;
  const READS_SIZE: usize = 1;

  // TODO: Make this safely parse arguments
  let args: Vec<String> = env::args().collect();

  // Get iterator to the sequences that will be aligned to the reference from the sequence genome file(s)
  // TODO: Handle single-end sequences, currently only does paired-end
  let sequences = fastq::Reader::from_file(path::Path::new(&args[2])).unwrap().records();
  let reverse_sequences = fastq::Reader::from_file(path::Path::new(&args[3])).unwrap().records();
  let mut sequences = sequences.zip(reverse_sequences);

  // Subset values for development
  // TODO: Remove/parameterize this on the console
  let mut subset_values = Vec::new();
  for _ in 0..READS_SIZE {
    subset_values.push(sequences.next().unwrap());
  }
  let sequences = subset_values.iter();

  // Parallelized map over the sequence iterator, producing an iterator of score vectors for each sequence
  // TODO: Configure parallelism from console arguments
  let reference_scores = sequences./*par_bridge().*/fold(vec![0.0; REFERENCE_GENOME_SIZE], 
    |mut acc, (sequence, reverse_sequence)| {
      // Get iterator to records in the reference genome file 
      // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
      // TODO: Replace unwraps() with error handling where appropriate
      let reference_genome = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap().records();

      // Create banded aligner
      let mut aligner = Aligner::with_scoring(Scoring {
        gap_open: GAP_OPEN,
        gap_extend: GAP_EXTEND,
        match_fn: |a: u8, b: u8| if a == b { 1 } else { -1 },
        match_scores: Some((1, -1)),
        xclip_prefix: 0,
        xclip_suffix: 0,
        yclip_prefix: 0,
        yclip_suffix: 0
      }, K, W);

      // Generate score vector for this read
      let sequence_scores = reference_genome.map(|reference| {
        // Get the reference and its hash
        let reference = reference.unwrap();
        let reference_kmers_hash = hash_kmers(reference.seq(), K);

        // Align the sequence and reverse sequence against the current reference
        let sequence_score = match sequence {
          Ok(seq) => aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score,
          Err(_) => 0
        };

        let reverse_sequence_score = match reverse_sequence {
          Ok(seq) => aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score,
          Err(_) => 0
        };

        // Return the highest score
        if sequence_score > reverse_sequence_score { sequence_score } else { reverse_sequence_score }
      });

      // If there's a score above the match threshold, get all matches >= to that score and mutate the acc vector
      // TODO: Attach name to reference so we can sort in O(ln) instead of iterating in O(n)
      let sequence_scores: Vec<i32> = sequence_scores.collect();
      let max_score = sequence_scores.iter().max().unwrap();
      if max_score >= &MATCH_THRESHOLD {
        for (i, v) in sequence_scores.iter().enumerate() {
          if v >= max_score {
            acc[i] += 1.0
          }
        }
      }

      acc
    }
  );

  let reference_scores: Vec<f32> = reference_scores.iter().map(|score| score / READS_SIZE as f32).collect();

  /* Create reference library iterator and filter results by lineage.
   * 
   * TODO: Make whether or not to do this configurable in the library.
   * 
   * TODO: This just works based on the grouping string e.g. lineage. There should be some handling
   * after this to heuristically generate names for groups without a name.
   */
  // Get iterator to the reference library
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(&args[4])).unwrap());
  let results = filter::collapse_results_by_lineage(reference_library.into_records(), reference_scores.iter());

  // Write filtered results to file
  let mut str_rep = String::new();

  str_rep += "lineage";
  str_rep += "\t";
  str_rep += "percent from locus\n";

  for (group, score) in results {
    if score > REPORT_THRESHOLD {
      str_rep += &group.to_string();
      str_rep += "\t";
      str_rep += &score.to_string();
      str_rep += "\n";
    }
  }

  let mut file = File::create("results.tsv").unwrap();
  file.write_all(str_rep.as_bytes()).unwrap();
}