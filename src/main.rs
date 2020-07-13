use std::env;
use std::path;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::pairwise::banded::*;
use bio::alignment::sparse::hash_kmers;
use bio::alignment::pairwise::Scoring;
use rayon::prelude::*;
use rayon::iter::ParallelBridge;

fn main() {
  // TODO: Make this safely parse arguments
  let args: Vec<String> = env::args().collect();

  // Get iterator to records in the reference genome file 
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Replace unwraps() with error handling where appropriate
  let reference_genome = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap().records();

  // Parallel map over the reference library
  reference_genome.par_bridge().for_each(|reference| { 
    // Get the reference and its hash
    let reference = reference.unwrap();
    let reference_kmers_hash = hash_kmers(reference.seq(), K);

    // Read the paired-end sequences that will be aligned
    // TODO: Handle single-end sequences
    let sequence = fastq::Reader::from_file(path::Path::new(&args[2])).unwrap().records();
    let reverse_sequence = fastq::Reader::from_file(path::Path::new(&args[3])).unwrap().records();

    // Define aligner scoring parameters
    /* TODO: We will want to make the choice of global vs local vs semiglobal alignment depend on the reference library
     * For now, we're running the local banded alignment with hashing, configured with this scoring struct.
     */
    const K: usize = 6; // kmer match length
    const W: usize = 20; // Window size for creating the band
    const MATCH: i32 = 1;
    const MISMATCH: i32 = -1;
    const GAP_OPEN: i32 = -3;
    const GAP_EXTEND: i32 = -1;

    let scoring = Scoring {
      gap_open: GAP_OPEN,
      gap_extend: GAP_EXTEND,
      match_fn: |a: u8, b: u8| if a == b { MATCH } else { MISMATCH },
      match_scores: Some((MATCH, MISMATCH)),
      xclip_prefix: 0,
      xclip_suffix: 0,
      yclip_prefix: 0,
      yclip_suffix: 0,
    };

    let mut aligner = Aligner::with_scoring(scoring, K, W);

    for record in sequence {
      if let Ok(seq) = record {
        let score = aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        print!("{}\n", score);
      }
    }

    for record in reverse_sequence {
      if let Ok(seq) = record {
        let score = aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        print!("{}\n", score);
      }
    }
  });
}
