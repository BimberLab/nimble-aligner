use std::env;
use std::path;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::pairwise::banded::*;
use bio::alignment::sparse::hash_kmers;
use bio::alignment::pairwise::Scoring;

fn main() {
  const ALIGN_SCORE_THRESHOLD: usize = 60;

  // TODO: Make this safely parse arguments
  let args: Vec<String> = env::args().collect();

  // Read the reference genome file 
  // TODO: Replace unwraps() with error handling where appropriate
  let reference_genome = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap();

  // Read the paired-end sequences that will be aligned
  // TODO: Handle single-end sequences
  let sequence = fastq::Reader::from_file(path::Path::new(&args[2])).unwrap();
  let reverse_sequence = fastq::Reader::from_file(path::Path::new(&args[3])).unwrap();

  /* TODO: We will want to make the choice of global vs local vs semiglobal alignment depend on the reference library
     For now, we're running the local banded alignment. */
  const K: usize = 6; // kmer match length
  const W: usize = 20; // Window size for creating the band
  const MATCH: i32 = 1; // Match score
  const MISMATCH: i32 = -1; // Mismatch score
  const GAP_OPEN: i32 = -3; // Gap open score
  const GAP_EXTEND: i32 = -1; // Gap extend score
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
  let mut aligner = Aligner::with_scoring(scoring, K, W);;

  // Get sequences from the input files and align them
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Parallelize alignment loop, should decrease alignment runtime proportional to the number of cores
  let reference_genome = reference_genome.records();
  let mut sequence = sequence.records();
  let mut reverse_sequence = reverse_sequence.records();

  let reference_genome_size = 1210;//&reference_genome.size_hint().1;
  let mut current_reference_num = 1;

  // Get top 1000 of iterators and put them in vectors
  // TODO: Remove these and use the iterators directly instead
  let mut sequence_vec = Vec::new();
  let mut reverse_sequence_vec = Vec::new();
  for _ in 1..100 {
    sequence_vec.push(sequence.next().unwrap());
  }

  for _ in 1..100 {
    reverse_sequence_vec.push(reverse_sequence.next().unwrap());
  }

  for reference in reference_genome {
    let reference = reference.unwrap();
    let reference_kmers_hash = hash_kmers(reference.seq(), K);
    
    for record in &sequence_vec {
      if let Ok(seq) = record {
        let score = aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        //let score = aligner.local(seq.seq(), reference.seq()).score;
        //let score = aligner.semiglobal_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        print!("Sequence score: {}\n", score)
      }
    }

    for reverse_record in &reverse_sequence_vec {
      if let Ok(seq) = reverse_record {
        let score = aligner.custom_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        //let score = aligner.local(seq.seq(), reference.seq()).score;
        //let score = aligner.semiglobal_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score;
        print!("Reverse sequence score: {}\n", score)
      }
    }

    print!("Reference genome state: {:?}/{:?}\n", current_reference_num, reference_genome_size);
    current_reference_num += 1;
  }
}
