use std::env;
use std::path;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::pairwise::banded::*;
use bio::alignment::sparse::hash_kmers;

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
  let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
  let k = 8;  // kmer match length
  let w = 6;  // Window size for creating the band
  let mut aligner = Aligner::new(-5, -1, score, k, w);

  // Get sequences from the input files and align them
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Parallelize alignment loop, should decrease alignment runtime proportional to the number of cores
  let reference_genome = reference_genome.records();
  let mut sequence = sequence.records();
  let mut reverse_sequence = reverse_sequence.records();

  let reference_genome_size = 1210;//&reference_genome.size_hint().1;
  let mut current_reference_num = 1;
  let mut curr_seq_num = 1;

  for reference in reference_genome {
    let reference = reference.unwrap();
    //let reference_kmers_hash = hash_kmers(reference.seq(), k);
    
    for record in sequence.by_ref() {
      if let Ok(seq) = record {
        //print!("Sequence score: {}\n", aligner.semiglobal_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score)
        print!("Sequence score: {}\n", aligner.local(seq.seq(), reference.seq()).score)
      }
      print!("{}\n", curr_seq_num);
      curr_seq_num += 1;
    }

    curr_seq_num = 1;

    for reverse_record in reverse_sequence.by_ref() {
      if let Ok(seq) = reverse_record {
        //print!("Reverse sequence score: {}\n", aligner.semiglobal_with_prehash(seq.seq(), reference.seq(), &reference_kmers_hash).score)
        print!("Reverse sequence score: {}\n", aligner.local(seq.seq(), reference.seq()).score)
      }
      print!("{}\n", curr_seq_num);
      curr_seq_num += 1;
    }

    print!("Reference genome state: {:?}/{:?}\n", current_reference_num, reference_genome_size);
    current_reference_num += 1;
  }
}
