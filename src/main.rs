use std::env;
use std::path;
use bio::io::fasta;
use bio::io::fastq;
use bio::alignment::pairwise::*;

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
     For now, we're running the local alignment. */
  let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
  let mut aligner = Aligner::new(-5, -1, &score);

  // Get sequences from the input files and align them
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Parallelize alignment loop, should decrease alignment runtime proportional to the number of cores
  let reference_genome = reference_genome.records();
  let mut sequence = sequence.records();
  let mut reverse_sequence = reverse_sequence.records();

  let reference_genome_size = 1210;//&reference_genome.size_hint().1;
  let mut current_reference_num = 1;

  for reference in reference_genome {
    let reference = reference.unwrap();
    let record = sequence.next();
    let reverse_record = reverse_sequence.next();

    print!("Reference genome state: {:?}/{:?}\n", current_reference_num, reference_genome_size);
    if let Some(seq) = record {
      print!("Sequence score: {}\n", aligner.local(seq.unwrap().seq(), reference.seq()).score)
    }

    if let Some(seq) = reverse_record {
      print!("Reverse sequence score: {}\n\n", aligner.local(seq.unwrap().seq(), reference.seq()).score)
    }

    current_reference_num += 1;
  }
}
