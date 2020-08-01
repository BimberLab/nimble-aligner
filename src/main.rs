mod filter;
mod utils;

use std::env;
use std::path;
use std::fs::File;
use std::io::Write;
use std::collections::HashMap;
use bio::io::fasta;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

fn main() {
  // Program parameters
  // TODO: Get these from the reference library or console
  const MATCH_THRESHOLD: usize = 60;
  const NUM_CORES: usize = 4;

  // Record sizes
  // TODO: Get these from reference library and/or data files
  const REFERENCE_GENOME_SIZE: usize = 1209;
  const READS_SIZE: usize = 2786342;

  // TODO: Make this safely parse arguments
  let args: Vec<String> = env::args().collect();

  // Get iterators to records in the reference genome file and library
  // TODO: Make handle FASTQ.gz rather than having to uncompress manually first
  // TODO: Replace unwraps() with error handling where appropriate
  let reference_genome = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap().records();
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(&args[4])).unwrap());

  // Load reference genome and library into memory
  let mut reference_seqs = Vec::new();
  for reference in reference_genome {
    let reference = reference.unwrap();
    reference_seqs.push(DnaString::from_acgt_bytes(reference.seq())); // Convert raw data to DNAString
  }
  
  let mut reference_names = Vec::new();
  for reference in reference_library.into_records() {
    let reference = reference.unwrap();
    reference_names.push(reference[0].to_string());
  }

  // Create debruijn-mapped index of the reference library
  // TODO: Check ways to make NUM_CORES more effective
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    NUM_CORES
  ).unwrap();

  // Get iterator to the sequences that will be aligned to the reference from the sequence genome file(s)
  // TODO: Handle single-end sequences, currently only does paired-end
  let sequences = fastq::Reader::from_file(path::Path::new(&args[2])).unwrap().records().map(|record| DnaString::from_acgt_bytes(record.unwrap().seq()));
  let reverse_sequences = fastq::Reader::from_file(path::Path::new(&args[3])).unwrap().records().map(|record| DnaString::from_acgt_bytes(record.unwrap().seq()));
  let mut sequences = sequences.zip(reverse_sequences);

  // Subset values for development
  // TODO: Remove/parameterize this on the console
  let mut subset_values = Vec::new();
  for _ in 0..READS_SIZE {
    subset_values.push(sequences.next().unwrap());
  }
  let sequences = subset_values.iter();

  // Map over the sequence iterator, producing an iterator of score vectors for each sequence
  let reference_scores = sequences.fold(vec![0.0; REFERENCE_GENOME_SIZE], 
    |mut acc, (sequence, reverse_sequence)| {

      /* Generate score and equivalence class for this read by aligning the sequence and reverse sequence against
       * the current reference. This alignment returns any scores that are greater than the match threshold. */
      let seq_score = pseduoalign(&sequence, &reference_index, MATCH_THRESHOLD);
      let rev_seq_score = pseduoalign(&reverse_sequence, &reference_index, MATCH_THRESHOLD);

      // Get the greater of the two scores and the associated equivalence class
      let mut max_score = 0;
      let mut max_eqv_class = Vec::new();
      if let Some((eqv_class, score)) = seq_score {
        max_score = score;
        max_eqv_class = eqv_class;
      } 

      if let Some((eqv_class, score)) = rev_seq_score {
        if score > max_score {
          max_eqv_class = eqv_class;
        }
      }

      // If there was a match, update the results accordingly
      if !max_eqv_class.is_empty() {
        for idx in max_eqv_class {
          acc[idx as usize] += 1.0;
        }
      }

      acc
    }
  );

  // Normalize the results vector down to a percentage
  let reference_scores: Vec<f32> = reference_scores.iter().map(|score| (score / READS_SIZE as f32) * 100.0).collect();

  /* Create reference library iterator and filter results by lineage.
   * 
   * TODO: Make whether or not to do this configurable in the library.
   * 
   * TODO: This just works based on the grouping string e.g. lineage. There should be some handling
   * after this to heuristically generate names for groups without a name.
   */
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(&args[4])).unwrap());
  let results = filter::collapse_results_by_lineage(reference_library.into_records(), reference_scores.iter());

  // Write filtered results to file
  let mut str_rep = String::new();

  str_rep += "lineage";
  str_rep += "\t";
  str_rep += "percent from locus\n";

  for (group, score) in results {
    if score > 0.0 {
      str_rep += &group.to_string();
      str_rep += "\t";
      str_rep += &score.to_string();
      str_rep += "\n";
    }
  }

  let mut file = File::create("results.tsv").unwrap();
  file.write_all(str_rep.as_bytes()).unwrap();
}

// Align the given sequence against the given reference with a score threshold
fn pseduoalign(sequence: &DnaString, reference_index: &PseudoAligner, match_threshold: usize) -> Option<(Vec<u32>, usize)> {
  match reference_index.map_read(&sequence) {
    Some((equiv_class, score)) => if score >= match_threshold && !equiv_class.is_empty() { Some((equiv_class, score)) } else { None }
    None => None
  }
}