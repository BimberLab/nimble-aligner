mod filter;
mod utils;

use std::path;
use std::fs::File;
use std::io::Write;
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

fn main() {
  // Parse arguments based on the yaml schema
  let yaml = load_yaml!("cli.yml");
  let matches = App::from_yaml(yaml).get_matches();

  let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
  let library = matches.value_of("library").unwrap();
  let library_fasta = matches.value_of("library_fasta").unwrap();
  let num_cores = matches.value_of("num_cores").unwrap_or("1").parse::<usize>().expect("Please provide an integer value for the number of cores.");

  // Parameters for alignment
  const MATCH_THRESHOLD: usize = 60;
  const REFERENCE_GENOME_SIZE: usize = 1209;
  const READS_SIZE: usize = 2786342;

  // Get iterators to records in the reference genome file and library
  let reference_genome = fasta::Reader::from_file(library_fasta).unwrap().records();
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).unwrap());

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
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    num_cores
  ).unwrap();

  // Get iterator to the sequences that will be aligned to the reference from the sequence genome file(s)
  let sequences = fastq::Reader::from_file(path::Path::new(input_files[0])).unwrap().records().map(|record| DnaString::from_acgt_bytes(record.unwrap().seq()));
  let mut reverse_sequences = None;

  if input_files.len() > 1 {
    reverse_sequences = Some(fastq::Reader::from_file(path::Path::new(input_files[1])).unwrap().records().map(|record| DnaString::from_acgt_bytes(record.unwrap().seq())));
  }

  // Map over the sequence iterator, producing an iterator of score vectors for each sequence
  let reference_scores = sequences.fold(vec![0.0; REFERENCE_GENOME_SIZE], 
    |mut acc, read| {

      /* Generate score and equivalence class for this read by aligning the sequence against
       * the current reference. This alignment returns any scores that are greater than the match threshold. */
      let seq_score = pseduoalign(&read, &reference_index, MATCH_THRESHOLD);
      let mut rev_seq_score = None;

      // If there's a reversed sequence, do the paired-end alignment
      if let Some(itr) = &mut reverse_sequences {
        rev_seq_score = Some(pseduoalign(&itr.next().unwrap(), &reference_index, MATCH_THRESHOLD));
      }

      // Get the score and the associated equivalence class
      let mut max_score = 0;
      let mut max_eqv_class = Vec::new();
      if let Some((eqv_class, score)) = seq_score {
        max_score = score;
        max_eqv_class = eqv_class;
      } 

      // If there's a reverse sequence and it matches, compare to the previous score and replace it if this score is higher
      if let Some(Some((eqv_class, score))) = rev_seq_score {
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

  // Create reference library iterator and filter results by lineage.
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).unwrap());
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