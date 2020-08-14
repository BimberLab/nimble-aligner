mod filter;
mod utils;

use std::path;
use std::fs::File;
use std::io::{Write, Error};
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

fn main() {
  print!("Loading and preprocessing data\n");
  // Parse arguments based on the yaml schema
  let yaml = load_yaml!("cli.yml");
  let matches = App::from_yaml(yaml).get_matches();

  let library = matches.value_of("library").unwrap();
  let library_fasta = matches.value_of("library_fasta").unwrap();
  let input_files: Vec<&str> = matches.values_of("input").unwrap().collect();
  let num_cores = matches.value_of("num_cores").unwrap_or("1").parse::<usize>().expect("Error -- please provide an integer value for the number of cores");

  // Parameters for alignment
  const MATCH_THRESHOLD: usize = 60;
  const REFERENCE_GENOME_SIZE: usize = 1209;
  const READS_SIZE: usize = 2786342;

  // Get iterators to records in the reference genome file and library
  let reference_genome: Vec<Result<fasta::Record, Error>> = fasta::Reader::from_file(library_fasta).expect(
    "Error -- could not read reference genome"
  ).records().collect();

  let mut reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library"
  ));
  let mut reference_library = reference_library.records();

  // Load reference genome and library into memory. Skip and log references that cannot be read.
  let mut reference_seqs = Vec::new();
  let mut reference_names = Vec::new();
  for (i, reference) in reference_genome.iter().enumerate() {
    let reference_name = reference_library.next();

    if let Ok(reference) = reference {
      if let Some(Ok(record)) = reference_name {
        reference_seqs.push(DnaString::from_acgt_bytes(reference.seq())); // Convert raw data to DNAString
        reference_names.push(record[0].to_string());
      } else {
        print!("Warning: Could not read library name #{}\n", i);
        continue;
      }
    } else {
      print!("Warning: Could not read library reference #{}\n", i);
      continue;
    }
  }

  // Create debruijn-mapped index of the reference library
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    num_cores
  ).expect("Error -- could not create pseudoaligner index of the reference library");

  // Get iterator to the sequences that will be aligned to the reference from the sequence genome file(s)
  // If a given sequence can be read, convert it to a DnaString. Otherwise, report an error.
  let sequences = fastq::Reader::from_file(path::Path::new(input_files[0]))
    .expect("Error -- cannot read sequence file: ")
    .records()
    .map(|record| match record { 
      Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
      _ => Err("Unable to read sequence")
    });

  let mut reverse_sequences = if input_files.len() > 1 {
    Some(fastq::Reader::from_file(path::Path::new(input_files[1]))
          .expect("Error -- cannot read reverse sequence file: ")
          .records()
          .map(|record| match record { 
            Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
            _ => Err("Unable to read sequence") 
          }))
  } else {
    None
  };


  print!("Pseudo-aligning reads to reference index\n");
  // Map over the sequence iterator, producing an iterator of score vectors for each sequence
  let reference_scores = sequences.fold(vec![0.0; REFERENCE_GENOME_SIZE], 
    |mut acc, read| {

      /* Generate score and equivalence class for this read by aligning the sequence against
       * the current reference. This alignment returns any scores that are greater than the match threshold. */
      let seq_score = pseduoalign(read, &reference_index, MATCH_THRESHOLD);
      let mut rev_seq_score = None;

      // If there's a reversed sequence, do the paired-end alignment
      if let Some(itr) = &mut reverse_sequences {
        let reverse_read = itr.next().expect("Error -- read and reverse read files do not have matching lengths: ");
        rev_seq_score = Some(pseduoalign(reverse_read, &reference_index, MATCH_THRESHOLD));
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
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library"
  ));


  print!("Filtering results by lineage\n");
  let results = filter::collapse_results_by_lineage(reference_library.into_records(), reference_scores.iter());


  print!("Writing results to file\n");
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

  let mut file = File::create("results.tsv").expect("Error -- could not create results file");
  file.write_all(str_rep.as_bytes()).expect("Error -- could not write results to file");

  print!("Output results written to ./results.tsv")
}

// Align the given sequence against the given reference with a score threshold
fn pseduoalign(sequence: Result<DnaString, &str>, reference_index: &PseudoAligner, match_threshold: usize) -> Option<(Vec<u32>, usize)> {
  if sequence.is_err() {
    return None;
  }

  match reference_index.map_read(&sequence.unwrap()) {
    Some((equiv_class, score)) => if score >= match_threshold && !equiv_class.is_empty() { Some((equiv_class, score)) } else { None }
    None => None
  }
}