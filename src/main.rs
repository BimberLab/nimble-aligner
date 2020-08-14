mod filter;
mod utils;

use std::path;
use std::fs::File;
use std::io::{Write, Error, ErrorKind};
use std::collections::HashMap;
use clap::{App, load_yaml};
use bio::io::fasta;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

fn main() {
  println!("Loading and preprocessing data");

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
  let reference_genome  = fasta::Reader::from_file(library_fasta).expect(
    "Error -- could not read reference genome"
  ).records();

  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library"
  )).into_records();

  // Parse out reference/reference name pairs that cannot be read for whatever reason
  let (reference_seqs, reference_names) = get_valid_reference_pairs(reference_genome, reference_library);

  // Create debruijn-mapped index of the reference library
  let reference_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
    &reference_seqs,
    &reference_names,
    &HashMap::new(),
    num_cores
  ).expect("Error -- could not create pseudoaligner index of the reference library");

  println!("Reading sequences");

  /* Get error-checked iterators to the sequences that will be aligned to the reference from the
   * sequence genome file(s) */
  let sequences = get_error_checked_fastq_reader(input_files[0]);

  // Only get reverse sequences if a file is provided
  let reverse_sequences = if input_files.len() > 1 {
    println!("Reading reverse sequences");
    Some(get_error_checked_fastq_reader(input_files[1]))
  } else {
    None
  };

  println!("Pseudo-aligning reads to reference index");

  let reference_scores = score(sequences, reverse_sequences, reference_index, REFERENCE_GENOME_SIZE, MATCH_THRESHOLD);

  // Normalize the results vector down to a percentage
  let reference_scores: Vec<f32> = reference_scores.iter().map(|score| (score / READS_SIZE as f32) * 100.0).collect();

  // Create reference library iterator to filter the scores by lineage
  let reference_library = utils::get_tsv_reader(File::open(path::Path::new(library)).expect(
    "Error -- could not read reference library"
  )).into_records();

  println!("Filtering results by lineage");

  let results = filter::collapse_results_by_lineage(reference_library, reference_scores.iter());

  println!("Writing results to file");

  // Write filtered results to file
  let mut str_rep = String::new();

  str_rep += "lineage";
  str_rep += "\t";
  str_rep += "percent of all reads\n";

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


/* Takes a vector of potential reference genome results and an iterator to the reference library TSV.
 * Produces 2 vectors of sequence-name pairs. Skip and log references that cannot be read.
 * If they can be read, converts the given sequence to a DnaString and get the associated name. */
fn get_valid_reference_pairs(reference_genome: bio::io::fasta::Records<File>, 
  mut reference_library: csv::StringRecordsIntoIter<std::fs::File>) -> (Vec<DnaString>, Vec<String>) {

  let mut reference_seqs = Vec::new();
  let mut reference_names = Vec::new();

  for (i, reference) in reference_genome.enumerate() {
    let reference_name = reference_library.next();

    if let Ok(reference) = reference {
      if let Some(Ok(record)) = reference_name {
        reference_seqs.push(DnaString::from_acgt_bytes(reference.seq())); // Convert raw data to DNAString
        reference_names.push(record[0].to_string());
      } else {
        println!("Warning: Could not read library name #{}", i);
        continue;
      }
    } else {
      println!("Warning: Could not read library reference #{}", i);
      continue;
    }
  }

  (reference_seqs, reference_names)
}


// Takes the path to a fastq file and returns an error-checked iterator of the DnaStrings of the file
fn get_error_checked_fastq_reader(file_path: &str) -> impl Iterator<Item = Result<DnaString, Error>> {
  fastq::Reader::from_file(path::Path::new(file_path))
    .expect("Error -- cannot read sequence file: ")
    .records()
    .map(|record| match record { 
      Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
      _ => Err(Error::new(ErrorKind::InvalidData, "Unable to read sequence"))
    })
}


/* Takes a set of sequences and optionally, reverse sequences, a debrujin map index of the reference
 * genome, the size of the reference genome, and the threshold to match against, and performs a
 * debrujin-graph based pseduoalignment, returning a score for each readable reference in the reference
 * genome. */
fn score(sequences: impl Iterator<Item = Result<DnaString, Error>>, 
  mut reverse_sequences: Option<impl Iterator<Item = Result<DnaString, Error>>>,
  index: PseudoAligner, reference_genome_size: usize, match_threshold: usize) -> Vec<f32> {
  sequences.fold(vec![0.0; reference_genome_size], 
    |mut acc, read| {
      /* Generate score and equivalence class for this read by aligning the sequence against
       * the current reference. This alignment returns any scores that are greater than the match threshold. */
      let seq_score = pseduoalign(read, &index, match_threshold);
      let mut rev_seq_score = None;

      // If there's a reversed sequence, do the paired-end alignment
      if let Some(itr) = &mut reverse_sequences {
        let reverse_read = itr.next().expect("Error -- read and reverse read files do not have matching lengths: ");
        rev_seq_score = Some(pseduoalign(reverse_read, &index, match_threshold));
      }

      // Get the score and the associated equivalence class of the forward sequence
      let mut max_score = 0;
      let mut max_eqv_class = Vec::new();
      if let Some((eqv_class, score)) = seq_score {
        max_score = score;
        max_eqv_class = eqv_class;
      } 

      // If there's a reverse sequence and it matches, compare to the forward score and replace it if this score is higher
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
  )
}


// Align the given sequence against the given reference with a score threshold
fn pseduoalign(sequence: Result<DnaString, Error>, reference_index: &PseudoAligner, match_threshold: usize) -> Option<(Vec<u32>, usize)> {
  if sequence.is_err() {
    return None;
  }

  match reference_index.map_read(&sequence.unwrap()) {
    Some((equiv_class, score)) => if score >= match_threshold && !equiv_class.is_empty() { Some((equiv_class, score)) } else { None }
    None => None
  }
}