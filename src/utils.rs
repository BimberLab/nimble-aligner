use crate::reference_library; 

use reference_library::ReferenceMetadata;
use std::path;
use std::io::{Write, Read, Error, ErrorKind};
use std::fs::File;
use csv::Reader;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

// Takes a reader and returns a csv reader that wraps it with tab delimiters
pub fn get_tsv_reader<R: Read>(reader: R) -> Reader<R>{
  csv::ReaderBuilder::new()
    .delimiter(b'\t')
    .from_reader(reader)
}


// Takes the path to a fastq file and returns an error-checked iterator of the DnaStrings of the file
pub fn get_error_checked_fastq_reader(file_path: &str) -> impl Iterator<Item = Result<DnaString, Error>> {
  fastq::Reader::from_file(path::Path::new(file_path))
    .expect("Error -- cannot read sequence file: ")
    .records()
    .map(|record| match record { 
      Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
      _ => Err(Error::new(ErrorKind::InvalidData, "Unable to read sequence"))
    })
}


/* Takes a vector of potential reference genome results and an iterator to the reference library TSV.
 * Produces 2 vectors of sequence-name pairs. Panics if there is a reference sequence that cannot be read.
 * If they can be read, converts the given sequence to a DnaString and get the associated name. */
pub fn validate_reference_pairs<'a>(reference_genome: bio::io::fasta::Records<File>, 
  mut reference_library: impl Iterator<Item = &'a String>) -> (Vec<DnaString>, Vec<String>) {

  let mut reference_seqs: Vec<DnaString> = Vec::new();
  let mut reference_names: Vec<String> = Vec::new();

  for (i, reference) in reference_genome.enumerate() {
    reference_seqs.push(DnaString::from_acgt_bytes(reference.expect(&format!("Error -- could not read reference sequence #{}", i)).seq()));
    reference_names.push(reference_library.next().expect(&format!("Error -- could not read library name #{} after JSON parse, corrupted internal state.", i)).clone());
  }

  (reference_seqs, reference_names)
}


// Takes a result from the filtration pipeline and converts all the scores to percent matches
pub fn convert_scores_to_percentage(scores: Vec<(String, i32)>, total_reads: usize) -> Vec<(String, f32)> {
  scores.iter().map(|(name, score)| (name.to_string(), (*score as f32 / total_reads as f32) * 100.0)).collect()
}


// Write the given vector of tuples to a TSV file
pub fn write_to_tsv(results: Vec<(String, f32)>, mut reference_metadata: ReferenceMetadata) {
  let mut str_rep = String::new();

  // If we have a group_on index that isn't nt_sequence, remove nt_sequence from the header list and drop the relevant column
  if reference_metadata.group_on != reference_metadata.nt_sequence_idx {
    reference_metadata.headers.retain(|header| header != "nt_sequence");
    reference_metadata.columns.remove(reference_metadata.nt_sequence_idx);
  }

  // Remove the unnecessary nt_length metadata from header and columns
  let nt_len_idx = reference_metadata.headers.iter().position(|header| header == "nt_length").expect("Error -- no header nt_length found when writing results to disk");
  reference_metadata.headers.retain(|header| header != "nt_length");
  reference_metadata.columns.remove(nt_len_idx);

  // Add the headers to the top of the string representation of the tsv file
  str_rep += &(reference_metadata.headers.join("\t") + "\tmatch percentage\n");

  for (group, score) in results {
    // get relevant row idx by matching on the group name
    let row_idx = reference_metadata.columns[reference_metadata.group_on].iter().position(|name| &group == name).expect(&format!("Error -- group {} not found when writing results to disk", group));

    // iterate the columns, adding the element at row_idx to the string
    for column in &reference_metadata.columns {
      str_rep += &column[row_idx];
      str_rep += "\t";
    }

    // append the score, since it isn't in the reference metadata
    str_rep += &score.to_string();
    str_rep += "\n";
  }

  let mut file = File::create("results.tsv").expect("Error -- could not create results file");
  file.write_all(str_rep.as_bytes()).expect("Error -- could not write results to file");
}


pub fn sort_score_vector<T>(mut scores: Vec<(String, T)>) -> Vec<(String, T)> {
  scores.sort_by(|a, b| a.0.cmp(&b.0));
  scores
}