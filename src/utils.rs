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
 * Produces 2 vectors of sequence-name pairs. Skip and log references that cannot be read.
 * If they can be read, converts the given sequence to a DnaString and get the associated name. */
pub fn get_valid_reference_pairs(reference_genome: bio::io::fasta::Records<File>, 
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


// Takes a result from the filtration pipeline and converts all the scores to percent matches
pub fn convert_scores_to_percentage(scores: Vec<(String, i32)>, total_reads: usize) -> Vec<(String, f32)> {
  scores.iter().map(|(name, score)| (name.to_string(), (*score as f32 / total_reads as f32) * 100.0)).collect()
}


// Write the given vector of tuples to a TSV file
pub fn write_to_tsv(results: Vec<(String, f32)>) {
  let mut str_rep = String::new();

  str_rep += "lineage\tmatch percentage\n";

  for (group, score) in results {
    str_rep += &group.to_string();
    str_rep += "\t";
    str_rep += &score.to_string();
    str_rep += "\n";
  }

  let mut file = File::create("results.tsv").expect("Error -- could not create results file");
  file.write_all(str_rep.as_bytes()).expect("Error -- could not write results to file");
}