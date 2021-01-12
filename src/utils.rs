use std::path;
use std::io::{Write, Read, Error, ErrorKind};
use std::fs::File;
use csv::Reader;
use unwrap::unwrap;
use bio::io::fastq;
use debruijn::dna_string::DnaString;

// Takes a reader and returns a csv reader that wraps it, configures to use tab delimiters
pub fn get_tsv_reader<R: Read>(reader: R) -> Reader<R>{
  csv::ReaderBuilder::new()
    .delimiter(b'\t')
    .from_reader(reader)
}


// Takes the path to a fastq.gz file and returns an error-checked iterator of the DnaStrings of the file
pub fn get_error_checked_fastq_reader(file_path: &str) -> impl Iterator<Item = Result<DnaString, Error>> {
  let (reader, _) = unwrap!(niffler::from_path(path::Path::new(file_path)),
    "Error -- could not determine compression format for {}", file_path);

  //let mut contents = vec![];
  //unwrap!(reader.read_to_end(&mut contents), "Error -- could not read full contents of {}", file_path);

  fastq::Reader::new(reader)
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
    reference_seqs.push(DnaString::from_acgt_bytes(unwrap!(reference, "Error -- could not read reference sequence #{}", i).seq()));
    reference_names.push(unwrap!(reference_library.next(), "Error -- could not read library name #{} after JSON parse, corrupted internal state.", i).clone());
  }

  (reference_seqs, reference_names)
}


// Takes a result from the filtration pipeline and appends match percentages to the score tuples
pub fn append_match_percent(scores: Vec<(Vec<String>, i32)>, total_hits: usize) -> Vec<(Vec<String>, i32, f32)> {
  scores.iter().map(|(names, score)| (names.clone(), *score, (*score as f32 / total_hits as f32) * 100.0)).collect()
}


// Write the given vector of scores to a TSV file
pub fn write_to_tsv(results: Vec<(Vec<String>, i32)>) {
  let mut str_rep = String::new();

  // Add the headers to the top of the string representation of the tsv file
  str_rep += "ambiguity class\tscore\n";

  // Append the results to the tsv string
  for (group, score) in results {
    str_rep += &group.join(",");
    str_rep += "\t";
    str_rep += &score.to_string();
    str_rep += "\n";
  }

  let mut file = File::create("results.tsv").expect("Error -- could not create results file");
  file.write_all(str_rep.as_bytes()).expect("Error -- could not write results to file");
}


// Take a score vector produced by utils::convert_scores_to_percentage() and sort them by name
pub fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
  scores.sort_by(|a, b| a.0.cmp(&b.0));
  scores
}