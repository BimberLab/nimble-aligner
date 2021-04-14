use crate::reference_library::ReferenceMetadata;
use bio::alphabets::rna::revcomp;
use bio::io::fastq;
use csv::Reader;
use debruijn::dna_string::DnaString;
use std::fs::File;
use std::io::{Error, ErrorKind, Read, Write};
use std::path;
use unwrap::unwrap;

pub type DNAStringIter = impl Iterator<Item = Result<DnaString, Error>>;

// Takes a reader and returns a csv reader that wraps it, configures to use tab delimiters
pub fn get_tsv_reader<R: Read>(reader: R) -> Reader<R> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader)
}

// Takes the path to a fastq.gz file and returns an error-checked iterator of the DnaStrings of the file
pub fn get_error_checked_fastq_readers(
    file_path: &str,
) -> (DNAStringIter, DNAStringIter) {
    (get_error_checked_fastq_reader(file_path), get_error_checked_fastq_reader(file_path))
}

fn get_error_checked_fastq_reader(file_path: &str) -> DNAStringIter {
    let (reader, _) = unwrap!(
        niffler::from_path(path::Path::new(file_path)),
        "Error -- could not determine compression format for {}",
        file_path
    );

    fastq::Reader::new(reader)
        .records()
        .map(|record| match record {
            Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
            _ => Err(Error::new(
                ErrorKind::InvalidData,
                "Unable to read sequence",
            )),
        })
}

/* Takes a reference to the ReferenceMetadata structure.
 * Produces 3 vectors of sequence-name pairs. Panics if there is a reference sequence that cannot be read.
 * If they can be read, converts the given sequence to a DnaString and get the associated name. */
pub fn validate_reference_pairs(
    reference: &ReferenceMetadata,
) -> (Vec<DnaString>, Vec<DnaString>, Vec<String>) {
    let reference_genome = reference.columns[reference.sequence_idx].iter();
    let mut reference_library = reference.columns[reference.sequence_name_idx].iter();

    let mut reference_seqs: Vec<DnaString> = Vec::new();
    let mut reference_seqs_rev: Vec<DnaString> = Vec::new();
    let mut reference_names: Vec<String> = Vec::new();

    for (i, reference) in reference_genome.enumerate() {
        reference_seqs.push(DnaString::from_acgt_bytes(reference.as_bytes()));
        reference_seqs_rev.push(DnaString::from_acgt_bytes(&revcomp(
            reference.clone().into_bytes(),
        )));
        reference_names.push(unwrap!(reference_library.next(), "Error -- could not read library name #{} after JSON parse, corrupted internal state.", i).clone());
    }

    (reference_seqs, reference_seqs_rev, reference_names)
}

// Takes a result from the filtration pipeline and appends match percentages to the score tuples
pub fn append_match_percent(
    scores: Vec<(Vec<String>, i32)>,
    total_hits: usize,
) -> Vec<(Vec<String>, i32, f32)> {
    scores
        .iter()
        .map(|(names, score)| {
            (
                names.clone(),
                *score,
                (*score as f32 / total_hits as f32) * 100.0,
            )
        })
        .collect()
}

// Write the given vector of scores to a TSV file
pub fn write_to_tsv(results: Vec<(Vec<String>, i32)>, output_path: &str) {
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

    let mut file = File::create(output_path).expect("Error -- could not create results file");
    file.write_all(str_rep.as_bytes())
        .expect("Error -- could not write results to file");
}

// Take a score vector produced by utils::convert_scores_to_percentage() and sort them by name
pub fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}
