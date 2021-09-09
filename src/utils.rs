use crate::reference_library::ReferenceMetadata;
use bio::alphabets::{dna, rna};
use csv::Reader;
use debruijn::dna_string::DnaString;
use std::fs::OpenOptions;
use std::io::{Read, Write};
use unwrap::unwrap;

// Takes a reader and returns a csv reader that wraps it, configures to use tab delimiters
pub fn get_tsv_reader<R: Read>(reader: R) -> Reader<R> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader)
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

    let revcomp = match reference.data_type.as_str() {
        "DNA" => dna::revcomp,
        "RNA" => rna::revcomp,
        _ => panic!(
            "Error -- cannot determine revcomp method to use -- ensure data_type is a valid type"
        ),
    };

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
pub fn write_to_tsv(results: Vec<(Vec<String>, i32)>, group_column: Option<String>, write_header: bool, output_path: &str) {
    let mut str_rep = String::new();

    // Add the headers to the top of the string representation of the tsv file
    if write_header {
        str_rep += "ambiguity class\tscore";

        match group_column {
            Some(ref _s) => {
                str_rep += "\t";
                str_rep += "cell barcode";
            },
            None => ()
        }

        str_rep += "\n";
    }

    // Append the results to the tsv string
    for (group, score) in results {
        str_rep += &group.join(",");
        str_rep += "\t";
        str_rep += &score.to_string();

        match group_column {
            Some(ref s) => {
                str_rep += "\t";
                str_rep += s;
            },
            None => ()
        }

        str_rep += "\n";
    }

    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(output_path)
        .expect("Error -- could not create results file");

    file.write_all(str_rep.as_bytes())
        .expect("Error -- could not write results to file");
}

// Take a score vector produced by utils::convert_scores_to_percentage() and sort them by name
pub fn sort_score_vector(mut scores: Vec<(Vec<String>, i32)>) -> Vec<(Vec<String>, i32)> {
    scores.sort_by(|a, b| a.0.cmp(&b.0));
    scores
}

// Determine if a file is a .bam or a .fasta based on file extension
pub fn is_fastq(file: &str) -> bool {
    let mut is_fasta = true;
    let components = file.split(".").skip(1);

    for component in components {
        if component == "bam" {
            is_fasta = false;
        }
    }

    is_fasta
}

#[cfg(test)]
mod tests {
    #[test]
    fn is_fasta_short() {
        let expected_results = true;
        let results = super::is_fastq("reference.fastq");
        assert_eq!(results, expected_results);
    }

    #[test]
    fn is_fasta_long() {
        let expected_results = true;
        let results = super::is_fastq("reference.bin.fastq.gz");
        assert_eq!(results, expected_results);
    }

    #[test]
    fn is_bam_short() {
        let expected_results = false;
        let results = super::is_fastq("reference.bam.gz");
        assert_eq!(results, expected_results);
    }

    #[test]
    fn is_bam_long() {
        let expected_results = false;
        let results = super::is_fastq("reference.bin.bam.zip");
        assert_eq!(results, expected_results);
    }
}
