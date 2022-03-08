use crate::reference_library::ReferenceMetadata;
use crate::align::AlignDebugInfo;
use bio::alphabets::{dna, rna};
use csv::Reader;
use debruijn::dna_string::DnaString;
use std::fs::{OpenOptions, File};
use std::io::{Read, Write};
use unwrap::unwrap;
use flate2::{Compression, GzBuilder};

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
pub fn write_to_tsv(results: Vec<(Vec<String>, i32)>, group_row: Option<Vec<String>>, write_header: bool, output_path: &str) {
    let mut str_rep = String::new();

    // Add the headers to the top of the string representation of the tsv file
    if write_header {
        str_rep += "ambiguity class\tscore";

        match group_row {
            Some(ref _s) => {
                str_rep += "\t";
                str_rep += "cell barcode";
            },
            None => ()
        }

        str_rep += "\n";
    }

    let group_row_iter = match group_row {
        Some(ref vec) => vec.clone(),
        None => Vec::new()
    };
    let mut group_row_iter = group_row_iter.iter();

    // Append the results to the tsv string
    for (group, score) in results {
        str_rep += &group.join(",");
        str_rep += "\t";
        str_rep += &score.to_string();

        match group_row {
            Some(ref _vec) => {
                str_rep += "\t";
                str_rep += group_row_iter.next().unwrap();
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

pub fn write_debug_info(info: AlignDebugInfo) {
    println!("Writing debug info");

    let mut str_rep = String::new();

    str_rep += "Read units aligned: "; str_rep += &info.read_units_aligned.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to being below score threshold: "; str_rep += &info.score_below_threshold.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to multiple match: "; str_rep += &info.discarded_multiple_match.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to non-zero mismatches: "; str_rep += &info.discarded_nonzero_mismatch.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to not matching the reference library: "; str_rep += &info.no_match.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to not matching the reference library and having a low score: "; str_rep += &info.no_match_and_score_below_threshold.to_string(); str_rep += "\n";
    str_rep += "Units filtered for different reasons between the forward and reverse read: "; str_rep += &info.different_filter_reasons.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to non-matching pair alignments: "; str_rep += &info.not_matching_pair.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to a failed force intersect: "; str_rep += &info.force_intersect_failure.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to completely disjoint alignments: "; str_rep += &info.disjoint_pair_intersection.to_string(); str_rep += "\n";
    str_rep += "Units filtered due to mangled empty scores: "; str_rep += &info.best_class_empty.to_string(); str_rep += "\n";
    str_rep += "Forward-running alignments discarded: "; str_rep += &info.forward_runs_discarded.to_string(); str_rep += "\n";
    str_rep += "Reverse-running alignments discarded: "; str_rep += &info.backward_runs_discarded.to_string(); str_rep += "\n";
    str_rep += "Reverse-read sets discarded due to mangled paired-end data: "; str_rep += &info.reverse_read_sets_discarded_noneven.to_string(); str_rep += "\n";
    str_rep += "Reads discarded due to being too short after processing: "; str_rep += &info.short_read.to_string(); str_rep += "\n";

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(info.debug_file)
        .expect("Error -- could not create debug file");

    file.write_all(str_rep.as_bytes())
        .expect("Error -- could not write debug info to file");
}

pub fn filter_scores(reference_scores: Vec<(Vec<String>, i32)>, score_filter: &i32) -> Vec<(Vec<String>, i32)> {
    // Remove scores below the score threshold
    let reference_scores: Vec<(Vec<String>, i32)> = reference_scores
        .into_iter()
        .filter(|(_, val)| val > score_filter)
        .collect();

    reference_scores
}


pub fn write_read_list(results: Vec<(Vec<String>, String)>, output_path: &str) {
    let mut str_rep = String::new();

    // Append the results to the tsv string
    for (group, score) in results {
        str_rep += &group.join(",");
        str_rep += "\t";
        str_rep += &score;
        str_rep += "\n";
    }

    let f = File::create(output_path).expect("Could not create output file path for alignment metadata");
    let mut gz = GzBuilder::new()
                            .filename(output_path)
                            .write(f, Compression::default());
    gz.write(str_rep.as_bytes()).expect("Could not write to alignment metadata file.");
    gz.finish().expect("Could not flush to alignment metadata file buffer.");
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
