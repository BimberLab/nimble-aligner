use crate::align::AlignDebugInfo;
use crate::align::BamData;
use crate::reference_library::ReferenceMetadata;
use csv::Reader;
use debruijn::dna_string::DnaString;
use flate2::{Compression, GzBuilder};
use std::fs::OpenOptions;
use std::io::{Read, Write};
use unwrap::unwrap;

fn truncate(s: &str, max_chars: usize) -> &str {
    match s.char_indices().nth(max_chars) {
        None => s,
        Some((idx, _)) => &s[..idx],
    }
}

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
        "DNA" => revcomp,
        "RNA" => revcomp,
        _ => panic!(
            "Error -- cannot determine revcomp method to use -- ensure data_type is a valid type"
        ),
    };

    for (i, reference) in reference_genome.enumerate() {
        reference_seqs.push(DnaString::from_acgt_bytes(reference.as_bytes()));
        reference_seqs_rev.push(DnaString::from_dna_string(&revcomp(reference)));
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

pub fn write_to_tsv(
    results: &Vec<(Vec<String>, i32)>,
    output_path: String
) {
    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(&output_path)
        .expect("Unable to open file");

    // Check if the file is empty to decide whether to write the header
    if file.metadata().expect("Unable to read file metadata").len() == 0 {
        // Write the header if the file is empty
        writeln!(file, "feature\tscore").expect("Unable to write header");
    }

    // Iterate over the results and write each as a TSV row
    for (features, score) in results {
        // Join all the features with a tab character
        let feature_str = features.join("\t");
        // Write the row to the file
        writeln!(file, "{}\t{}", feature_str, score).expect("Unable to write row");
    }
}


// Take a score vector produced by utils::convert_scores_to_percentage() and sort them by name
pub fn sort_score_vector(
    mut scores: Vec<(Vec<String>, (i32, BamData, BamData))>,
) -> Vec<(Vec<String>, (i32, BamData, BamData))> {
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

    str_rep += "Read units aligned: ";
    str_rep += &info.read_units_aligned.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.read_units_aligned as f64 / info.get_total_attempted_reads() as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to being below score threshold: ";
    str_rep += &info.score_below_threshold.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.score_below_threshold as f64 / info.get_total_attempted_reads() as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to multiple match: ";
    str_rep += &info.discarded_multiple_match.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.discarded_multiple_match as f64 / info.get_total_attempted_reads() as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to non-zero mismatches: ";
    str_rep += &info.discarded_nonzero_mismatch.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.discarded_nonzero_mismatch as f64 / info.get_total_attempted_reads() as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to not matching the reference library: ";
    str_rep += &info.no_match.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.no_match as f64 / info.get_total_attempted_reads() as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to non-matching pair alignments: ";
    str_rep += &info.not_matching_pair.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.not_matching_pair as f64 / info.get_total_attempted_reads() as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Reads discarded due to being too short after processing: ";
    str_rep += &info.short_read.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.short_read as f64 / info.get_total_attempted_reads() as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads attempted)";
    str_rep += "\n";

    str_rep += "Units filtered due to not matching the reference library and having a low score: ";
    str_rep += &info.no_match_and_score_below_threshold.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.no_match_and_score_below_threshold as f64 / info.get_total_attempted_reads() as f64
            * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Units filtered for different reasons between the forward and reverse read: ";
    str_rep += &info.different_filter_reasons.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.different_filter_reasons as f64 / info.read_units_aligned as f64 * 100.0)
            .to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Units filtered due to a failed force intersect: ";
    str_rep += &info.force_intersect_failure.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.force_intersect_failure as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Reads discarded due to equivalence class exceeding the max hits to report: ";
    str_rep += &info.max_hits_exceeded.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.max_hits_exceeded as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction FF reported: ";
    str_rep += &info.ff_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.ff_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction RR reported: ";
    str_rep += &info.rr_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.rr_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction UU reported: ";
    str_rep += &info.uu_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.uu_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction FR reported: ";
    str_rep += &info.fr_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.fr_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction FU reported: ";
    str_rep += &info.fu_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.fu_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction RF reported: ";
    str_rep += &info.rf_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.rf_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction RU reported: ";
    str_rep += &info.ru_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.ru_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction UF reported: ";
    str_rep += &info.uf_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.uf_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Read direction UR reported: ";
    str_rep += &info.ur_reported.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.ur_reported as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";
    str_rep += "\n";

    str_rep += "Reads dropped due to lacking corrected cell barcode: ";
    str_rep += &info.number_cr_skipped.to_string();
    str_rep += " (";
    str_rep += truncate(
        &(info.number_cr_skipped as f64 / info.read_units_aligned as f64 * 100.0).to_string(),
        6,
    );
    str_rep += "% of reads aligned)";

    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(false)
        .open(info.debug_file)
        .expect("Error -- could not create debug file");

    file.write_all(str_rep.as_bytes())
        .expect("Error -- could not write debug info to file");
}

pub fn filter_scores(
    reference_scores: Vec<(Vec<String>, i32)>,
    score_filter: &i32,
) -> Vec<(Vec<String>, i32)> {
    // Remove scores below the score threshold
    let reference_scores: Vec<(Vec<String>, i32)> = reference_scores
        .into_iter()
        .filter(|(_, val)| val > score_filter)
        .collect();

    reference_scores
}

pub struct PseudoalignerData {
    pub reference_names: Vec<Vec<String>>,
    pub read_umi_name: Vec<String>,
    pub barcode_sample_name: Vec<String>,
    pub score: Vec<f64>,
    pub raw_score: Vec<usize>,
    pub pair: Vec<String>,
    pub sequence: Vec<String>,
    pub strand_filter_reason: Vec<String>,
}

pub struct BamSpecificAlignMetadata {
    pub mapq: Vec<u8>,
    pub orientation: Vec<String>,
    pub hits: Vec<String>,
    pub qnames: Vec<Vec<u8>>,
    pub quals: Vec<Vec<u8>>,
    pub txs: Vec<String>,
}

pub fn write_read_list(
    pseudoaligner_data: &PseudoalignerData,
    bam_data: Option<&BamSpecificAlignMetadata>,
    output_path: &str,
) {
    let mut str_rep = String::new();

    // Append the results to the tsv string
    for (i, _) in pseudoaligner_data.reference_names.iter().enumerate() {
        if !(i < pseudoaligner_data.reference_names.len()
            && i < pseudoaligner_data.read_umi_name.len()
            && i < pseudoaligner_data.barcode_sample_name.len()
            && i < pseudoaligner_data.score.len()
            && i < pseudoaligner_data.pair.len()
            && i < pseudoaligner_data.sequence.len()
            && i < pseudoaligner_data.strand_filter_reason.len())
        {
            println!("Debug data truncated due to indexing error");
            return;
        }

        str_rep += &pseudoaligner_data.reference_names[i].join(",");
        str_rep += "\t";
        str_rep += &pseudoaligner_data.read_umi_name[i];
        str_rep += "\t";
        str_rep += &pseudoaligner_data.barcode_sample_name[i];
        str_rep += "\t";
        str_rep += &pseudoaligner_data.score[i].to_string();
        str_rep += "\t";
        str_rep += &pseudoaligner_data.raw_score[i].to_string();
        str_rep += "\t";
        str_rep += &pseudoaligner_data.pair[i];
        str_rep += "\t";
        str_rep += &pseudoaligner_data.sequence[i];
        str_rep += "\t";
        str_rep += &pseudoaligner_data.strand_filter_reason[i];

        if let Some(ref metadata) = bam_data {
            if metadata.mapq.len() > 0 && i < metadata.mapq.len() {
                str_rep += "\t";
                str_rep += &metadata.mapq[i].to_string();
            }

            if metadata.orientation.len() > 0 && i < metadata.orientation.len() {
                str_rep += "\t";
                str_rep += &metadata.orientation[i].to_string();
            }

            if metadata.hits.len() > 0 && i < metadata.hits.len() {
                str_rep += "\t";
                str_rep += &metadata.hits[i].to_string();
            }
        }

        str_rep += "\n";
    }

    let f = OpenOptions::new()
        .append(true)
        .create(true)
        .open(output_path)
        .expect("Could not create output file path for alignment metadata");

    let mut gz = GzBuilder::new()
        .filename(output_path)
        .write(f, Compression::default());
    gz.write(str_rep.as_bytes())
        .expect("Could not write to alignment metadata file.");
    gz.finish()
        .expect("Could not flush to alignment metadata file buffer.");
}

pub fn revcomp(dna: &str) -> String {
    // result vector
    let mut rdna: String = String::with_capacity(dna.len());

    // iterate through the input &str
    for c in dna.chars().rev() {
        // test the input
        match is_dna(c) {
            false => panic!("Input sequence base is not DNA: {}", dna),
            true => rdna.push(switch_base(c)),
        }
    }
    rdna
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N',
    }
}

pub fn shannon_entropy(dna: &str) -> f64 {
    let total_length = dna.len() as f64;

    let mut frequencies = [0.0; 4]; // A, T, C, G

    for char in dna.chars() {
        match char {
            'A' => frequencies[0] += 1.0,
            'T' => frequencies[1] += 1.0,
            'C' => frequencies[2] += 1.0,
            'G' => frequencies[3] += 1.0,
            _ => (),
        }
    }

    frequencies.iter_mut().for_each(|f| *f /= total_length);

    let entropy: f64 = frequencies.iter()
                                  .filter(|&&f| f > 0.0) // Exclude zero frequencies to avoid NaN
                                  .map(|&f| f * f.log2())
                                  .sum();

    -entropy
}

fn is_dna(dna: char) -> bool {
    match dna {
        'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'U' | 'u' | 'N' | 'n' => true,
        _ => false,
    }
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
