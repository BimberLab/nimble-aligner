use crate::reference_library::Reference;
use csv::Reader;
use debruijn::dna_string::DnaString;
use flate2::{Compression, GzBuilder};
use std::fs::OpenOptions;
use std::io::{Read, Write};
use unwrap::unwrap;

// Takes a reader and returns a csv reader that wraps it, configures to use tab delimiters
pub fn get_tsv_reader<R: Read>(reader: R) -> Reader<R> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader)
}

/* Using a Reference structure, produce 2 vectors of sequence data (normal and reverse comp), as well as a list of sequence names. */
pub fn get_reference_sequence_data(
    reference: &Reference,
) -> (Vec<DnaString>, Vec<DnaString>, Vec<String>) {
    let sequences = reference.columns[reference.sequence_idx].iter();
    let mut sequence_names = reference.columns[reference.sequence_name_idx].iter();

    let mut sequence_dnastrings: Vec<DnaString> = Vec::new();
    let mut sequence_dnastrings_revcomp: Vec<DnaString> = Vec::new();
    let mut reference_names: Vec<String> = Vec::new();

    for (i, feature_sequence) in sequences.enumerate() {
        sequence_dnastrings.push(DnaString::from_acgt_bytes(feature_sequence.as_bytes()));
        sequence_dnastrings_revcomp.push(DnaString::from_dna_string(&revcomp(feature_sequence)));
        reference_names.push(unwrap!(sequence_names.next(), "Error -- could not read library name #{} after JSON parse, corrupted internal state.", i).clone());
    }

    (sequence_dnastrings, sequence_dnastrings_revcomp, reference_names)
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
    mut scores: Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
) -> Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))> {
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

pub fn revcomp(sequence: &str) -> String {
    let mut revcomp_sequence_buffer: String = String::with_capacity(sequence.len());

    for bp in sequence.chars().rev() {
        match is_valid_base_pair(bp) {
            false => panic!("Input sequence base is not DNA: {}", bp),
            true => revcomp_sequence_buffer.push(revcomp_base_pair(bp)),
        }
    }
    revcomp_sequence_buffer
}

fn revcomp_base_pair(bp: char) -> char {
    match bp {
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

fn is_valid_base_pair(bp: char) -> bool {
    match bp {
        'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'U' | 'u' | 'N' | 'n' => true,
        _ => false,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_reference_sequence_data() {
        let reference = Reference {
            group_on: 0,
            headers: vec!["id".to_string(), "sequence_name".to_string(), "sequence".to_string()],
            columns: vec![
                vec!["1".to_string(), "2".to_string()],
                vec!["gene1".to_string(), "gene2".to_string()],
                vec!["ATGC".to_string(), "CGTA".to_string()]
            ],
            sequence_name_idx: 1,
            sequence_idx: 2,
        };

        let (dnastrings, dnastrings_revcomp, names) = get_reference_sequence_data(&reference);

        assert_eq!(dnastrings.len(), 2);
        assert_eq!(dnastrings_revcomp.len(), 2);
        assert_eq!(names, vec!["gene1", "gene2"]);
        assert_eq!(dnastrings[0].to_string(), "ATGC");
        assert_eq!(dnastrings[1].to_string(), "CGTA");
        assert_eq!(dnastrings_revcomp[0].to_string(), "GCAT");
        assert_eq!(dnastrings_revcomp[1].to_string(), "TACG");
    }

    #[test]
    #[should_panic(expected = "Error -- could not read library name")]
    fn test_get_reference_sequence_data_panic_on_missing_name() {
        let reference = Reference {
            group_on: 0,
            headers: vec!["id".to_string(), "sequence_name".to_string(), "sequence".to_string()],
            columns: vec![
                vec!["1".to_string()],
                vec!["gene1".to_string()],  // Missing second name
                vec!["ATGC".to_string(), "CGTA".to_string()]
            ],
            sequence_name_idx: 1,
            sequence_idx: 2,
        };

        get_reference_sequence_data(&reference);
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp("ATGC"), "GCAT");
        assert_eq!(revcomp("CCGGTTAA"), "TTAACCGG");
    }

    #[test]
    #[should_panic(expected = "Input sequence base is not DNA")]
    fn test_revcomp_invalid_input() {
        revcomp("ATGX");
    }

    #[test]
    fn test_revcomp_base_pair() {
        assert_eq!(revcomp_base_pair('A'), 'T');
        assert_eq!(revcomp_base_pair('T'), 'A');
        assert_eq!(revcomp_base_pair('C'), 'G');
        assert_eq!(revcomp_base_pair('G'), 'C');
        assert_eq!(revcomp_base_pair('U'), 'A');

        assert_eq!(revcomp_base_pair('a'), 't');
        assert_eq!(revcomp_base_pair('t'), 'a');
        assert_eq!(revcomp_base_pair('c'), 'g');
        assert_eq!(revcomp_base_pair('g'), 'c');
        assert_eq!(revcomp_base_pair('u'), 'a');

        assert_eq!(revcomp_base_pair('N'), 'N');
        assert_eq!(revcomp_base_pair('n'), 'N');
    }

    #[test]
    fn test_revcomp_base_pair_invalid() {
        assert_eq!(revcomp_base_pair('X'), 'N');
        assert_eq!(revcomp_base_pair('Y'), 'N');
        assert_eq!(revcomp_base_pair('Z'), 'N');
        assert_eq!(revcomp_base_pair('1'), 'N');
        assert_eq!(revcomp_base_pair('#'), 'N');
    }

    #[test]
    fn test_is_valid_base_pair() {
        // Valid DNA base pairs
        assert!(is_valid_base_pair('A'));
        assert!(is_valid_base_pair('T'));
        assert!(is_valid_base_pair('C'));
        assert!(is_valid_base_pair('G'));

        // Valid DNA base pairs, lower case
        assert!(is_valid_base_pair('a'));
        assert!(is_valid_base_pair('t'));
        assert!(is_valid_base_pair('c'));
        assert!(is_valid_base_pair('g'));

        // Valid RNA base pair
        assert!(is_valid_base_pair('U'));
        assert!(is_valid_base_pair('u'));

        // Nucleotide ambiguity code
        assert!(is_valid_base_pair('N'));
        assert!(is_valid_base_pair('n'));

        // Invalid characters
        assert!(!is_valid_base_pair('X'));
        assert!(!is_valid_base_pair('Y'));
        assert!(!is_valid_base_pair('Z'));
        assert!(!is_valid_base_pair('1'));
        assert!(!is_valid_base_pair('#'));
    }

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
