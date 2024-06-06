use crate::reference_library::Reference;
use debruijn::dna_string::DnaString;
use std::fs::OpenOptions;
use std::io::Write;

/* Using a Reference structure, produce a vector of sequence data, as well as a list of sequence names. Contains both normal and revcomp versions of the features*/
pub fn get_reference_sequence_data(
    reference: &Reference,
) -> (Vec<DnaString>, Vec<String>) {
    let sequences = reference.columns[reference.sequence_idx].iter();
    let mut sequence_names = reference.columns[reference.sequence_name_idx].iter();

    let mut sequence_dnastrings: Vec<DnaString> = Vec::new();
    let mut reference_names: Vec<String> = Vec::new();

    for feature_sequence in sequences {
        sequence_dnastrings.push(DnaString::from_acgt_bytes(feature_sequence.as_bytes()));

        let sequence_name = sequence_names.next().expect("Error -- could not read library name after JSON parse, corrupted internal state.");
        reference_names.push(sequence_name.clone());
    }

    (sequence_dnastrings, reference_names)
}

// Given a set of results from score() pipeline, write a TSV of the data to the given file
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
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use tempfile::NamedTempFile;

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

        let (dnastrings, names) = get_reference_sequence_data(&reference);

        assert_eq!(dnastrings.len(), 2);
        assert_eq!(names, vec!["gene1", "gene2"]);
        assert_eq!(dnastrings[0].to_string(), "ATGC");
        assert_eq!(dnastrings[1].to_string(), "CGTA");
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
    fn test_write_to_tsv_with_header() {
        let results = vec![
            (vec!["feature1".to_string(), "feature2".to_string()], 10),
            (vec!["feature3".to_string(), "feature4".to_string()], 20)
        ];

        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        write_to_tsv(&results, temp_path.clone());

        let file = File::open(&temp_path).unwrap();
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader.lines().collect::<Result<_, _>>().unwrap();

        assert_eq!(lines[0], "feature\tscore");
        assert_eq!(lines[1], "feature1\tfeature2\t10");
        assert_eq!(lines[2], "feature3\tfeature4\t20");
    }

    #[test]
    fn test_write_to_tsv_without_header() {
        let results = vec![
            (vec!["feature5".to_string(), "feature6".to_string()], 30)
        ];

        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        {
            let mut file = File::create(&temp_path).unwrap();
            writeln!(file, "feature\tscore").unwrap();
        }

        write_to_tsv(&results, temp_path.clone());

        let file = File::open(&temp_path).unwrap();
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader.lines().collect::<Result<_, _>>().unwrap();

        assert_eq!(lines[0], "feature\tscore");
        assert_eq!(lines[1], "feature5\tfeature6\t30");
    }

    #[test]
    fn test_sort_score_vector() {
        let scores = vec![
            (vec!["Charlie".to_string()], (90, vec!["A".to_string()], vec!["Fail".to_string()])),
            (vec!["Alice".to_string()], (95, vec!["A".to_string()], vec!["Pass".to_string()])),
            (vec!["Bob".to_string()], (85, vec!["B".to_string()], vec!["Pass".to_string()])),
        ];

        let sorted_scores = sort_score_vector(scores);

        assert_eq!(sorted_scores[0].0[0], "Alice");
        assert_eq!(sorted_scores[1].0[0], "Bob");
        assert_eq!(sorted_scores[2].0[0], "Charlie");

        assert_eq!(sorted_scores[0].1 .0, 95);
        assert_eq!(sorted_scores[1].1 .0, 85);
        assert_eq!(sorted_scores[2].1 .0, 90);
    }

    #[test]
    fn test_sort_score_vector_empty() {
        let scores: Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))> = vec![];
        let sorted_scores = sort_score_vector(scores.clone());
        assert_eq!(sorted_scores, scores);
    }

    #[test]
    fn test_sort_score_vector_single_element() {
        let scores = vec![(vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()]))];
        let sorted_scores = sort_score_vector(scores.clone());
        assert_eq!(sorted_scores, scores);
    }

    #[test]
    fn test_sort_score_vector_multiple_elements_sorted() {
        let scores = vec![
            (vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()])),
            (vec!["b".to_string()], (2, vec!["x2".to_string()], vec!["y2".to_string()])),
            (vec!["c".to_string()], (3, vec!["x3".to_string()], vec!["y3".to_string()])),
        ];
        let sorted_scores = sort_score_vector(scores.clone());
        assert_eq!(sorted_scores, scores);
    }

    #[test]
    fn test_sort_score_vector_multiple_elements_unsorted() {
        let scores = vec![
            (vec!["c".to_string()], (3, vec!["x3".to_string()], vec!["y3".to_string()])),
            (vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()])),
            (vec!["b".to_string()], (2, vec!["x2".to_string()], vec!["y2".to_string()])),
        ];
        let expected_sorted_scores = vec![
            (vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()])),
            (vec!["b".to_string()], (2, vec!["x2".to_string()], vec!["y2".to_string()])),
            (vec!["c".to_string()], (3, vec!["x3".to_string()], vec!["y3".to_string()])),
        ];
        let sorted_scores = sort_score_vector(scores);
        assert_eq!(sorted_scores, expected_sorted_scores);
    }

    #[test]
    fn test_sort_score_vector_multiple_elements_with_same_key() {
        let scores = vec![
            (vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()])),
            (vec!["a".to_string()], (2, vec!["x2".to_string()], vec!["y2".to_string()])),
            (vec!["b".to_string()], (3, vec!["x3".to_string()], vec!["y3".to_string()])),
        ];
        let expected_sorted_scores = vec![
            (vec!["a".to_string()], (1, vec!["x".to_string()], vec!["y".to_string()])),
            (vec!["a".to_string()], (2, vec!["x2".to_string()], vec!["y2".to_string()])),
            (vec!["b".to_string()], (3, vec!["x3".to_string()], vec!["y3".to_string()])),
        ];
        let sorted_scores = sort_score_vector(scores);
        assert_eq!(sorted_scores, expected_sorted_scores);
    }

    fn almost_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    #[test]
    fn test_shannon_entropy_empty_string() {
        let dna = "";
        let entropy = shannon_entropy(dna);
        assert!(almost_equal(entropy, 0.0, 1e-10));
    }

    #[test]
    fn test_shannon_entropy_single_nucleotide() {
        let dna = "A";
        let entropy = shannon_entropy(dna);
        assert!(almost_equal(entropy, 0.0, 1e-10));
    }

    #[test]
    fn test_shannon_entropy_two_nucleotides() {
        let dna = "AT";
        let entropy = shannon_entropy(dna);
        assert!(almost_equal(entropy, 1.0, 1e-10));
    }

    #[test]
    fn test_shannon_entropy_equal_frequencies() {
        let dna = "ATCG";
        let entropy = shannon_entropy(dna);
        assert!(almost_equal(entropy, 2.0, 1e-10));
    }

    #[test]
    fn test_shannon_entropy_unequal_frequencies() {
        let dna = "AAAT";
        let entropy = shannon_entropy(dna);
        let expected_entropy = -(0.75 * 0.75_f64.log2() + 0.25 * 0.25_f64.log2());
        assert!(almost_equal(entropy, expected_entropy, 1e-10));
    }

    #[test]
    fn test_shannon_entropy_realistic_sequence() {
        let dna = "ATCGATCGATCG";
        let entropy = shannon_entropy(dna);
        assert!(almost_equal(entropy, 2.0, 1e-10));
    }
}
