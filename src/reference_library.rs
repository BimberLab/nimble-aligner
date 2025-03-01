use crate::align::{self, LibraryChemistry};
use crate::utils::revcomp;
use serde_json::Value;
use std::fs::read_to_string;
use std::path::Path;
use unwrap::unwrap;

pub const SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR: &str = "§";

#[derive(Debug)]
pub struct Reference {
    pub group_on: usize,            // Index of a column defining related "feature families." All hits within a family count as one hit for the group.
    pub headers: Vec<String>,       // List of column names
    pub columns: Vec<Vec<String>>,  // Column data
    pub sequence_name_idx: usize,   // Index of the column containing feature names
    pub sequence_idx: usize,        // Index of the column containing sequence data
}

// Parses a .json that contains a reference library. Returns a tuple of the library's aligner config values and the library data
pub fn get_reference_library(path: &Path, strand_filter: LibraryChemistry) -> (align::AlignFilterConfig, Reference) {
    // Parse raw JSON to serde_json value
    let raw_json_string = read_to_string(path).expect("Error -- could not read reference library");

    let v: Value = serde_json::from_str(&raw_json_string)
        .expect("Error -- could not parse reference library JSON");

    // Get aligner configuration from the first JSON object in the file
    let aligner_config_obj = &v[0];
    let score_percent = aligner_config_obj["score_percent"]
        .as_f64()
        .expect("Error -- could not parse score_percent as f64")
        as f64;
    let score_filter = aligner_config_obj["score_filter"]
        .as_i64()
        .expect("Error -- could not parse score_filter as int64");
    let score_threshold = aligner_config_obj["score_threshold"]
        .as_i64()
        .expect("Error -- could not parse score_threshold as int64")
        as usize;
    let num_mismatches = aligner_config_obj["num_mismatches"]
        .as_i64()
        .expect("Error -- could not parse num_mismatches as int64")
        as usize;
    let discard_multiple_matches = aligner_config_obj["discard_multiple_matches"]
        .as_bool()
        .expect("Error -- could not parse discard_multiple_mismatches as boolean");
    let require_valid_pair = aligner_config_obj["require_valid_pair"]
        .as_bool()
        .expect("Error -- could not parse require_valid_pair as boolean");
    let discard_multi_hits = aligner_config_obj["discard_multi_hits"]
        .as_i64()
        .expect("Error -- could not parse discard_multi_hits as int64")
        as usize;
    let intersect_level = aligner_config_obj["intersect_level"]
        .as_i64()
        .expect("Error -- could not parse intersect_level as int64");
    let max_hits_to_report = aligner_config_obj["max_hits_to_report"]
        .as_i64()
        .expect("Error -- could not parse max_hits_to_report as int64")
        as usize;
    let intersect_level = match intersect_level {
        0 => align::IntersectLevel::NoIntersect,
        1 => align::IntersectLevel::IntersectWithFallback,
        2 => align::IntersectLevel::ForceIntersect,
        _ => panic!("Error -- invalid intersect level in config file. Please choose intersect level 0, 1, or 2.")
    };
    let group_on = aligner_config_obj["group_on"]
        .as_str()
        .expect("Error -- could not parse group_on as string")
        .to_string();
    let trim_target_length = aligner_config_obj["trim_target_length"]
        .as_i64()
        .expect("Error -- could not parse trim_target_length as usize")
        as usize;
    let trim_strictness = aligner_config_obj["trim_strictness"]
        .as_f64()
        .expect("Error -- could not parse trim_strictness as f64")
        as f64;


    // Get reference library from the second JSON object in the file
    let reference_obj = &v[1];
    let headers = to_string_vec(&reference_obj["headers"], "headers");
    let columns = &reference_obj["columns"];
    let sequence_name_idx =
        get_column_index(&headers, "sequence_name").expect("Could not find header sequence_name");

    // If there is no defined column for grouping, each feature is its own "feature family." Otherwise, use the defined column
    let group_on = if group_on == "" {
        sequence_name_idx
    } else {
        unwrap!(
            get_column_index(&headers, &group_on),
            "Error -- could not find column for group_on {}",
            &group_on
        )
    };

    let sequence_idx =
        get_column_index(&headers, "sequence").expect("Error -- could not find sequences column");
    let columns = columns
        .as_array()
        .expect("Error -- could not parse columns as array");
    let columns: Vec<Vec<String>> = columns
        .iter()
        .map(|column| to_string_vec(column, "column"))
        .collect();


    // Wrap all the parsed values into structs
    let align_config = align::AlignFilterConfig {
        reference_genome_size: columns[sequence_name_idx].len(),
        score_percent,
        score_threshold,
        num_mismatches,
        discard_nonzero_mismatch: false,
        discard_multiple_matches,
        score_filter: score_filter as i32,
        require_valid_pair,
        discard_multi_hits,
        intersect_level,
        max_hits_to_report,
        strand_filter,
        trim_target_length,
        trim_strictness
    };

    /* For each feature in the reference genome, we add a reverse-complemented version. This is in order to compute alignment orientation,
     * which informs several layers of filtration. */
    let mut new_columns = Vec::new();
    let num_rows = columns[0].len();

    for row_idx in 0..num_rows {
        let mut row = Vec::new();
        let mut revcomp_row = Vec::new();

        for (col_idx, col) in columns.iter().enumerate() {
            let mut value = col[row_idx].clone();

            if col_idx == sequence_idx {
                value = value.replace('U', "T").replace('u', "t");
            }
    
            row.push(value.clone());
            revcomp_row.push(value.clone());
        }

        revcomp_row[sequence_name_idx] = revcomp_row[sequence_name_idx].clone() + SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev";
        revcomp_row[sequence_idx] = revcomp(&revcomp_row[sequence_idx]);

        new_columns.push(row);
        new_columns.push(revcomp_row);
    }

    // Transpose new_columns back to the original column format
    let mut final_columns: Vec<Vec<String>> = vec![vec![]; columns.len()];
    for row in new_columns {
        for (i, val) in row.iter().enumerate() {
            final_columns[i].push(val.clone());
        }
    }

    let reference_metadata = Reference {
        group_on,
        headers,
        columns: final_columns,
        sequence_name_idx,
        sequence_idx,
    };

    sanity_check_align_config(&align_config);

    (align_config, reference_metadata)
}

// Given a column header, find the index of the corresponding column if it exists
fn get_column_index(headers: &[String], search_header: &str) -> Option<usize> {
    for (i, header) in headers.iter().enumerate() {
        if header == search_header {
            return Some(i);
        }
    }

    None
}

// Convert a given serde_json value into a string array if possible, and crash otherwise
fn to_string_vec(v: &Value, array_name: &str) -> Vec<String> {
    let result: Vec<String> = unwrap!(
        v.as_array(),
        "Error -- could not parse {} as array",
        array_name
    )
    .iter()
    .map(|string| {
        unwrap!(
            string.as_str(),
            "Error -- could not parse {} element \"{}\" as a string",
            array_name,
            string
        )
        .to_string()
    })
    .collect();

    result
}

pub fn sanity_check_align_config(
    align_config: &align::AlignFilterConfig
) {
    // score_percent should be between 0 and 1
    if !(0.0..=1.0).contains(&align_config.score_percent) {
        panic!("Error -- score_percent must be between 0 and 1");
    }

    // score_filter should be positive
    if align_config.score_filter < 0 {
        panic!("Error -- score_filter must be positive");
    }

    // trim_strictness should be between 0 and 1
    if !(0.0..=1.0).contains(&align_config.trim_strictness) {
        panic!("Error -- trim_strictness must be between 0 and 1");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    // Test for valid conversion of JSON array to Vec<String>
    #[test]
    fn test_to_string_vec_valid_conversion() {
        let json_value = json!(["apple", "banana", "cherry"]);
        let expected = vec!["apple".to_string(), "banana".to_string(), "cherry".to_string()];
        assert_eq!(to_string_vec(&json_value, "test array"), expected);
    }

    // Test handling when array contains non-string elements
    #[test]
    #[should_panic(expected = "Error -- could not parse test array element")]
    fn test_to_string_vec_invalid_elements() {
        let json_value = json!([true, "banana", 42]);
        to_string_vec(&json_value, "test array");
    }

    // Test with an empty JSON array
    #[test]
    fn test_to_string_vec_empty_array() {
        let json_value = json!([]);
        let expected: Vec<String> = Vec::new();
        assert_eq!(to_string_vec(&json_value, "empty array"), expected);
    }

    // Test for error handling when input is not an array
    #[test]
    #[should_panic(expected = "Error -- could not parse not an array as array")]
    fn test_to_string_vec_not_an_array() {
        let json_value = json!({"key": "value"});
        to_string_vec(&json_value, "not an array");
    }

    // Test for successful retrieval of header index
    #[test]
    fn test_get_column_index_found() {
        let headers = vec!["id".to_string(), "name".to_string(), "age".to_string()];
        let search_header = "name";
        let expected = Some(1);  // index of "name" in the headers
        assert_eq!(get_column_index(&headers, search_header), expected);
    }

    // Test for header not found in the list
    #[test]
    fn test_get_column_index_not_found() {
        let headers = vec!["id".to_string(), "name".to_string(), "age".to_string()];
        let search_header = "email";  // "email" is not in the headers
        let expected = None;
        assert_eq!(get_column_index(&headers, search_header), expected);
    }

    // Test with an empty headers list
    #[test]
    fn test_get_column_index_empty_list() {
        let headers: Vec<String> = Vec::new();  // empty list
        let search_header = "anything";
        let expected = None;
        assert_eq!(get_column_index(&headers, search_header), expected);
    }

    // Test for duplicate headers in the list
    #[test]
    fn test_get_column_index_duplicate_headers() {
        let headers = vec!["id".to_string(), "name".to_string(), "name".to_string(), "age".to_string()];
        let search_header = "name";
        let expected = Some(1);  // should return the index of the first occurrence of "name"
        assert_eq!(get_column_index(&headers, search_header), expected);
    }

    #[test]
    fn test_get_reference_library_valid_json() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-correct.json");
        let strand_filter = LibraryChemistry::None;
        let (align_config, reference_metadata) = get_reference_library(path, strand_filter);

        assert_eq!(align_config.score_percent, 0.85);
        assert_eq!(align_config.score_filter, 200);
        assert_eq!(align_config.score_threshold, 300);
        assert_eq!(align_config.num_mismatches, 2);
        assert_eq!(align_config.discard_multiple_matches, true);
        assert_eq!(align_config.require_valid_pair, false);
        assert_eq!(align_config.discard_multi_hits, 1);
        assert_eq!(align_config.intersect_level, align::IntersectLevel::IntersectWithFallback);
        assert_eq!(align_config.max_hits_to_report, 10);
        assert_eq!(align_config.trim_target_length, 40);
        assert_eq!(align_config.trim_strictness, 0.9);
        assert_eq!(reference_metadata.group_on, 1);
        assert_eq!(reference_metadata.headers, vec!["id".to_string(), "feature_id".to_string(), "sequence_name".to_string(), "sequence".to_string()]);
        assert_eq!(reference_metadata.columns[0], vec!["1".to_string(), "1".to_string(), "2".to_string(), "2".to_string()]);
        assert_eq!(reference_metadata.columns[1], vec!["fid1".to_string(), "fid1".to_string(), "fid2".to_string(), "fid2".to_string()]);
        assert_eq!(reference_metadata.columns[2], vec!["seq_name1".to_string(), "seq_name1§rev".to_string(), "seq_name2".to_string(), "seq_name2§rev".to_string()]);
        assert_eq!(reference_metadata.columns[3], vec!["ATGC".to_string(), "GCAT".to_string(), "CGTA".to_string(), "TACG".to_string()]);
        assert_eq!(reference_metadata.sequence_name_idx, 2);
        assert_eq!(reference_metadata.sequence_idx, 3);
    }

    // Test that missing fields fail to parse
    #[test]
    #[should_panic(expected = "Error -- could not parse score_percent as f64")]
    fn test_get_reference_library_missing_fields() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-missing-fields.json");
        let strand_filter = LibraryChemistry::None; 
        get_reference_library(path, strand_filter);
    }

    // Test that broken types fail to parse
    #[test]
    #[should_panic(expected = "Error -- could not parse score_percent as f64")]
    fn test_get_reference_library_incorrect_field_types() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-types-broken.json");
        let strand_filter = LibraryChemistry::None; 
        get_reference_library(path, strand_filter);
    }

    // Test that broken file format fails to parse
    #[test]
    #[should_panic(expected = "Error -- could not parse reference library JSON")]
    fn test_get_reference_library_corrupted_json() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-broken-format.json");
        let strand_filter = LibraryChemistry::None;
        get_reference_library(path, strand_filter);
    }

    #[test]
    #[should_panic(expected = "Error -- score_percent must be between 0 and 1")]
    fn test_invalid_score_percent() {
        let align_config = align::AlignFilterConfig {
            score_percent: 1.5, // invalid
            score_threshold: 100,
            num_mismatches: 2,
            score_filter: 50,
            discard_multi_hits: 1,
            max_hits_to_report: 10,
            trim_target_length: 40,
            trim_strictness: 0.9,
            discard_multiple_matches: true,
            require_valid_pair: true,
            intersect_level: align::IntersectLevel::NoIntersect,
            strand_filter: LibraryChemistry::None,
            discard_nonzero_mismatch: true,
            reference_genome_size: 1,
        };

        sanity_check_align_config(&align_config);
    }

    #[test]
    #[should_panic(expected = "Error -- score_filter must be positive")]
    fn test_negative_score_filter() {
        let align_config = align::AlignFilterConfig {
            score_percent: 0.9,
            score_threshold: 100,
            num_mismatches: 2,
            score_filter: -10, // invalid
            discard_multi_hits: 1,
            max_hits_to_report: 10,
            trim_target_length: 40,
            trim_strictness: 0.9,
            discard_multiple_matches: true,
            require_valid_pair: true,
            intersect_level: align::IntersectLevel::NoIntersect,
            strand_filter: LibraryChemistry::None,
            discard_nonzero_mismatch: true,
            reference_genome_size: 1,
        };

        sanity_check_align_config(&align_config);
    }

    #[test]
    #[should_panic(expected = "Error -- trim_strictness must be between 0 and 1")]
    fn test_invalid_trim_strictness() {
        let align_config = align::AlignFilterConfig {
            score_percent: 0.9,
            score_threshold: 100,
            num_mismatches: 2,
            score_filter: 50,
            discard_multi_hits: 1,
            max_hits_to_report: 10,
            trim_target_length: 40,
            trim_strictness: 1.5, // invalid
            discard_multiple_matches: true,
            require_valid_pair: true,
            intersect_level: align::IntersectLevel::NoIntersect,
            strand_filter: LibraryChemistry::None,
            discard_nonzero_mismatch: true,
            reference_genome_size: 1,
        };

        sanity_check_align_config(&align_config);
    }

    #[test]
    fn test_valid_config() {
        let align_config = align::AlignFilterConfig {
            score_percent: 0.85,
            score_threshold: 100,
            num_mismatches: 2,
            score_filter: 50,
            discard_multi_hits: 1,
            max_hits_to_report: 10,
            trim_target_length: 40,
            trim_strictness: 0.9,
            discard_multiple_matches: true,
            require_valid_pair: true,
            intersect_level: align::IntersectLevel::NoIntersect,
            strand_filter: LibraryChemistry::None,
            discard_nonzero_mismatch: true,
            reference_genome_size: 1,
        };

        sanity_check_align_config(&align_config);  // Should pass without panic
    }

    #[test]
    fn test_get_reference_library_rna_to_dna_conversion() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-rna.json");
        let strand_filter = LibraryChemistry::None;
        let (_align_config, reference_metadata) = get_reference_library(path, strand_filter);

        assert_eq!(reference_metadata.columns[3][0], "ATGCTT".to_string());
        assert_eq!(reference_metadata.columns[3][1], "AAGCAT".to_string());
        assert_eq!(reference_metadata.columns[3][2], "tTgcAT".to_string());
        assert_eq!(reference_metadata.columns[3][3], "ATgcAa".to_string());
    }

    #[test]
    fn test_get_reference_library_mixed_case_rna_to_dna_conversion() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-mixed-case-rna.json");
        let strand_filter = LibraryChemistry::None;
        let (_align_config, reference_metadata) = get_reference_library(path, strand_filter);

        assert_eq!(reference_metadata.columns[3][0], "atGcTt".to_string());
        assert_eq!(reference_metadata.columns[3][1], "aAgCat".to_string());
        assert_eq!(reference_metadata.columns[3][2], "TtgCAt".to_string());
        assert_eq!(reference_metadata.columns[3][3], "aTGcaA".to_string());
    }

    #[test]
    fn test_get_reference_library_sequence_without_rna_bases() {
        let path = Path::new("tests/test-sequences/libraries/reference-library-no-rna-bases.json");
        let strand_filter = LibraryChemistry::None;
        let (_align_config, reference_metadata) = get_reference_library(path, strand_filter);

        assert_eq!(reference_metadata.columns[3][0], "ATGCGT".to_string());
        assert_eq!(reference_metadata.columns[3][1], "ACGCAT".to_string());
        assert_eq!(reference_metadata.columns[3][2], "CGTACG".to_string());
        assert_eq!(reference_metadata.columns[3][3], "CGTACG".to_string());
    }
}