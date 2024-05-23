use crate::align::{self, LibraryChemistry};
use serde_json::Value;
use std::fs::read_to_string;
use std::path::Path;
use unwrap::unwrap;

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
        strand_filter
    };

    let reference_metadata = Reference {
        group_on,
        headers,
        columns,
        sequence_name_idx,
        sequence_idx,
    };

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
        assert_eq!(reference_metadata.group_on, 1);
        assert_eq!(reference_metadata.headers, vec!["id".to_string(), "feature_id".to_string(), "sequence_name".to_string(), "sequence".to_string()]);
        assert_eq!(reference_metadata.columns[0], vec!["1".to_string(), "2".to_string()]);
        assert_eq!(reference_metadata.columns[1], vec!["fid1".to_string(), "fid2".to_string()]);
        assert_eq!(reference_metadata.columns[2], vec!["seq_name1".to_string(), "seq_name2".to_string()]);
        assert_eq!(reference_metadata.columns[3], vec!["ATGC".to_string(), "CGTA".to_string()]);
        assert_eq!(reference_metadata.sequence_name_idx, 2);
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
}