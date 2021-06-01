use crate::align;
use serde_json::Value;
use std::fs::read_to_string;
use std::path::Path;
use unwrap::unwrap;

#[derive(Debug)]
pub struct ReferenceMetadata {
    pub group_on: usize,
    pub headers: Vec<String>,
    pub columns: Vec<Vec<String>>,
    pub sequence_name_idx: usize,
    pub sequence_idx: usize,
    pub data_type: String,
}

// Parses a .json that contains a reference library. Returns a tuple of the library's config information and the library data
pub fn get_reference_library(path: &Path) -> (align::AlignFilterConfig, ReferenceMetadata) {
    // Parse raw JSON to serde_json value
    let raw_json_string = read_to_string(path).expect("Error -- could not read reference library");

    let v: Value = serde_json::from_str(&raw_json_string)
        .expect("Error -- could not parse reference library JSON");

    // Get aligner configuration from the first JSON object in the file
    let config_obj = &v[0];
    let score_threshold = config_obj["score_threshold"]
        .as_i64()
        .expect("Error -- could not parse score_threshold as int64")
        as usize;
    let score_filter = config_obj["score_filter"]
        .as_i64()
        .expect("Error -- could not parse percent_threshold as int64");
    let num_mismatches = config_obj["num_mismatches"]
        .as_i64()
        .expect("Error -- could not parse num_mismatches as int64")
        as usize;
    let discard_multiple_matches = config_obj["discard_multiple_matches"]
        .as_bool()
        .expect("Error -- could not parse discard_multiple_mismatches as boolean");
    let require_valid_pair = config_obj["require_valid_pair"]
        .as_bool()
        .expect("Error -- could not parse require_valid_pair as boolean");
    let discard_multi_hits = config_obj["discard_multi_hits"]
        .as_i64()
        .expect("Error -- could not parse discard_multi_hits as int64")
        as usize;
    let intersect_level = config_obj["intersect_level"]
        .as_i64()
        .expect("Error -- could not parse intersect_level as int64");
    let intersect_level = match intersect_level {
    0 => align::IntersectLevel::NoIntersect,
    1 => align::IntersectLevel::IntersectWithFallback,
    2 => align::IntersectLevel::ForceIntersect,
    _ => panic!("Error -- invalid intersect level in config file. Please choose intersect level 0, 1, or 2.")
  };

    let group_on = config_obj["group_on"]
        .as_str()
        .expect("Error -- could not parse group_on as string")
        .to_string();

    // Get reference library metadata from the second JSON object in the file
    let reference = &v[1];
    let headers = to_string_vec(&reference["headers"], "headers");
    let columns = &reference["columns"];
    let sequence_name_idx =
        get_column_index(&headers, "sequence_name").expect("Could not find header sequence_name");
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

    // Parse columns into a matrix of strings
    let columns = columns
        .as_array()
        .expect("Error -- could not parse columns as array");
    let columns: Vec<Vec<String>> = columns
        .iter()
        .map(|column| to_string_vec(column, "column"))
        .collect();
    let data_type = config_obj["data_type"]
        .as_str()
        .expect("Error -- could not parse data_type as string")
        .to_string();

    let align_config = align::AlignFilterConfig {
        reference_genome_size: columns[sequence_name_idx].len(),
        score_threshold,
        num_mismatches,
        discard_nonzero_mismatch: false,
        discard_multiple_matches,
        score_filter: score_filter as i32,
        require_valid_pair,
        discard_multi_hits,
        intersect_level,
    };

    let reference_metadata = ReferenceMetadata {
        group_on,
        headers,
        columns,
        sequence_name_idx,
        sequence_idx,
        data_type
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
