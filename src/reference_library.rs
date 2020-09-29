use std::path::Path;
use std::fs::read_to_string;
use serde_json::Value;
use crate::align;


#[derive(Debug)]
pub struct ReferenceMetadata {
  pub group_on: Vec<usize>,
  pub headers: Vec<String>,
  pub columns: Vec<Vec<String>>,
  pub nt_sequence_idx: usize
}


// Reads the reference library and returns it in the iterator format consumed by the aligner library
pub fn get_reference_library(path: &Path) -> (align::AlignFilterConfig, ReferenceMetadata) {
  // Parse raw JSON to serde_json value
  let raw_json_string = read_to_string(path).expect(
    "Error -- could not read reference library for filtration"
  );

  let v: Value = serde_json::from_str(&raw_json_string).expect(
    "Error -- could not parse reference library JSON"
  );


  // Get aligner configuration from the first JSON object in the file
  let config_obj = &v[0];
  let score_threshold = config_obj["score_threshold"].as_i64().expect("Error -- could not parse score_threshold as int64") as usize;
  let percent_threshold = config_obj["percent_threshold"].as_f64().expect("Error -- could not parse percent_threshold as float64");
  let num_mismatches = config_obj["num_mismatches"].as_i64().expect("Error -- could not parse num_mismatches as int64") as usize;
  let discard_multiple_matches = config_obj["discard_multiple_matches"].as_bool().expect("Error -- could not parse discard_multiple_mismatches as boolean");
  let group_on = to_string_vec(&config_obj["group_on"], "group_on");

  let align_config = align::AlignFilterConfig {
    reference_genome_size: 1209,
    score_threshold,
    num_mismatches,
    discard_differing_read_pairs: false,
    discard_nonzero_mismatch: false,
    discard_multiple_matches
  };


  // Get reference library metadata from the second JSON object in the file
  let reference = &v[1];
  let headers = to_string_vec(&reference["headers"], "headers");
  let columns = &reference["columns"];
  let group_on: Vec<usize> = group_on.into_iter().map(|group| get_column_index(&headers, &group).expect(&format!("Error -- could not parse group {}", group))).collect();
  let nt_sequence_idx = get_column_index(&headers, "nt_sequence").expect("Could not find header nt_sequence");

  // Parse columns into a matrix of strings
  let columns = columns.as_array().expect("Error -- could not parse columns as array");
  let columns: Vec<Vec<String>> = columns.into_iter().map(|column| to_string_vec(column, "column")).collect();

  let reference_metadata = ReferenceMetadata {
    group_on,
    headers,
    columns,
    nt_sequence_idx 
  };
  
  (align_config, reference_metadata)
}


fn get_column_index(headers: &Vec<String>, search_header: &str) -> Option<usize> {
  for (i, header) in headers.into_iter().enumerate() {
    if header == search_header {
      return Some(i);
    }
  }

  None
}


fn to_string_vec(v: &Value, array_name: &str) -> Vec<String> {
  let result: Vec<String> = v.as_array().expect(&format!("Error -- could not parse {} as array", array_name)).into_iter().map(|string| {
    string.as_str().expect(&format!("Error -- could not parse {} element \"{}\" as a string", array_name, string)).to_string()
  }).collect();

  result
}