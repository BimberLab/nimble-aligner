use std::path::Path;
use std::fs::read_to_string;
use serde_json::Value;
use crate::align;


// Reads the reference library and returns it in the iterator format consumed by the aligner library
pub fn get_reference_library(path: &Path) -> align::AlignFilterConfig {
  let raw_json_string = read_to_string(path).expect(
    "Error -- could not read reference library for filtration"
  );

  let v: Value = serde_json::from_str(&raw_json_string).expect(
    "Error -- could not parse reference library JSON"
  );

  let config = &v[0];
  let reference = &v[1];

  let score_threshold = config["score_threshold"].as_i64().expect("Error -- could not parse score_threshold as int64") as usize;
  let percent_threshold = config["percent_threshold"].as_f64().expect("Error -- could not parse percent_threshold as float64");
  let num_mismatches = config["num_mismatches"].as_i64().expect("Error -- could not parse num_mismatches as int64") as usize;
  let discard_multiple_matches = config["discard_multiple_matches"].as_bool().expect("Error -- could not parse discard_multiple_Mismatches as boolean");
  let group_on = to_string_vec(&config["group_on"], "group_on");
  let headers = to_string_vec(&reference[0], "headers");
  let columns = &reference[1];

  align::AlignFilterConfig {
    reference_genome_size: 1209,
    score_threshold,
    num_mismatches,
    discard_differing_read_pairs: false,
    discard_nonzero_mismatch: false,
    discard_multiple_matches
  }
}


fn get_column_index(headers: Vec<String>, search_header: String) -> Option<usize> {
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