use std::path::PathBuf;
use tempfile::NamedTempFile;

extern crate nimble;
extern crate debruijn;
extern crate debruijn_mapping;

#[path = "./utils.rs"]
mod utils;

fn run_pipeline_on_sample(filename: &str) {
    let lib_filename = "strandedness.json";

    let (_, reference_index, reference_metadata, align_config) =
        utils::get_data(filename, lib_filename, nimble::align::LibraryChemistry::Unstranded);

    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test-sequences");
    data_path.push("reads");
    data_path.push(filename);

    let sequence_path = data_path
        .into_os_string()
        .into_string()
        .expect("Failed to convert path to string");

    let outfile = NamedTempFile::new().unwrap();
    let outfile_path = outfile.path().to_string_lossy().to_string();

    nimble::process::bam::process(
        vec![sequence_path],
        vec![reference_index],
        vec![reference_metadata],
        vec![align_config],
        vec![outfile_path],
        2,
        false,
    );
}

#[test]
fn test_pipeline_sample_bam() {
    run_pipeline_on_sample("sample.bam");
}

#[test]
fn test_pipeline_sample_no_r1_bam() {
    run_pipeline_on_sample("sample_no_r1.bam");
}