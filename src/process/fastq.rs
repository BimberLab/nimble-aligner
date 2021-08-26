use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::fastq::get_error_checked_fastq_readers;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::write_to_tsv;

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
) {
    /* Get error-checked iterators to the sequences that will be aligned to the reference from the
     * sequence genome file(s) */
    let sequences = get_error_checked_fastq_readers(input_files[0]);

    // Only get reverse sequences if a reverse sequence file is provided
    let reverse_sequences = if input_files.len() > 1 {
        println!("Reading reverse sequences");
        Some(get_error_checked_fastq_readers(input_files[1]))
    } else {
        None
    };

    println!("Pseudo-aligning reads to reference index");

    // Perform alignment and filtration using the score package
    let results = score(
        sequences,
        reverse_sequences,
        reference_index,
        reference_metadata,
        align_config,
    );

    println!("Writing results to file");

    write_to_tsv(results, None, true, output_path);

    print!("Output results written to output path");
}
