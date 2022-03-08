use crate::align::{AlignFilterConfig, AlignDebugInfo, PseudoAligner};
use crate::parse::fastq::get_error_checked_fastq_readers;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::{write_to_tsv, write_debug_info, filter_scores, write_read_list};

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
    debug_file: Option<String>,
    alignment_file: Option<String>
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

    let owned_debug_file = if debug_file.is_some() {
        debug_file.unwrap()
    } else {
        "".to_owned()
    };

    let owned_alignment_file = if alignment_file.is_some() {
        alignment_file.unwrap()
    } else {
        "".to_owned()
    };

    let mut debug_info: AlignDebugInfo = Default::default();

    if owned_debug_file != "".to_owned() {
        debug_info.debug_file = owned_debug_file.clone();
    };

    // Perform alignment and filtration using the score package
    let (results, alignment_metadata) = if owned_debug_file.clone() != "" {
        score(
            sequences,
            reverse_sequences,
            reference_index,
            reference_metadata,
            align_config,
            Some(&mut debug_info)
        )
    } else {
        score(
            sequences,
            reverse_sequences,
            reference_index,
            reference_metadata,
            align_config,
            None
        )
    };

    println!("Writing results to file");

    write_to_tsv(filter_scores(results, &align_config.score_filter), None, true, output_path);

    if owned_debug_file != "".to_owned() {
        write_debug_info(debug_info);
    }

    if owned_alignment_file != "".to_owned() {
        write_read_list(alignment_metadata, &owned_alignment_file);
    }

    print!("Output results written to output path");
}
