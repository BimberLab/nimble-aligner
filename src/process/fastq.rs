use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::fastq::get_error_checked_fastq_readers;
use crate::reference_library::Reference;
use crate::score::score;
use crate::utils::write_to_tsv;

pub fn process(
    input_files: Vec<String>,
    reference_indices: Vec<(PseudoAligner, PseudoAligner)>,
    references: Vec<Reference>,
    aligner_configs: Vec<AlignFilterConfig>,
    output_paths: Vec<String>
) {
    // For each reference, pass iterators of one or both fastq files to the scoring pipeline. Then, write the fastq results to one TSV per input.
    for (i, index) in reference_indices.into_iter().enumerate() {
        let (results, _alignment_metadata, _) = score(
            get_error_checked_fastq_readers(input_files[0].clone()),
            if input_files.len() > 1 { Some(get_error_checked_fastq_readers(input_files[1].clone())) } else { None },
            &Vec::new(),
            &index,
            &references[i],
            &aligner_configs[i],
            None,
        );

        write_to_tsv(
            &results.into_iter().map(|(features, (score, _, _))| (features, score)).collect(),
            output_paths[i].clone()
        );
    }
}
