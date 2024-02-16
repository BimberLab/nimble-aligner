use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::fastq::get_error_checked_fastq_readers;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::write_to_tsv;

pub fn process(
    input_files: Vec<String>,
    reference_indices: Vec<(PseudoAligner, PseudoAligner)>,
    reference_metadata: Vec<ReferenceMetadata>,
    align_configs: Vec<AlignFilterConfig>,
    output_paths: Vec<String>,
    _num_cores: usize,
) {
    for (i, index) in reference_indices.into_iter().enumerate() {
        let (results, _alignment_metadata, _) = score(
            get_error_checked_fastq_readers(input_files[0].clone()),
            if input_files.len() > 1 { Some(get_error_checked_fastq_readers(input_files[1].clone())) } else { None },
            &Vec::new(),
            &index,
            &reference_metadata[i],
            &align_configs[i],
            None,
        );

        write_to_tsv(
            &results.into_iter().map(|(features, (score, _, _))| (features, score)).collect(),
            output_paths[i].clone()
        );
    }
}
