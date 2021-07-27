use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::reference_library::ReferenceMetadata;

pub fn process(
    input_files: Vec<&str>,
    reference_index: (PseudoAligner, PseudoAligner),
    reference_metadata: ReferenceMetadata,
    align_config: AlignFilterConfig,
    output_path: &str,
) {
}
