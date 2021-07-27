use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::reference_library::ReferenceMetadata;
use crate::parse::bam;

pub fn process(
    input_files: Vec<&str>,
    reference_index: (PseudoAligner, PseudoAligner),
    reference_metadata: ReferenceMetadata,
    align_config: AlignFilterConfig,
    output_path: &str,
) {
    let mut reader = bam::UMIReader::new(input_files[0]);

    loop {
    let final_umi = reader.next();

        
    }
}