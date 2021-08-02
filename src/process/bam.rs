use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::bam;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;

pub fn process(
    input_files: Vec<&str>,
    reference_index: (PseudoAligner, PseudoAligner),
    reference_metadata: ReferenceMetadata,
    align_config: AlignFilterConfig,
    output_path: &str,
) {
    let mut reader = bam::UMIReader::new(input_files[0]);

    loop {
        let (final_umi, sequences, reverse_sequences) = reader.next();

        // Perform alignment and filtration using the score package
        let results = score(
            sequences,
            Some(reverse_sequences),
            reference_index,
            &reference_metadata,
            align_config,
        );

        if final_umi {
            return;
        };
    }
}
