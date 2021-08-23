use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::bam;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use debruijn::dna_string::DnaString;
use std::io::Error;

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
) {
    let mut reader = bam::UMIReader::new(input_files[0]);

    loop {
        let final_umi = reader.next();

        let current_umi_group = reader.current_umi_group.clone();

        let scores = get_score(
            &current_umi_group,
            reference_index,
            reference_metadata,
            align_config,
        );

        if final_umi {
            return;
        };
    }
}

fn get_score<'a>(
    current_umi_group: &'a Vec<DnaString>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
) -> Vec<(Vec<String>, i32)> {
    let sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .step_by(2)
            .map(|rec| Ok(rec.to_owned())),
    );

    let sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .step_by(2)
            .map(|rec| Ok(rec.to_owned())),
    );

    let reverse_sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(rec.to_owned())),
    );

    let reverse_sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(rec.to_owned())),
    );

    // Perform alignment and filtration using the score package
    score(
        (sequences, sequences_clone),
        Some((reverse_sequences, reverse_sequences_clone)),
        reference_index,
        &reference_metadata,
        align_config,
    )
}
