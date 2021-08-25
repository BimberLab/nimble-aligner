use std::io::Error;
use debruijn::dna_string::DnaString;
use array_tool::vec::Intersect;

use crate::align::{AlignFilterConfig, PseudoAligner};
use crate::parse::bam;
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
    let mut reader = bam::UMIReader::new(input_files[0]);

    loop {
        let final_umi = reader.next();

        if final_umi {
            return;
        };

        let current_umi_group = reader.current_umi_group.clone();

        let s = get_score(
            &current_umi_group,
            reference_index,
            reference_metadata,
            align_config,
        );
        
        if s.len() == 0 {
            continue;
        }

        let mut scores = s.iter();
        
        let first_score = scores.next().unwrap();
        let mut group = first_score.0.clone();
        let mut score = first_score.1;

        loop {
            let next = scores.next();

            if next.is_none() {
                break;
            }

            group = group.intersect(next.unwrap().0.clone());
            score += next.unwrap().1;
        }
        
        if group.len() > 0 {
            write_to_tsv(vec![(group, score)], output_path);
        }
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

    // TODO: Check on this logic. Currently, if the data isn't all read-pair, we throw out ALL of
    // the reverse pairs.
    let reverse_sequence_pair = if current_umi_group.len() % 2 == 0 {
        Some((reverse_sequences, reverse_sequences_clone))
    } else {
        None
    };

    // Perform alignment and filtration using the score package
    score(
        (sequences, sequences_clone),
        reverse_sequence_pair,
        reference_index,
        &reference_metadata,
        align_config,
    )
}
