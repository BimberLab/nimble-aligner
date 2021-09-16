use array_tool::vec::Intersect;
use debruijn::dna_string::DnaString;
use std::io::Error;

use crate::align::{AlignFilterConfig, AlignDebugInfo, PseudoAligner};
use crate::parse::bam;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::{write_to_tsv, write_debug_info};

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
    debug_file: Option<String>
) {
    let mut reader = bam::UMIReader::new(input_files[0]);
    let mut has_written_headers = false;
    
    let owned_debug_file = if debug_file.is_some() {
        debug_file.unwrap()
    } else {
        "".to_owned()
    };

    let mut debug_info: AlignDebugInfo = Default::default();

    if owned_debug_file != "".to_owned() {
        debug_info.debug_file = owned_debug_file.clone();
    };

    loop {
        let final_umi = reader.next();

        if final_umi {
            if owned_debug_file != "".to_owned() {
                write_debug_info(debug_info);
            }

            return;
        };

        let current_umi_group = reader.current_umi_group.clone();

        let s = if owned_debug_file.clone() != "" {
            get_score(
                &current_umi_group,
                reference_index,
                reference_metadata,
                align_config,
                Some(&mut debug_info))
        } else {
            get_score(
                &current_umi_group,
                reference_index,
                reference_metadata,
                align_config,
                None)
        };

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
            write_to_tsv(vec![(group, score)], Some(reader.current_cell_barcode.clone()), !has_written_headers, output_path);
            has_written_headers = true;
        }
    }
}

fn get_score<'a>(
    current_umi_group: &'a Vec<DnaString>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
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
    let (reverse_sequence_pair, debug_info) = if current_umi_group.len() % 2 == 0 {
        (Some((reverse_sequences, reverse_sequences_clone)), debug_info)
    } else {
        if debug_info.is_some() {
            let debug_info = debug_info.unwrap();
            debug_info.reverse_read_sets_discarded_noneven += 1;
            (None, Some(debug_info))
        } else {
            (None, None)
        }
    };

    // Perform alignment and filtration using the score package
    score(
        (sequences, sequences_clone),
        reverse_sequence_pair,
        reference_index,
        &reference_metadata,
        align_config,
        debug_info
    )
}
