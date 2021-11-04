use array_tool::vec::Intersect;
use debruijn::dna_string::DnaString;
use std::io::Error;
use std::collections::HashMap;

use crate::align::{AlignFilterConfig, AlignDebugInfo, PseudoAligner};
use crate::parse::bam;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::{write_to_tsv, write_debug_info, filter_scores};

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
    debug_file: Option<String>
) {
    let mut reader = bam::UMIReader::new(input_files[0]);
    let mut score_map: HashMap<(Vec<String>, String), i32> = HashMap::new();
    let mut cell_barcodes: Vec<String> = Vec::new();

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

            let mut results = Vec::new();
            for (key, value) in score_map.into_iter() {
                let group = key.0;
                let cell_barcode = key.1;
                let score = value;

                results.push((group, score));
                cell_barcodes.push(cell_barcode);
            }

            write_to_tsv(filter_scores(results, &align_config.score_filter), Some(cell_barcodes), false, output_path);

            return;
        };

        let mut current_umi_group = reader.current_umi_group.clone();

        let extra_read = if current_umi_group.len() % 2 != 0 {
            current_umi_group.pop()
        } else {
            None
        };

        let mut s = if owned_debug_file.clone() != "" {
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

        let mut s_extra = match extra_read {
            Some(read) => {
                let read_f = vec![Ok(read.clone())];
                let read_r = vec![Ok(read.clone())];
                score(
                    (Box::new(read_f.into_iter()), Box::new(read_r.into_iter())),
                    None,
                    reference_index,
                    &reference_metadata,
                    align_config,
                    None
                )
            },
            None => Vec::new()
        };

        s.append(&mut s_extra);

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
            let accessor = score_map.entry((group, reader.current_cell_barcode.clone())).or_insert(0);
            *accessor = *accessor + score;
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

    let reverse_sequence_pair = Some((reverse_sequences, reverse_sequences_clone));

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
