use array_tool::vec::Intersect;
use debruijn::dna_string::DnaString;
use std::io::Error;
use std::collections::HashMap;

use crate::align::{AlignFilterConfig, AlignDebugInfo, PseudoAligner};
use crate::parse::bam;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::{write_to_tsv, write_debug_info, write_read_list, PseudoalignerData, BamSpecificAlignMetadata, revcomp};

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
    debug_file: Option<String>,
    alignment_file: Option<String>
) -> Vec<(Vec<String>, i32)> {
    let mut reader = bam::UMIReader::new(input_files[0]);
    let mut has_aligned = false;
    let mut score_map: HashMap<(Vec<String>, String), i32> = HashMap::new();

    let mut cell_barcodes: Vec<String> = Vec::new();

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

    let mut alignment_metadata: PseudoalignerData = PseudoalignerData {
        reference_names: Vec::new(),
        read_umi_name: Vec::new(),
        barcode_sample_name: Vec::new(),
        score: Vec::new(),
        pair: Vec::new(),
        sequence: Vec::new()
    };

    let mut bam_specific_alignment_metadata: BamSpecificAlignMetadata = BamSpecificAlignMetadata {
        mapq: Vec::new(),
        orientation: Vec::new()
    };

    let mut debug_info: AlignDebugInfo = Default::default();

    if owned_debug_file != "".to_owned() {
        debug_info.debug_file = owned_debug_file.clone();
    };

    loop {
        let final_umi = reader.next();

        if final_umi && has_aligned {
            println!("Finished aligning. Writing output files.");

            if owned_debug_file != "".to_owned() {
                println!("Writing debug table.");
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

            println!("Writing results.");
            write_to_tsv(&results, Some(cell_barcodes), false, output_path);

            return results;
        };
        has_aligned = true;

        let mut current_umi_group = reader.current_umi_group.clone();
        let current_metadata_table = metadata_to_sequence_hashmap(reader.current_metadata_group.clone(), current_umi_group.clone());

        let (extra_read, extra_metadata) = if current_umi_group.len() % 2 != 0 {
            (current_umi_group.pop(), reader.current_metadata_group.pop())
        } else {
            (None, None)
        };

        let (mut s, mut res) = if owned_debug_file.clone() != "" {
            get_score(
                &current_umi_group,
                reference_index,
                reference_metadata,
                align_config,
                Some(&mut debug_info),
                &reader.current_metadata_group.clone().into_iter().map(|(_, _, _, r)| r).collect::<Vec<bool>>())
        } else {
            get_score(
                &current_umi_group,
                reference_index,
                reference_metadata,
                align_config,
                None,
                &reader.current_metadata_group.clone().into_iter().map(|(_, _, _, r)| r).collect::<Vec<bool>>())
        };

        let (mut s_extra, mut res_extra) = match extra_read {
            Some(read) => {

                let (read_f, read_r) = if extra_metadata.unwrap().3 {
                    (vec![Ok(check_reverse_comp((&read.clone(), &true)))], vec![Ok(check_reverse_comp((&read.clone(), &true)))])
                } else {
                    (vec![Ok(read.clone())], vec![Ok(read.clone())])
                };

                score(
                    (Box::new(read_f.into_iter()), Box::new(read_r.into_iter())),
                    None,
                    reference_index,
                    &reference_metadata,
                    align_config,
                    None
                )
            },
            None => (Vec::new(), Vec::new())
        };

        s.append(&mut s_extra);
        res.append(&mut res_extra);
        
        bam_specific_alignment_metadata.mapq.append(&mut get_mapq_scores(get_sequence_list_from_metadata(&res), &current_metadata_table));
        bam_specific_alignment_metadata.orientation.append(&mut get_orientation(get_sequence_list_from_metadata(&res), &current_metadata_table));
        alignment_metadata.reference_names.append(&mut res.clone().into_iter().map(|(group, _, _)| group).collect::<Vec<Vec<String>>>());
        alignment_metadata.sequence.append(&mut res.clone().into_iter().map(|(_, seq, _)| seq).collect::<Vec<String>>());
        alignment_metadata.score.append(&mut res.clone().into_iter().map(|(_, _, score)| score).collect::<Vec<usize>>());
        alignment_metadata.barcode_sample_name.append(&mut res.clone().into_iter().map(|_| (&reader).current_cell_barcode.clone()).collect::<Vec<String>>());
        alignment_metadata.read_umi_name.append(&mut res.clone().into_iter().map(|_| (&reader).current_umi.clone()).collect::<Vec<String>>());
        alignment_metadata.pair.append(&mut get_pair(get_sequence_list_from_metadata(&res), &current_metadata_table));

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

        if owned_alignment_file != "".to_owned() {
            write_read_list(&alignment_metadata, Some(&bam_specific_alignment_metadata), &owned_alignment_file);

            bam_specific_alignment_metadata.mapq.clear();
            bam_specific_alignment_metadata.orientation.clear();
            alignment_metadata.reference_names.clear();
            alignment_metadata.sequence.clear();
            alignment_metadata.score.clear();
            alignment_metadata.barcode_sample_name.clear();
            alignment_metadata.read_umi_name.clear();
            alignment_metadata.pair.clear();
        }
    }
}

fn metadata_to_sequence_hashmap(metadata: Vec<(u8, String, String, bool)>, sequences: Vec<DnaString>) -> HashMap<String, (u8, String, String, bool)> {
    let mut ret: HashMap<String, (u8, String, String, bool)> = HashMap::new();

    for (i, seq) in sequences.iter().enumerate() {
        ret.insert(seq.to_string(), metadata[i].clone());
    }

    ret
}

fn get_mapq_scores(sequences: Vec<String>, metadata_table: &HashMap<String, (u8, String, String, bool)>) -> Vec<u8> {
    let mut ret: Vec<u8> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.0),
            None => if let Some(v) = metadata_table.get(&rc.to_string()) {
                ret.push(v.0)
            } else {
                ret.push(0)
            }
        }
    }

    ret
}

fn get_orientation(sequences: Vec<String>, metadata_table: &HashMap<String, (u8, String, String, bool)>) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.1.clone()),
            None => if let Some(v) = metadata_table.get(&rc.to_string()) {
                ret.push(v.1.clone())
            } else {
                ret.push(String::new())
            }
        }
    }

    ret
}

fn get_pair(sequences: Vec<String>, metadata_table: &HashMap<String, (u8, String, String, bool)>) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.2.clone()),
            None => if let Some(v) = metadata_table.get(&rc.to_string()) {
                ret.push(v.2.clone())
            } else {
                ret.push(String::new())
            }
        }
    }

    ret
}

fn get_sequence_list_from_metadata(metadata: &Vec<(Vec<String>, String, usize)>) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for (_, seq, _) in metadata.iter() {
        ret.push(seq.to_string());
    }

    ret
}

fn get_score<'a>(
    current_umi_group: &'a Vec<DnaString>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
    reverse_comp_read: &'a Vec<bool>
) -> (Vec<(Vec<String>, i32)>, Vec<(Vec<String>, String, usize)>) {
    let sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group.iter().zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group.iter().zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let reverse_sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group.iter().zip(reverse_comp_read.iter())
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let reverse_sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group.iter().zip(reverse_comp_read.iter())
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
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

fn check_reverse_comp(rec: (&DnaString, &bool)) -> DnaString {
    let (seq, reverse_comp) = rec;

    if *reverse_comp {
        DnaString::from_dna_string(&revcomp(&seq.to_string())).to_owned()
    } else {
        seq.to_owned()
    }
}