use crate::align::{AlignDebugInfo, AlignFilterConfig, PseudoAligner};
use crate::parse::bam::UMIReader;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::{revcomp, write_debug_info, BamSpecificAlignMetadata, PseudoalignerData};
use array_tool::vec::Intersect;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::scores::blosum62;
use debruijn::dna_string::DnaString;
use rayon::prelude::*;
use rust_htslib::{bam, bam::header::HeaderRecord, bam::record::Aux, bam::record::Record};

use std::sync::{Arc, Mutex};

use std::collections::HashMap;
use std::io::Error;
use std::sync::mpsc::channel;

pub fn process(
    input_files: Vec<&str>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    output_path: &str,
    debug_file: Option<String>,
    alignment_file: Option<String>,
    num_cores: usize,
    reference_seqs: Vec<String>,
    reference_names: Vec<String>,
) -> Vec<(Vec<String>, i32)> {
    let mut reader = UMIReader::new(input_files[0], debug_file.is_none());
    let mut has_aligned = false;
    //let mut score_map: HashMap<(Vec<String>, String), i32> = HashMap::new();
    let mut header = bam::Header::new();
    let mut hd_record = HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", &"1.6");
    hd_record.push_tag(b"SO", &"unknown");
    header.push_record(&hd_record);
    let mut bam_writer = bam::Writer::from_path(&output_path, &header, bam::Format::Bam).unwrap();

    let cell_barcodes: Vec<String> = Vec::new();

    //let program_start = Instant::now();

    let owned_debug_file = if debug_file.is_some() {
        debug_file.unwrap()
    } else {
        "".to_owned()
    };

    let _owned_alignment_file = if alignment_file.is_some() {
        alignment_file.unwrap()
    } else {
        "".to_owned()
    };

    let mut alignment_metadata: PseudoalignerData = PseudoalignerData {
        reference_names: Vec::new(),
        read_umi_name: Vec::new(),
        barcode_sample_name: Vec::new(),
        score: Vec::new(),
        raw_score: Vec::new(),
        pair: Vec::new(),
        sequence: Vec::new(),
        strand_filter_reason: Vec::new(),
    };

    let mut bam_specific_alignment_metadata: BamSpecificAlignMetadata = BamSpecificAlignMetadata {
        mapq: Vec::new(),
        orientation: Vec::new(),
        hits: Vec::new(),
        quals: Vec::new(),
        qnames: Vec::new(),
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
                debug_info.number_cr_skipped = reader.number_cr_skipped;
                write_debug_info(debug_info);
            }

            return Vec::new();
        };
        has_aligned = true;

        //let start = Instant::now();
        let mut current_umi_group = reader.current_umi_group.clone();
        let current_metadata_table = metadata_to_sequence_hashmap(
            reader.current_metadata_group.clone(),
            current_umi_group.clone(),
        );

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
                &reader
                    .current_metadata_group
                    .clone()
                    .into_iter()
                    .map(|(_, _, _, r, _, _, _)| r)
                    .collect::<Vec<bool>>(),
            )
        } else {
            get_score(
                &current_umi_group,
                reference_index,
                reference_metadata,
                align_config,
                None,
                &reader
                    .current_metadata_group
                    .clone()
                    .into_iter()
                    .map(|(_, _, _, r, _, _, _)| r)
                    .collect::<Vec<bool>>(),
            )
        };

        //let duration = start.elapsed();
        //println!("total time to get score: {:?}", duration);

        //let start = Instant::now();
        let (mut s_extra, mut res_extra) = match extra_read {
            Some(read) => {
                let (read_f, read_r) = if extra_metadata.unwrap().3 {
                    (
                        vec![Ok(check_reverse_comp((&read.clone(), &true)))],
                        vec![Ok(check_reverse_comp((&read.clone(), &true)))],
                    )
                } else {
                    (vec![Ok(read.clone())], vec![Ok(read.clone())])
                };

                score(
                    (Box::new(read_f.into_iter()), Box::new(read_r.into_iter())),
                    None,
                    reference_index,
                    &reference_metadata,
                    align_config,
                    None,
                )
            }
            None => (Vec::new(), Vec::new()),
        };

        s.append(&mut s_extra);
        res.append(&mut res_extra);

        if s.len() == 0 {
            //println!("\n");
            continue;
        }

        bam_specific_alignment_metadata
            .mapq
            .append(&mut get_mapq_scores(
                get_sequence_list_from_metadata(&res),
                &current_metadata_table,
            ));
        bam_specific_alignment_metadata
            .orientation
            .append(&mut get_orientation(
                get_sequence_list_from_metadata(&res),
                &current_metadata_table,
            ));
        bam_specific_alignment_metadata.hits.append(&mut get_hits(
            get_sequence_list_from_metadata(&res),
            &current_metadata_table,
        ));

        bam_specific_alignment_metadata
            .qnames
            .append(&mut get_qnames(
                get_sequence_list_from_metadata(&res),
                &current_metadata_table,
            ));
        bam_specific_alignment_metadata.quals.append(&mut get_quals(
            get_sequence_list_from_metadata(&res),
            &current_metadata_table,
        ));
        alignment_metadata.reference_names.append(
            &mut res
                .clone()
                .into_iter()
                .map(|(group, _, _, _, _)| group)
                .collect::<Vec<Vec<String>>>(),
        );
        alignment_metadata.sequence.append(
            &mut res
                .clone()
                .into_iter()
                .map(|(_, seq, _, _, _)| seq)
                .collect::<Vec<String>>(),
        );
        alignment_metadata.score.append(
            &mut res
                .clone()
                .into_iter()
                .map(|(_, _, score, _, _)| score)
                .collect::<Vec<f64>>(),
        );
        alignment_metadata.raw_score.append(
            &mut res
                .clone()
                .into_iter()
                .map(|(_, _, _, raw_score, _)| raw_score)
                .collect::<Vec<usize>>(),
        );
        alignment_metadata.strand_filter_reason.append(
            &mut res
                .clone()
                .into_iter()
                .map(|(_, _, _, _, reason)| reason)
                .collect::<Vec<String>>(),
        );
        alignment_metadata.barcode_sample_name.append(
            &mut res
                .clone()
                .into_iter()
                .map(|_| (&reader).current_cell_barcode.clone())
                .collect::<Vec<String>>(),
        );
        alignment_metadata.read_umi_name.append(
            &mut res
                .clone()
                .into_iter()
                .map(|_| (&reader).current_umi.clone())
                .collect::<Vec<String>>(),
        );
        alignment_metadata.pair.append(&mut get_pair(
            get_sequence_list_from_metadata(&res),
            &current_metadata_table,
        ));

        // Perfom positional alignment
        let pos_alignments = positional_alignment(
            &reference_seqs,
            &alignment_metadata.sequence,
            &reference_names,
        );

        // Filter out all reads that didn't have positional alignments
        let mut bam_specific_alignment_metadata_filtered = BamSpecificAlignMetadata {
            mapq: Vec::new(),
            orientation: Vec::new(),
            hits: Vec::new(),
            quals: Vec::new(),
            qnames: Vec::new(),
        };

        let mut alignment_metadata_filtered = PseudoalignerData {
            reference_names: Vec::new(),
            read_umi_name: Vec::new(),
            barcode_sample_name: Vec::new(),
            score: Vec::new(),
            raw_score: Vec::new(),
            pair: Vec::new(),
            sequence: Vec::new(),
            strand_filter_reason: Vec::new(),
        };

        for (i, alns) in pos_alignments.iter().enumerate() {
            if !alns.is_empty() {
                alignment_metadata_filtered
                    .reference_names
                    .push(alignment_metadata.reference_names[i].clone());
                alignment_metadata_filtered
                    .read_umi_name
                    .push(alignment_metadata.read_umi_name[i].clone());
                alignment_metadata_filtered
                    .barcode_sample_name
                    .push(alignment_metadata.barcode_sample_name[i].clone());
                alignment_metadata_filtered
                    .score
                    .push(alignment_metadata.score[i]);
                alignment_metadata_filtered
                    .raw_score
                    .push(alignment_metadata.raw_score[i]);
                alignment_metadata_filtered
                    .pair
                    .push(alignment_metadata.pair[i].clone());
                alignment_metadata_filtered
                    .sequence
                    .push(alignment_metadata.sequence[i].clone());

                bam_specific_alignment_metadata_filtered
                    .mapq
                    .push(bam_specific_alignment_metadata.mapq[i]);
                bam_specific_alignment_metadata_filtered
                    .orientation
                    .push(bam_specific_alignment_metadata.orientation[i].clone());
                bam_specific_alignment_metadata_filtered
                    .hits
                    .push(bam_specific_alignment_metadata.hits[i].clone());
            }
        }

        let mut pos_alignments_filtered = Vec::new();

        for (i, alns) in pos_alignments.iter().enumerate() {
            if !alns.is_empty() {
                pos_alignments_filtered.push(alns.clone());
            }
        }

        write_to_bam(
            &mut bam_writer,
            &header,
            &reference_seqs,
            &alignment_metadata_filtered,
            &bam_specific_alignment_metadata_filtered,
            &pos_alignments_filtered,
            reader.current_cell_barcode.clone(),
        );

        /*write_read_list(
            &alignment_metadata,
            Some(&bam_specific_alignment_metadata),
            &owned_alignment_file,
        );*/

        bam_specific_alignment_metadata.mapq.clear();
        bam_specific_alignment_metadata.orientation.clear();
        bam_specific_alignment_metadata.hits.clear();
        alignment_metadata.reference_names.clear();
        alignment_metadata.sequence.clear();
        alignment_metadata.score.clear();
        alignment_metadata.raw_score.clear();
        alignment_metadata.strand_filter_reason.clear();
        alignment_metadata.barcode_sample_name.clear();
        alignment_metadata.read_umi_name.clear();
        alignment_metadata.pair.clear();

        //let duration = loop_start.elapsed();
        //println!("time to complete loop: {:?}", duration);

        //let duration = start.elapsed();
        //println!("time to perform misc tasks: {:?}", duration);

        //let start = Instant::now();
        /*let mut scores = s.iter();

        let first_score = scores.next().unwrap();
        let mut group = first_score.0.clone();
        let mut score = first_score.1;

        // A UMI-worth of reads is collapsed into a single feature set, determined by the intersection of
        // all of the matched features across all the UMI sequences, and their scores added
        loop {
            let next = scores.next();

            if next.is_none() {
                break;
            }

            group = group.intersect(next.unwrap().0.clone());
            score += next.unwrap().1;
        }

        if group.len() > 0 {
            let accessor = score_map
                .entry((group, reader.current_cell_barcode.clone()))
                .or_insert(0);
            *accessor = *accessor + score;
        }*/

        //let duration = start.elapsed();
        //println!("time to write score to table: {:?}\n\n", duration);
    }
}

fn metadata_to_sequence_hashmap(
    metadata: Vec<(u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
    sequences: Vec<DnaString>,
) -> HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)> {
    let mut ret: HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)> =
        HashMap::new();

    for (i, seq) in sequences.iter().enumerate() {
        ret.insert(seq.to_string(), metadata[i].clone());
    }

    ret
}

fn get_mapq_scores(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<u8> {
    let mut ret: Vec<u8> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.0),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.0)
                } else {
                    ret.push(0)
                }
            }
        }
    }

    ret
}

fn get_hits(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.4.clone()),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.4.clone())
                } else {
                    ret.push(String::new())
                }
            }
        }
    }

    ret
}

fn get_qnames(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<Vec<u8>> {
    let mut ret: Vec<Vec<u8>> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.5.clone()),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.5.clone())
                } else {
                    ret.push(Vec::new())
                }
            }
        }
    }

    ret
}

fn get_quals(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<Vec<u8>> {
    let mut ret: Vec<Vec<u8>> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.6.clone()),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.6.clone())
                } else {
                    ret.push(Vec::new())
                }
            }
        }
    }

    ret
}

fn get_orientation(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.1.clone()),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.1.clone())
                } else {
                    ret.push(String::new())
                }
            }
        }
    }

    ret
}

fn get_pair(
    sequences: Vec<String>,
    metadata_table: &HashMap<String, (u8, String, String, bool, String, Vec<u8>, Vec<u8>)>,
) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for seq in sequences {
        let rc = check_reverse_comp((&DnaString::from_dna_string(&seq), &true));
        let entry = metadata_table.get(&seq);

        match entry {
            Some(v) => ret.push(v.2.clone()),
            None => {
                if let Some(v) = metadata_table.get(&rc.to_string()) {
                    ret.push(v.2.clone())
                } else {
                    ret.push(String::new())
                }
            }
        }
    }

    ret
}

fn get_sequence_list_from_metadata(
    metadata: &Vec<(Vec<String>, String, f64, usize, String)>,
) -> Vec<String> {
    let mut ret: Vec<String> = Vec::new();

    for (_, seq, _, _, _) in metadata.iter() {
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
    reverse_comp_read: &'a Vec<bool>,
) -> (
    Vec<(Vec<String>, i32)>,
    Vec<(Vec<String>, String, f64, usize, String)>,
) {
    let sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let reverse_sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .zip(reverse_comp_read.iter())
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(check_reverse_comp(rec))),
    );

    let reverse_sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        current_umi_group
            .iter()
            .zip(reverse_comp_read.iter())
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
        debug_info,
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

fn positional_alignment(
    reference: &Vec<String>,
    sequence_list: &Vec<String>,
    reference_names: &Vec<String>,
) -> Vec<Vec<(Alignment, String)>> {
    let aligner = Aligner::with_capacity(
        reference[0].len(),
        sequence_list[0].len(),
        -5,
        -1,
        &blosum62,
    );
    let scoring = Scoring::new(-5, -1, &|a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

    let (sender, receiver) = channel();
    let sender = Arc::new(Mutex::new(sender));

    sequence_list
        .par_iter()
        .enumerate()
        .for_each(|(i, sequence)| {
            let mut aligner = Aligner::with_scoring(scoring.clone());
            let mut alignments = Vec::new();
            for (j, ref_seq) in reference.iter().enumerate() {
                let alignment = aligner.local(sequence.as_bytes(), ref_seq.as_bytes());
                alignments.push((alignment, reference_names[j].clone()));
            }
            alignments.sort_by(|a, b| b.0.score.cmp(&a.0.score));
            sender
                .lock()
                .unwrap()
                .send((i, alignments[0].clone()))
                .unwrap();
        });

    drop(sender);

    let mut ordered_results = Vec::with_capacity(sequence_list.len());
    ordered_results.resize(sequence_list.len(), Vec::new());

    for (i, alignment) in receiver.iter() {
        ordered_results[i].push(alignment);
    }

    ordered_results
}

fn write_to_bam(
    bam_writer: &mut bam::Writer,
    header: &bam::Header,
    reference_seqs: &[String],
    alignment_metadata_filtered: &PseudoalignerData,
    bam_specific_alignment_metadata_filtered: &BamSpecificAlignMetadata,
    pos_alignments_filtered: &[Vec<(Alignment, String)>],
    cb: String,
) {
    for (i, alns) in pos_alignments_filtered.iter().enumerate() {
        if !alns.is_empty() {
            let ref_name = alignment_metadata_filtered.reference_names[i][0].as_bytes();
            let sequence = alignment_metadata_filtered.sequence[i].as_bytes();

            let qname = bam_specific_alignment_metadata_filtered.qnames[i].as_slice();
            let qual = bam_specific_alignment_metadata_filtered.quals[i].as_slice();

            for (alignment, feature_name) in alns {
                let mut rec = Record::new();

                //header.tid(ref_name).unwrap();

                rec.set_pos(alignment.ystart as i64);
                rec.set_bin(0);

                // compute AS from the alignment score
                rec.push_aux(b"AS", Aux::I32(alignment.score as i32))
                    .unwrap();
                rec.push_aux(b"CB", Aux::String(&cb)).unwrap();
                rec.push_aux(b"XF", Aux::String(&feature_name)).unwrap();

                // set variable length data
                rec.set(
                    qname,
                    Some(&bam::record::CigarString::from_alignment(alignment, true)),
                    sequence,
                    qual,
                );

                rec.push_aux(b"NH", Aux::U8(0)).unwrap();
                rec.push_aux(b"HI", Aux::U8(0)).unwrap();
                rec.push_aux(b"uT", Aux::Char(b'P')).unwrap();

                let mut nm_count = 0;
                for operation in alignment.operations.iter() {
                    match operation {
                        bio::alignment::AlignmentOperation::Subst => nm_count += 1,
                        bio::alignment::AlignmentOperation::Del => nm_count += 1,
                        bio::alignment::AlignmentOperation::Ins => nm_count += 1,
                        _ => (),
                    }
                }
                rec.push_aux(b"nM", Aux::I32(nm_count as i32)).unwrap();
                rec.push_aux(b"xf", Aux::I32(0)).unwrap();

                let mut flags = 0;
                // Set the read paired flag
                flags |= 0x1; // set the 1st bit (0x1) to indicate paired-end read

                // Set the read mapped flag
                if alignment.score > 0 {
                    flags |= 0x2; // set the 2nd bit (0x2) to indicate read is mapped
                }

                // Set the read reverse strand flag
                if alignment.ystart > alignment.yend {
                    flags |= 0x10; // set the 5th bit (0x10) to indicate reverse strand
                }

                // Set the read first/last mate flags
                if alignment.xstart < alignment.ystart {
                    flags |= 0x40; // set the 7th bit (0x40) to indicate first mate
                } else {
                    flags |= 0x80; // set the 8th bit (0x80) to indicate last mate
                }

                // Set the flags in the BAM record
                rec.set_flags(flags as u16);

                bam_writer.write(&rec).unwrap();
            }
        }
    }
}
