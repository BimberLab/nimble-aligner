use crate::align::{AlignDebugInfo, AlignFilterConfig, PseudoAligner, FilterReason, AlignmentDirection};
use crate::parse::bam::BAM_FIELDS_TO_REPORT;
use crate::parse::bam::UMIReader;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::revcomp;
use debruijn::dna_string::DnaString;
use flate2::write::GzEncoder;
use flate2::Compression;

use std::collections::HashMap;
use std::sync::mpsc::{self};
use std::sync::{Arc, Mutex};
use std::thread;

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Error, Write};

const MAX_UMIS_IN_CHANNEL: usize = 50;

fn bam_data_values(bam_data: &Vec<String>) -> String {
    bam_data
        .iter()
        .enumerate() 
        .filter(|&(index, _)| index != 15) 
        .map(|(_, value)| value.as_str())
        .collect::<Vec<&str>>()
        .join("\t")
}

fn bam_data_header(prefix: &str) -> String {
    BAM_FIELDS_TO_REPORT
        .iter()
        .enumerate() 
        .filter(|&(index, _)| index != 15) 
        .map(|(_, &field)| format!("{}_{}", prefix, field))
        .collect::<Vec<String>>()
        .join("\t")
}

pub fn process(
    input_files: Vec<String>,
    reference_indices: Vec<(PseudoAligner, PseudoAligner)>,
    reference_metadata: Vec<ReferenceMetadata>,
    align_configs: Vec<AlignFilterConfig>,
    output_paths: Vec<String>,
    num_cores: usize,
) {
    // associated with r2, r1, r2 forward, r1 forward, r2 reverse, r1 reverse, both
    let (log_sender, log_receiver) =
        mpsc::channel::<((Vec<String>, (i32, Vec<String>, Vec<String>, (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentDirection)), usize)>();

    let log_thread = thread::spawn(move || {
        println!("Spawning logging thread.");
        let mut log_files: Vec<GzEncoder<BufWriter<File>>> = output_paths
            .iter()
            .map(|path| {
                let file = OpenOptions::new()
                    .write(true)
                    .create(true)
                    .open(path)
                    .unwrap();
                let buffered_writer = BufWriter::new(file);
                GzEncoder::new(buffered_writer, Compression::default())
            })
            .collect();
        let mut first_write: Vec<bool> = vec![true; log_files.len()];

        loop {
            match log_receiver.recv() {
                Ok((msg, index)) => {
                    let file_handle = &mut log_files[index];
    
                    if first_write[index] {
                        writeln!(
                            file_handle,
                            "nimble_features\tnimble_score\t{}\t{}\t{}",
                            bam_data_header("r1"),
                            bam_data_header("r2"),
                            "r1_filter_forward\tr1_forward_score\tr1_filter_reverse\tr1_reverse_score\tr2_filter_forward\tr2_forward_score\tr2_filter_reverse\tr2_reverse_score\ttriage_reason\taligndirection",
                        )
                        .unwrap();
                        first_write[index] = false;
                    }
    
                    writeln!(
                        file_handle,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        msg.0.join(","),
                        msg.1 .0,
                        bam_data_values(&msg.1 .2), // r1
                        bam_data_values(&msg.1 .1), // r2. these are reversed from the order they take in the data structures, hence the weird numbering
                        &msg.1.4.0,
                        &msg.1.4.1,
                        &msg.1.6.0,
                        &msg.1.6.1,
                        &msg.1.3.0,
                        &msg.1.3.1,
                        &msg.1.5.0,
                        &msg.1.5.1,
                        &msg.1.7,
                        &msg.1.8,
                    )
                    .unwrap();
                }
                Err(_) => {
                    println!("Log thread received termination signal");
                    break;
                }
            }
        }

        // Explicitly finish GzEncoders to ensure all data is flushed
        for encoder in log_files.into_iter() {
            let _ = encoder.finish();
        }
    });

    let (sender, receiver) = mpsc::sync_channel(MAX_UMIS_IN_CHANNEL);
    let receiver = Arc::new(Mutex::new(receiver));

    let reference_indices = Arc::new(reference_indices);
    let reference_metadata = Arc::new(reference_metadata);
    let align_configs = Arc::new(align_configs);

    let producer_handle = thread::spawn(move || {
        println!("Spawning reader thread.");
        let mut reader = UMIReader::new(&input_files[0], false);
        let mut has_aligned = false;

        loop {
            let final_umi = reader.next();

            if final_umi && has_aligned {
                println!("Finished reading UMIs from input file.");
                break;
            } else {
                sender
                    .send((
                        reader.current_umi_group.clone(),
                        reader.current_metadata_group.clone(),
                    ))
                    .unwrap();
            }

            has_aligned = true;
        }
    });

    let num_consumers = if num_cores > 1 {
        num_cores - 1
    } else {
        num_cores
    };

    let mut consumer_handles = Vec::with_capacity(num_consumers);
    for thread_num in 0..num_consumers {
        println!("Spawning consumer thread {:?}", thread_num);
        let receiver_clone = Arc::clone(&receiver);

        let reference_indices = Arc::clone(&reference_indices);
        let reference_metadata = Arc::clone(&reference_metadata);
        let align_configs = Arc::clone(&align_configs);
        let log_sender_clone = log_sender.clone();

        let handle = thread::spawn(move || loop {
            let data = receiver_clone.lock().unwrap().recv();

            match data {
                Ok((umi, current_metadata_group)) => {
                    let results = align_umi_to_libraries(
                        umi,
                        current_metadata_group,
                        &*reference_indices,
                        &*reference_metadata,
                        &*align_configs,
                        thread_num,
                    );

                    for (i, library_scores) in results.into_iter().enumerate() {
                        for score in library_scores {
                            log_sender_clone.send((score.clone(), i)).unwrap();
                        }
                    }
                }
                Err(_) => break,
            }
        });
        consumer_handles.push(handle);
    }

    producer_handle.join().unwrap();
    println!("Joined on producer.");

    for handle in consumer_handles {
        println!("Joined on consumer.");
        handle.join().unwrap();
    }

    drop(log_sender);

    log_thread.join().unwrap();
    println!("Joined on logging; terminating.");
}

fn get_score<'a>(
    current_umi_group: &'a Vec<DnaString>,
    current_metadata_group: &'a Vec<Vec<String>>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
    reverse_comp_read: &'a Vec<bool>,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentDirection)>
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
        current_metadata_group,
        reference_index,
        &reference_metadata,
        align_config,
        debug_info,
    )
}

fn align_umi_to_libraries(
    umi: Vec<DnaString>,
    current_metadata_group: Vec<Vec<String>>,
    reference_indices: &Vec<(PseudoAligner, PseudoAligner)>,
    reference_metadata: &Vec<ReferenceMetadata>,
    align_configs: &Vec<AlignFilterConfig>,
    _thread_num: usize,
) -> Vec<Vec<(Vec<String>, (i32, Vec<String>, Vec<String>, (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentDirection))>> {
    let mut results = vec![];

    for (i, reference_index) in reference_indices.iter().enumerate() {
        let (mut s, _, filter_reasons) = get_score(
            &umi,
            &current_metadata_group,
            reference_index,
            &reference_metadata[i],
            &align_configs[i],
            None,
            &current_metadata_group
                .clone()
                .into_iter()
                .map(|v| parse_str_as_bool(&v[3]))
                .collect::<Vec<bool>>(),
        );

        if s.len() == 0 {
            results.push(vec![]);
        } else {
            let mut non_matching_reads: Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))> = Vec::new();
            let mut scored_qnames = Vec::new();

            for score in &s {
                let qname = score.1.1[0].clone();
                scored_qnames.push(qname);
            }

            for i in (0..current_metadata_group.len()).step_by(2) {
                if i + 1 < current_metadata_group.len() {
                    let pair = (current_metadata_group[i].clone(), current_metadata_group[i + 1].clone());
                    let qname = &pair.1[0];

                    if scored_qnames.contains(qname) {
                        continue;
                    }

                    non_matching_reads.push((vec![], (0, pair.0, pair.1)));
                }
            }

            s.extend(non_matching_reads);

            let mut transformed_scores = Vec::new();
            for score in s.into_iter() {
                let r1_key_part = check_reverse_comp((&DnaString::from_dna_string(&score.1.1[15]), &parse_str_as_bool(&score.1.1[3]))).to_string();
                let r2_key_part = check_reverse_comp((&DnaString::from_dna_string(&score.1.2[15]), &parse_str_as_bool(&score.1.2[3]))).to_string();
                let filter_result = filter_reasons.get(&(r1_key_part.clone() + &r2_key_part));

                let new_score = match filter_result {
                    Some(v) => {
                        (
                            score.0,
                            (
                                score.1.0,
                                score.1.1,
                                score.1.2,
                                v.0,
                                v.1,
                                v.2,
                                v.3,
                                v.4,
                                v.5,
                            ),
                        )
                    },
                    None => {
                        (
                            score.0,
                            (
                                score.1.0,
                                score.1.1,
                                score.1.2,
                                (FilterReason::None, 0),
                                (FilterReason::None, 0),
                                (FilterReason::None, 0),
                                (FilterReason::None, 0),
                                FilterReason::None,
                                AlignmentDirection::None,
                            ),
                        )
                    },
                };
                transformed_scores.push(new_score);
            }

            let s = transformed_scores;

            results.push(s);
        }
    }

    results
}

fn check_reverse_comp(rec: (&DnaString, &bool)) -> DnaString {
    let (seq, reverse_comp) = rec;

    if *reverse_comp {
        DnaString::from_dna_string(&revcomp(&seq.to_string())).to_owned()
    } else {
        seq.to_owned()
    }
}

fn parse_str_as_bool(v: &str) -> bool {
    match v {
        "true" => true,
        "false" => false,
        _ => panic!("Could not parse revcomp field \"{}\" as boolean", v)
    }
}