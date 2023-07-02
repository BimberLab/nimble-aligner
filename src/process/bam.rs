use crate::align::BamData;
use crate::align::{AlignDebugInfo, AlignFilterConfig, PseudoAligner};
use crate::parse::bam::UMIReader;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::revcomp;
use debruijn::dna_string::DnaString;

use std::sync::mpsc::{self};
use std::sync::{Arc, Mutex};
use std::thread;

use std::fs::{File, OpenOptions};
use std::io::{Error, Write};

const MAX_UMIS_IN_CHANNEL: usize = 100;

pub fn process(
    input_files: Vec<String>,
    reference_indices: Vec<(PseudoAligner, PseudoAligner)>,
    reference_metadata: Vec<ReferenceMetadata>,
    align_configs: Vec<AlignFilterConfig>,
    output_paths: Vec<String>,
    num_cores: usize,
) {
    let (log_sender, log_receiver) =
        mpsc::channel::<((Vec<String>, (i32, BamData, BamData)), usize)>();

    let log_thread = thread::spawn(move || {
        println!("Spawning logging thread.");
        let mut log_files: Vec<File> = output_paths
            .iter()
            .map(|path| {
                OpenOptions::new()
                    .write(true)
                    .create(true)
                    .open(path)
                    .unwrap()
            })
            .collect();
        let mut first_write: Vec<bool> = vec![true; log_files.len()];

        loop {
            match log_receiver.recv() {
                Ok((msg, index)) => {
                    let file_handle: &mut File = &mut log_files[index];

                    if first_write[index] {
                        write!(file_handle, "features\tscore\tumi\tcb\n").unwrap();
                        first_write[index] = false;
                    }

                    write!(
                        file_handle,
                        "{}\t{}\t{}\t{}\n",
                        msg.0.join(","),
                        msg.1 .0,
                        msg.1 .1.umi,
                        msg.1 .1.cb
                    )
                    .unwrap();
                }
                Err(_) => {
                    println!("Log thread received termination signal");
                    break;
                }
            }
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
    current_metadata_group: &'a Vec<(
        u8,
        String,
        String,
        bool,
        String,
        Vec<u8>,
        Vec<u8>,
        String,
        String,
        String,
        String,
    )>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
    reverse_comp_read: &'a Vec<bool>,
) -> (
    Vec<(Vec<String>, (i32, BamData, BamData))>,
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
        current_metadata_group,
        reference_index,
        &reference_metadata,
        align_config,
        debug_info,
    )
}

fn align_umi_to_libraries(
    umi: Vec<DnaString>,
    current_metadata_group: Vec<(
        u8,
        String,
        String,
        bool,
        String,
        Vec<u8>,
        Vec<u8>,
        String,
        String,
        String,
        String,
    )>,
    reference_indices: &Vec<(PseudoAligner, PseudoAligner)>,
    reference_metadata: &Vec<ReferenceMetadata>,
    align_configs: &Vec<AlignFilterConfig>,
    _thread_num: usize,
) -> Vec<Vec<(Vec<String>, (i32, BamData, BamData))>> {
    let mut results = vec![];

    for (i, reference_index) in reference_indices.iter().enumerate() {
        let (s, _) = get_score(
            &umi,
            &current_metadata_group,
            reference_index,
            &reference_metadata[i],
            &align_configs[i],
            None,
            &current_metadata_group
                .clone()
                .into_iter()
                .map(|(_, _, _, r, _, _, _, _, _, _, _)| r)
                .collect::<Vec<bool>>(),
        );

        if s.len() == 0 {
            results.push(vec![]);
        } else {
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
