use crate::align::BamData;
//use crate::ALLOCATOR;
use crate::align::{AlignDebugInfo, AlignFilterConfig, PseudoAligner};
use crate::parse::bam::UMIReader;
use crate::reference_library::ReferenceMetadata;
use crate::score::score;
use crate::utils::revcomp;
use debruijn::dna_string::DnaString;
use std::time::Duration;

use std::sync::mpsc::{self};
use std::sync::{Arc, Mutex};
use std::thread;

use std::fs::{File, OpenOptions};
use std::io::{Error, Write};
use std::time::Instant;

const MAX_UMIS_IN_CHANNEL: usize = 50;
const SAFETY_BUFFER: f64 = 1.20;
const WAIT_TIMEOUT: Duration = Duration::from_secs(5);
const SLEEP_DURATION: Duration = Duration::from_millis(50);

fn bam_data_values(bam_data: &BamData) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        bam_data.sequence,
        bam_data.mapq,
        bam_data.orientation,
        bam_data.pair,
        bam_data.rev_comp,
        bam_data.hit,
        bam_data.qname,
        //bam_data.qual.iter().map(|&q| q as char).collect::<String>(),
        bam_data.tx,
        bam_data.umi,
        bam_data.cb,
        bam_data.an
    )
}

fn bam_data_header(prefix: &str) -> String {
    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        format!("{}_sequence", prefix),
        format!("{}_mapq", prefix),
        format!("{}_orientation", prefix),
        format!("{}_pair", prefix),
        format!("{}_rev_comp", prefix),
        format!("{}_hit", prefix),
        format!("{}_qname", prefix),
        format!("{}_tx", prefix),
        format!("{}_umi", prefix),
        format!("{}_cb", prefix),
        format!("{}_an", prefix))
}

struct OwnedMetadataWrapper {
    pub sequence: Box<str>,
    pub mapq: Box<u8>,
    pub orientation: Box<str>,
    pub pair: Box<str>,
    pub rev_comp: Box<bool>,
    pub hit: Box<str>,
    pub qname: Box<[u8]>,
    pub qual: Box<[u8]>,
    pub tx: Box<str>,
    pub umi: Box<str>,
    pub cb: Box<str>,
    pub an: Box<str>,
}

fn convert_message_to_alignable(metadata: Vec<Box<OwnedMetadataWrapper>>) -> Vec<(
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
    )> {
    metadata.into_iter().map(|wrapper| {
        (
            *wrapper.mapq,
            (*wrapper.orientation).to_string(),
            (*wrapper.pair).to_string(),
            *wrapper.rev_comp,
            (*wrapper.hit).to_string(),
            (*wrapper.qname).to_vec(),
            (*wrapper.qual).to_vec(),
            (*wrapper.tx).to_string(),
            (*wrapper.umi).to_string(),
            (*wrapper.cb).to_string(),
            (*wrapper.an).to_string(),
        )
    }).collect()
}

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
                        write!(
                            file_handle,
                            "features\tscore\t{}\t{}\n",
                            bam_data_header("r1"),
                            bam_data_header("r2")
                        )
                        .unwrap();
                        first_write[index] = false;
                    }
            
                    write!(
                        file_handle,
                        "{}\t{}\t{}\t{}\n",
                        msg.0.join(","),
                        msg.1 .0,
                        bam_data_values(&msg.1 .2), // r1
                        bam_data_values(&msg.1 .1) // r2
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
                // Wrapper that stops endless memory leak from strings and vecs getting lost in the
                // sender queue. Seems like there's still some odd behavior around qual/tx/an, so
                // they're getting dropped from the output for the time being.
                let metadata: Vec<Box<OwnedMetadataWrapper>> = reader.current_metadata_group.iter().map(|tuple| {
                    Box::new(OwnedMetadataWrapper {
                        sequence: String::new().into_boxed_str(),
                        mapq: Box::new(tuple.0),
                        rev_comp: Box::new(tuple.3),
                        qual: Vec::new().into_boxed_slice(),
                        tx: String::new().into_boxed_str(),
                        an: String::new().into_boxed_str(),
                        orientation: tuple.1.clone().into_boxed_str(),
                        pair: tuple.2.clone().into_boxed_str(),
                        hit: tuple.4.clone().into_boxed_str(),
                        qname: tuple.5.clone().into_boxed_slice(),
                        //qual: tuple.6.clone().into_boxed_slice(),
                        //tx: tuple.7.clone().into_boxed_str(),
                        umi: tuple.8.clone().into_boxed_str(),
                        cb: tuple.9.clone().into_boxed_str(),
                        //an: tuple.10.clone().into_boxed_str(),
                    })
                }).collect();

                // Sender uses a wrapped object w/ Boxes instead of the clone of the data
                sender
                    .send((
                        reader.current_umi_group.clone(),
                        //reader.current_metadata_group.clone(),
                        metadata,
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

            //let _safe_to_allocate = block_on_memory_headroom(num_consumers);

            match data {
                Ok((umi, owned_metadata_group)) => {

                    // Use new message structure that reduces memory leakage from the BAM flags
                    let current_metadata_group = convert_message_to_alignable(owned_metadata_group);
                    //let current_metadata_group = owned_metadata_group;

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

fn block_on_memory_headroom(num_consumers: usize) -> bool {
    /*let start_time = Instant::now();
    let total_memory = ALLOCATOR.limit();

    let num_cores_for_avg = if num_consumers > 1 {
        num_consumers - 1
    } else {
        num_consumers
    };

    while start_time.elapsed() < WAIT_TIMEOUT {
        let current_memory = ALLOCATOR.allocated();
        let average_memory_per_thread = current_memory / num_cores_for_avg;
        let remaining_headroom = total_memory - current_memory;
        let predicted_memory_use = (average_memory_per_thread as f64) * SAFETY_BUFFER;

        if (remaining_headroom as f64) > predicted_memory_use {
            return true; // Hopefully safe to proceed
        }

        // Sleep for a short duration to prevent busy waiting
        thread::sleep(SLEEP_DURATION);
    }*/

    false // Timeout reached, moving forward may cause an OOM abort
}
