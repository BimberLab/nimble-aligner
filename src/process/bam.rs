use crate::align::{ AlignFilterConfig, PseudoAligner, FilterReason, AlignmentOrientation };
use crate::parse::bam::BAM_FIELDS_TO_REPORT;
use crate::parse::bam::UMIReader;
use crate::reference_library::Reference;
use crate::score::call;
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

// Multithreaded UMI-scoring pipeline for processing .bam files
pub fn process(
    input_files: Vec<String>,
    reference_indices: Vec<(PseudoAligner, PseudoAligner)>,
    references: Vec<Reference>,
    aligner_configs: Vec<AlignFilterConfig>,
    output_paths: Vec<String>,
    num_cores: usize,
) {
    // The logger mpsc channels takes the result of alignments and write them out to the count file
    let (log_sender, log_receiver) =
        mpsc::channel::<((Vec<String>, (i32, Vec<String>, Vec<String>, (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)), usize)>();

    // Spawn the logging thread
    let log_thread = thread::spawn(move || {
        println!("Spawning logging thread.");

        // For each provied output file path, create a new file with a buffered writer for writing gzipped data
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

        // List of first_writes, for determining whether to write the header
        let mut first_write: Vec<bool> = vec![true; log_files.len()];

        // Log receiver loop
        loop {
            // When a message is recieved, write the header if it's the first write. Then, write the feature scores, r1/r2 data, and all included bam headers
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
                        msg.0.join(","), // feature calls
                        msg.1 .0, // scores
                        bam_data_values(&msg.1 .2), // r1 bam data
                        bam_data_values(&msg.1 .1), // r2 bam data
                        &msg.1.4.0, // r1_filter_forward
                        &msg.1.4.1, // r1_forward_score
                        &msg.1.6.0, // r1_filter_reverse
                        &msg.1.6.1, // r1_reverse_score
                        &msg.1.3.0, // r2_filter_forward
                        &msg.1.3.1, // r2_forward_score
                        &msg.1.5.0, // r2_filter_reverse
                        &msg.1.5.1, // r2_reverse_score
                        &msg.1.7,   // triagereason
                        &msg.1.8,   // aligndirection
                    )
                    .unwrap();
                }
                Err(_) => {
                    println!("Log thread received termination signal");
                    break;
                }
            }
        }

        // Flush buffers once termination signal is received
        for encoder in log_files.into_iter() {
            let _ = encoder.finish();
        }
    });

    // Sender and reciever manage sending UMIs to the scoring pipeline, then sending those scores to the logger
    let (sender, receiver) = mpsc::sync_channel(MAX_UMIS_IN_CHANNEL);
    let receiver = Arc::new(Mutex::new(receiver));

    let reference_indices = Arc::new(reference_indices);
    let references = Arc::new(references);
    let aligner_configs = Arc::new(aligner_configs);

    // Spawn the producer thread, which reads sequences from the .bam file grouped on UMI 
    let producer_handle = thread::spawn(move || {
        println!("Spawning reader thread.");
        let mut reader = UMIReader::new(&input_files[0], false);
        let mut has_aligned = false;

        // Iterate the reader, sending the UMI sequence data and all associated .bam metadata to the aligner
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

    // We reserve one thread for reading and logging -- the rest are used for alignment in the consumer threads
    let num_consumers = if num_cores > 1 {
        num_cores - 1
    } else {
        num_cores
    };

    // Spawn num_consumers alignment threads
    let mut consumer_handles = Vec::with_capacity(num_consumers);
    for thread_num in 0..num_consumers {
        println!("Spawning consumer thread {:?}", thread_num);
        let receiver_clone = Arc::clone(&receiver);

        let reference_indices = Arc::clone(&reference_indices);
        let reference_metadata = Arc::clone(&references);
        let align_configs = Arc::clone(&aligner_configs);
        let log_sender_clone = log_sender.clone();

        let handle = thread::spawn(move || loop {
            // Loop and read data from the .bam reader
            let data = receiver_clone.lock().unwrap().recv();

            // If there's data, send it to the alignment pipeline
            match data {
                Ok((umi, current_metadata_group)) => {
                    let results = align_umi_to_libraries(
                        umi,
                        current_metadata_group,
                        &*reference_indices,
                        &*reference_metadata,
                        &*align_configs
                    );

                    // Send the resulting data to the logging thread, one send per library result
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

    // Once there's no more data to produce, join
    producer_handle.join().unwrap();
    println!("Joined on producer.");

    // Once the producer has joined, join the consumers once they're finished with the final alignments
    for handle in consumer_handles {
        println!("Joined on consumer.");
        handle.join().unwrap();
    }

    // Finally, drop the log sender and join the log_thraed once it's done aligning
    drop(log_sender);

    log_thread.join().unwrap();
    println!("Joined on logging; terminating.");
}

fn get_calls<'a>(
    umi: &'a Vec<DnaString>,
    umi_metadata: &'a Vec<Vec<String>>,
    reference_index: &(PseudoAligner, PseudoAligner),
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
    reverse_comp_read: &'a Vec<bool>,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)>
) {
    // Get iterators to sequences and their corresponding read-mates from the current UMI 
    // We send 2 iterators per set of r1/r2 sequences
    // The sequences are reverse complemented in accorance with the provided reverse_comp_read list, which contains true/false per input sequence
    let sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
       umi 
            .iter()
            .zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(reverse_comp_if_needed(rec))),
    );

    let sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        umi
            .iter()
            .zip(reverse_comp_read.iter())
            .step_by(2)
            .map(|rec| Ok(reverse_comp_if_needed(rec))),
    );

    let mate_sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        umi
            .iter()
            .zip(reverse_comp_read.iter())
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(reverse_comp_if_needed(rec))),
    );

    let mate_sequences_clone: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a> = Box::new(
        umi
            .iter()
            .zip(reverse_comp_read.iter())
            .skip(1)
            .step_by(2)
            .map(|rec| Ok(reverse_comp_if_needed(rec))),
    );

    // Using the set of sequences and their mates, alongside the UMI's metadata, the library index and library data, and the aligner configuration, produce a set of calls
    call(
        (sequences, sequences_clone),
        Some((mate_sequences, mate_sequences_clone)),
        umi_metadata,
        reference_index,
        &reference,
        aligner_config,
    )
}

fn align_umi_to_libraries(
    umi: Vec<DnaString>,
    umi_metadata: Vec<Vec<String>>,
    reference_indices: &Vec<(PseudoAligner, PseudoAligner)>,
    references: &Vec<Reference>,
    aligner_configs: &Vec<AlignFilterConfig>,
) -> Vec<Vec<(Vec<String>, (i32, Vec<String>, Vec<String>, (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation))>> {
    let mut results = vec![];

    // For each library, send the data into the alignment pipeline for scoring
    for (i, reference_index) in reference_indices.iter().enumerate() {
        let (mut s, _, filter_reasons) = get_calls(
            &umi,
            &umi_metadata,
            reference_index,
            &references[i],
            &aligner_configs[i],
            &umi_metadata
                .clone()
                .into_iter()
                .map(|v| parse_str_as_bool(&v[3])) // set the sequence to be reverse complemented if it has been reverse complemented in the bam record
                .collect::<Vec<bool>>(),
        );

        

        // TODO refactor below, reliant on get_score being refactored first, probably. I think it's responsible for adding additional rows per-UMI
        if s.len() == 0 {
            results.push(vec![]);
        } else {
            let mut non_matching_reads: Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))> = Vec::new();
            let mut scored_qnames = Vec::new();

            for score in &s {
                let qname = score.1.1[0].clone();
                scored_qnames.push(qname);
            }

            for i in (0..umi_metadata.len()).step_by(2) {
                if i + 1 < umi_metadata.len() {
                    let pair = (umi_metadata[i].clone(), umi_metadata[i + 1].clone());
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
                let r1_key_part = reverse_comp_if_needed((&DnaString::from_dna_string(&score.1.1[15]), &parse_str_as_bool(&score.1.1[3]))).to_string();
                let r2_key_part = reverse_comp_if_needed((&DnaString::from_dna_string(&score.1.2[15]), &parse_str_as_bool(&score.1.2[3]))).to_string();
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
                                AlignmentOrientation::None,
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

fn reverse_comp_if_needed(rec: (&DnaString, &bool)) -> DnaString {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_comp_if_needed() {
        let dna = DnaString::from_dna_string("ATGC");
        let reverse_needed = true;
        let no_reverse_needed = false;

        let result_with_reverse = reverse_comp_if_needed((&dna, &reverse_needed));
        assert_eq!(result_with_reverse.to_string(), "GCAT");

        let result_without_reverse = reverse_comp_if_needed((&dna, &no_reverse_needed));
        assert_eq!(result_without_reverse.to_string(), "ATGC");
    }

    #[test]
    fn test_parse_str_as_bool_true() {
        let input = "true";
        assert!(parse_str_as_bool(input));
    }

    #[test]
    fn test_parse_str_as_bool_false() {
        let input = "false";
        assert!(!parse_str_as_bool(input));
    }

    #[test]
    #[should_panic(expected = "Could not parse revcomp field \"invalid\" as boolean")]
    fn test_parse_str_as_bool_panic() {
        let input = "invalid";
        parse_str_as_bool(input);
    }
}