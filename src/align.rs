use crate::filter;
use crate::reference_library;

use std::collections::HashMap;
use std::fmt;
use std::io::Error;

use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::Path;

use std::fmt::{Display, Formatter};
use std::time::{Duration, Instant};

use array_tool::vec::Intersect;
use array_tool::vec::Uniq;
use debruijn::dna_string::DnaString;
use lexical_sort::{natural_lexical_cmp, StringSort};
use reference_library::ReferenceMetadata;

const MIN_READ_LENGTH: usize = 12;

pub type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::Kmer30>;

#[derive(Debug)]
pub enum IntersectLevel {
    NoIntersect,
    IntersectWithFallback,
    ForceIntersect,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum FilterReason {
    ScoreBelowThreshold,
    DiscardedMultipleMatch,
    DiscardedNonzeroMismatch,
    NoMatch,
    NoMatchAndScoreBelowThreshold,
    DifferentFilterReasons,
    NotMatchingPair,
    ForceIntersectFailure,
    ShortRead,
    MaxHitsExceeded,
}

impl Display for FilterReason {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match *self {
            FilterReason::ScoreBelowThreshold => write!(f, "Score Below Threshold"),
            FilterReason::DiscardedMultipleMatch => write!(f, "Discarded Multiple Match"),
            FilterReason::DiscardedNonzeroMismatch => write!(f, "Discarded Nonzero Mismatch"),
            FilterReason::NoMatch => write!(f, "No Match"),
            FilterReason::NoMatchAndScoreBelowThreshold => {
                write!(f, "No Match and Score Below Threshold")
            }
            FilterReason::DifferentFilterReasons => write!(f, "Different Filter Reasons"),
            FilterReason::NotMatchingPair => write!(f, "Not Matching Pair"),
            FilterReason::ForceIntersectFailure => write!(f, "Force Intersect Failure"),
            FilterReason::ShortRead => write!(f, "Short Read"),
            FilterReason::MaxHitsExceeded => write!(f, "Max Hits Exceeded"),
        }
    }
}

#[derive(Debug)]
pub struct AlignFilterConfig {
    pub reference_genome_size: usize,
    pub score_percent: f64,
    pub score_threshold: usize,
    pub num_mismatches: usize,
    pub discard_nonzero_mismatch: bool,
    pub discard_multiple_matches: bool,
    pub score_filter: i32,
    pub intersect_level: IntersectLevel,
    pub require_valid_pair: bool,
    pub discard_multi_hits: usize,
    pub max_hits_to_report: usize,
    pub strand_filter: StrandFilter,
}

#[derive(Debug, Copy, Clone)]
pub enum StrandFilter {
    Unstranded,
    FivePrime,
    ThreePrime,
    None,
}

#[derive(Default)]
pub struct AlignDebugInfo {
    pub debug_file: String,
    pub read_units_aligned: usize,
    pub score_below_threshold: usize,
    pub discarded_multiple_match: usize,
    pub discarded_nonzero_mismatch: usize,
    pub no_match: usize,
    pub no_match_and_score_below_threshold: usize,
    pub different_filter_reasons: usize,
    pub not_matching_pair: usize,
    pub force_intersect_failure: usize,
    pub short_read: usize,
    pub max_hits_exceeded: usize,
    pub ff_reported: usize,
    pub rr_reported: usize,
    pub uu_reported: usize,
    pub fr_reported: usize,
    pub fu_reported: usize,
    pub rf_reported: usize,
    pub ru_reported: usize,
    pub uf_reported: usize,
    pub ur_reported: usize,
    pub number_cr_skipped: usize,
}

#[derive(Debug, Copy, Clone)]
pub enum AlignmentDirection {
    FF,
    RR,
    UU,
    FR,
    FU,
    RF,
    RU,
    UF,
    UR,
}
impl fmt::Display for AlignmentDirection {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentDirection::FF => write!(f, "FF"),
            AlignmentDirection::RR => write!(f, "RR"),
            AlignmentDirection::UU => write!(f, "UU"),
            AlignmentDirection::FR => write!(f, "FR"),
            AlignmentDirection::FU => write!(f, "FU"),
            AlignmentDirection::RF => write!(f, "RF"),
            AlignmentDirection::RU => write!(f, "RU"),
            AlignmentDirection::UF => write!(f, "UF"),
            AlignmentDirection::UR => write!(f, "UR"),
        }
    }
}

fn replace_key_with_strand_dir(
    key: String,
    matched_seqs: &mut Vec<(Vec<String>, String, f64, usize, String)>,
    dir: AlignmentDirection,
) {
    for elem in matched_seqs {
        if elem.4 == key {
            elem.4 = dir.to_string()
        } else {
            elem.4 = String::new()
        }
    }
}

impl AlignmentDirection {
    fn get_alignment_dir(
        forward_pair_state: PairState,
        reverse_pair_state: PairState,
        debug_info: &mut AlignDebugInfo,
    ) -> AlignmentDirection {
        let dir = match (forward_pair_state, reverse_pair_state) {
            (PairState::First, PairState::First) => AlignmentDirection::FF,
            (PairState::First, PairState::Second) => AlignmentDirection::FR,
            (PairState::Second, PairState::First) => AlignmentDirection::RF,
            (PairState::Second, PairState::Second) => AlignmentDirection::RR,
            (PairState::First, PairState::None) => AlignmentDirection::FU,
            (PairState::Second, PairState::None) => AlignmentDirection::RU,
            (PairState::None, PairState::First) => AlignmentDirection::UF,
            (PairState::None, PairState::Second) => AlignmentDirection::UR,

            (PairState::Both, PairState::First) => AlignmentDirection::UU,
            (PairState::First, PairState::Both) => AlignmentDirection::UU,
            (PairState::Both, PairState::Second) => AlignmentDirection::UU,
            (PairState::Second, PairState::Both) => AlignmentDirection::UU,
            (PairState::Both, PairState::Both) => AlignmentDirection::UU,
            (PairState::Both, PairState::None) => AlignmentDirection::UU,
            (PairState::None, PairState::Both) => AlignmentDirection::UU,
            (PairState::None, PairState::None) => AlignmentDirection::UU,
            /*(PairState::Both, PairState::First) => AlignmentDirection::FF,
            (PairState::First, PairState::Both) => AlignmentDirection::FF,
            (PairState::Both, PairState::Second) => AlignmentDirection::RR,
            (PairState::Second, PairState::Both) => AlignmentDirection::RR,
            (PairState::Both, PairState::Both) => AlignmentDirection::FR,
            (PairState::Both, PairState::None) => AlignmentDirection::FU,
            (PairState::None, PairState::Both) => AlignmentDirection::UF,
            (PairState::None, PairState::None) => AlignmentDirection::UU,*/
        };

        match dir {
            AlignmentDirection::FF => debug_info.ff_reported += 1,
            AlignmentDirection::RR => debug_info.rr_reported += 1,
            AlignmentDirection::UU => debug_info.uu_reported += 1,
            AlignmentDirection::FR => debug_info.fr_reported += 1,
            AlignmentDirection::FU => debug_info.fu_reported += 1,
            AlignmentDirection::RF => debug_info.rf_reported += 1,
            AlignmentDirection::RU => debug_info.ru_reported += 1,
            AlignmentDirection::UF => debug_info.uf_reported += 1,
            AlignmentDirection::UR => debug_info.ur_reported += 1,
        };

        dir
    }

    fn filter_hits(
        mut forward_hits: Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
        mut reverse_hits: Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
        results: &mut HashMap<Vec<String>, i32>,
        reference_metadata: &ReferenceMetadata,
        config: &AlignFilterConfig,
        debug_info: &mut AlignDebugInfo,
        matched_sequences: &mut Vec<(Vec<String>, String, f64, usize, String)>,
        key: String,
        debug: bool,
    ) {
        if let Some((f_pair_state, _, _)) = forward_hits {
            if let Some((r_pair_state, _, _)) = reverse_hits {
                let dir =
                    AlignmentDirection::get_alignment_dir(f_pair_state, r_pair_state, debug_info);
                if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                    println!(
                        "filtered read-pair: {:?}, due to alignment dir {:?}",
                        &key, dir
                    );

                    if debug {
                        replace_key_with_strand_dir(key, matched_sequences, dir);
                    }
                    return;
                }
            }

            let dir =
                AlignmentDirection::get_alignment_dir(f_pair_state, PairState::None, debug_info);
            if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                println!(
                    "filtered read-pair: {:?}, due to alignment dir {:?}",
                    &key, dir
                );

                if debug {
                    replace_key_with_strand_dir(key, matched_sequences, dir);
                }
                return;
            }
        } else if let Some((r_pair_state, _, _)) = reverse_hits {
            let dir =
                AlignmentDirection::get_alignment_dir(PairState::None, r_pair_state, debug_info);
            if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                println!(
                    "filtered read-pair: {:?}, due to alignment dir {:?}",
                    &key, dir
                );

                if debug {
                    replace_key_with_strand_dir(key, matched_sequences, dir);
                }
                return;
            }
        }

        // Take the "best" alignment. The specific behavior is determined by the intersect level set in the aligner config
        let match_eqv_class = match config.intersect_level {
            IntersectLevel::NoIntersect => get_all_reads(&mut forward_hits, &mut reverse_hits),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(&mut forward_hits, &mut reverse_hits, true, debug_info)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(&mut forward_hits, &mut reverse_hits, false, debug_info)
            }
        };

        let mut key = get_score_map_key(&match_eqv_class, reference_metadata, &config); // Process the equivalence class into a score key
        key.string_sort_unstable(natural_lexical_cmp); // Sort for deterministic names

        // Max hits filter
        if key.len() > config.max_hits_to_report {
            debug_info.update(Some((FilterReason::MaxHitsExceeded, 0.0, 0)));
            return;
        }

        // Ensure we don't add empty keys, if any got to this point
        if key.len() == 0 {
            return;
        }

        let accessor = results.entry(key).or_insert(0);
        *accessor += 1;
    }

    fn filter_read(dir: AlignmentDirection, lib_type: &StrandFilter) -> bool {
        match lib_type {
            StrandFilter::Unstranded => AlignmentDirection::filter_unstranded(dir),
            StrandFilter::FivePrime => AlignmentDirection::filter_fiveprime(dir),
            StrandFilter::ThreePrime => AlignmentDirection::filter_threeprime(dir),
            StrandFilter::None => false,
        }
    }

    fn filter_unstranded(dir: AlignmentDirection) -> bool {
        match dir {
            AlignmentDirection::FF => true,
            AlignmentDirection::RR => true,
            AlignmentDirection::UU => true,
            AlignmentDirection::FR => false,
            AlignmentDirection::FU => false,
            AlignmentDirection::RF => false,
            AlignmentDirection::RU => false,
            AlignmentDirection::UF => false,
            AlignmentDirection::UR => false,
        }
    }

    fn filter_fiveprime(dir: AlignmentDirection) -> bool {
        match dir {
            AlignmentDirection::FF => true,
            AlignmentDirection::RR => true,
            AlignmentDirection::UU => true,
            AlignmentDirection::FR => false,
            AlignmentDirection::FU => false,
            AlignmentDirection::RF => false,
            AlignmentDirection::RU => false,
            AlignmentDirection::UF => false,
            AlignmentDirection::UR => false,
        }
    }

    fn filter_threeprime(dir: AlignmentDirection) -> bool {
        match dir {
            AlignmentDirection::FF => true,
            AlignmentDirection::RR => true,
            AlignmentDirection::UU => true,
            AlignmentDirection::FR => false,
            AlignmentDirection::FU => false,
            AlignmentDirection::RF => false,
            AlignmentDirection::RU => false,
            AlignmentDirection::UF => false,
            AlignmentDirection::UR => false,
        }
    }
}

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum PairState {
    First,
    Second,
    Both,
    None,
}

impl AlignDebugInfo {
    fn update(&mut self, reason: Option<(FilterReason, f64, usize)>) -> (f64, usize) {
        match reason {
            Some((FilterReason::ScoreBelowThreshold, _, _)) => self.score_below_threshold += 1,
            Some((FilterReason::DiscardedMultipleMatch, _, _)) => {
                self.discarded_multiple_match += 1
            }
            Some((FilterReason::DiscardedNonzeroMismatch, _, _)) => {
                self.discarded_nonzero_mismatch += 1
            }
            Some((FilterReason::NoMatch, _, _)) => self.no_match += 1,
            Some((FilterReason::NoMatchAndScoreBelowThreshold, _, _)) => {
                self.no_match_and_score_below_threshold += 1
            }
            Some((FilterReason::DifferentFilterReasons, _, _)) => {
                self.different_filter_reasons += 1
            }
            Some((FilterReason::NotMatchingPair, _, _)) => self.not_matching_pair += 1,
            Some((FilterReason::ForceIntersectFailure, _, _)) => self.force_intersect_failure += 1,
            Some((FilterReason::ShortRead, _, _)) => self.short_read += 1,
            Some((FilterReason::MaxHitsExceeded, _, _)) => self.max_hits_exceeded += 1,
            None => (),
        };

        if let Some((_, s, ns)) = reason {
            (s, ns)
        } else {
            (0.0, 0)
        }
    }

    fn merge(&mut self, info: AlignDebugInfo) {
        self.read_units_aligned += info.read_units_aligned;
        self.score_below_threshold += info.score_below_threshold;
        self.discarded_multiple_match += info.discarded_multiple_match;
        self.discarded_nonzero_mismatch += info.discarded_nonzero_mismatch;
        self.no_match += info.no_match;
        self.no_match_and_score_below_threshold += info.no_match_and_score_below_threshold;
        self.different_filter_reasons += info.different_filter_reasons;
        self.not_matching_pair += info.not_matching_pair;
        self.force_intersect_failure += info.force_intersect_failure;
        self.short_read += info.short_read;
        self.ff_reported += info.ff_reported;
        self.rr_reported += info.rr_reported;
        self.uu_reported += info.uu_reported;
        self.fr_reported += info.fr_reported;
        self.fu_reported += info.fu_reported;
        self.rf_reported += info.rf_reported;
        self.ru_reported += info.ru_reported;
        self.uf_reported += info.uf_reported;
        self.ur_reported += info.ur_reported;
    }

    pub fn get_total_attempted_reads(&self) -> usize {
        return self.read_units_aligned
            + self.score_below_threshold
            + self.discarded_multiple_match
            + self.discarded_nonzero_mismatch
            + self.no_match
            + self.no_match_and_score_below_threshold
            + self.different_filter_reasons
            + self.not_matching_pair
            + self.force_intersect_failure
            + self.short_read;
    }
}

/* Takes a set of sequences and optionally, reverse sequences, a debrujin map index of the reference
 * genome, the reference library metadata object, and the aligner configuration, and performs a
 * debrujin-graph based pseduoalignment, returning a score for each readable reference in the reference
 * genome.
 * This function does some alignment-time filtration based on the provided configuration. */
pub fn score<'a>(
    sequence_iter_pair: (
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    ),
    reverse_sequence_iter_pair: Option<(
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    )>,
    current_metadata_group: &Vec<(
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
    index_pair: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
    output_path: &str,
) -> (
    Vec<(Vec<String>, i32)>,
    Vec<(Vec<String>, String, f64, usize, String)>,
) {
    let (index_forward, index_backward) = index_pair;
    let (sequences, sequences_2) = sequence_iter_pair;
    let (reverse_sequences, reverse_sequences_2) = match reverse_sequence_iter_pair {
        Some((l, r)) => (Some(l), Some(r)),
        None => (None, None),
    };

    let (forward_score, mut forward_matched_sequences, mut forward_align_debug_info) =
        generate_score(
            sequences,
            reverse_sequences,
            current_metadata_group,
            index_forward,
            reference_metadata,
            config,
            String::from("forward"),
            output_path,
        );
    let (mut backward_score, mut backward_matched_sequences, backward_align_debug_info) =
        generate_score(
            sequences_2,
            reverse_sequences_2,
            current_metadata_group,
            index_backward,
            reference_metadata,
            config,
            String::from("reverse"),
            output_path,
        );

    forward_matched_sequences.append(&mut backward_matched_sequences);

    forward_align_debug_info.merge(backward_align_debug_info);
    let mut results: HashMap<Vec<String>, i32> = HashMap::new();

    for (key, f) in forward_score.into_iter() {
        let r = match backward_score.get(&key) {
            Some(_) => backward_score.remove(&key),
            None => None,
        };

        AlignmentDirection::filter_hits(
            Some(f),
            r,
            &mut results,
            reference_metadata,
            config,
            &mut forward_align_debug_info,
            &mut forward_matched_sequences,
            key,
            debug_info.is_some(),
        );
    }

    for (key, r) in backward_score.into_iter() {
        AlignmentDirection::filter_hits(
            None,
            Some(r),
            &mut results,
            reference_metadata,
            config,
            &mut forward_align_debug_info,
            &mut forward_matched_sequences,
            key,
            debug_info.is_some(),
        );
    }

    // Update the results map
    let mut ret = Vec::new();
    for (key, value) in results.into_iter() {
        ret.push((key, value));
    }

    if let Some(debug_info) = debug_info {
        debug_info.merge(forward_align_debug_info);
    }

    (ret, forward_matched_sequences)
}

fn generate_score<'a>(
    sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    mut reverse_sequences: Option<Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>>,
    current_metadata_group: &Vec<(
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
    index: &PseudoAligner,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
    reference_orientation: String,
    output_path: &str,
) -> (
    HashMap<String, (PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    AlignDebugInfo,
) {
    // HashMap of the alignment results. The keys are either strong hits or equivalence classes of hits
    let mut score_map: HashMap<
        String,
        (PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>),
    > = HashMap::new();
    let mut debug_info: AlignDebugInfo = Default::default();
    let mut read_matches: Vec<(Vec<String>, String, f64, usize, String)> = Vec::new();

    //let score_start = Instant::now();

    // Iterate over every read/reverse read pair and align it, incrementing scores for the matching references/equivalence classes
    let mut metadata_iter = current_metadata_group.iter();
    for read in sequences {
        let forward_metadata_raw = metadata_iter.next().unwrap();
        let reverse_metadata_raw = metadata_iter.next().unwrap();

        let (mapq, orientation, pair, rev_comp, hit, qname, qual, tx, umi, cb, an) =
            forward_metadata_raw;
        let (mapq2, orientation2, pair2, rev_comp2, hit2, qname2, qual2, tx2, umi2, cb2, an2) =
            reverse_metadata_raw;

        let qname_ascii = match String::from_utf8(qname.clone()) {
            Ok(string) => string,
            Err(e) => {
                eprintln!("Error: {}", e);
                String::new()
            }
        };

        let qname_ascii2 = match String::from_utf8(qname2.clone()) {
            Ok(string) => string,
            Err(e) => {
                eprintln!("Error: {}", e);
                String::new()
            }
        };

        let forward_metadata: (
            &u8,
            &String,
            &String,
            &bool,
            &String,
            String,
            &String,
            &String,
            &String,
            &String,
        ) = (
            mapq,
            orientation,
            pair,
            rev_comp,
            hit,
            qname_ascii,
            tx,
            umi,
            cb,
            an,
        );
        let reverse_metadata = (
            mapq2,
            orientation2,
            pair2,
            rev_comp2,
            hit2,
            qname_ascii2,
            tx2,
            umi2,
            cb2,
            an2,
        );

        let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
        let mut read_rev = None;

        /* Generate score and equivalence class for this read by aligning the sequence against
         * the current reference, if there is a match.*/
        //let pseudo_start = Instant::now();
        let (seq_score, forward_filter_reason) = pseudoalign(&read, index, &config);
        //let duration = pseudo_start.elapsed();
        //println!("time to pseudoalign: {:?}", duration);

        //let start = Instant::now();

        // If there's a reversed sequence, do the paired-end alignment
        let mut rev_seq_score = None;
        let mut rev_filter_reason: Option<(FilterReason, f64, usize)> = None;
        let mut reverse_read_t = DnaString::blank(0);
        if let Some(itr) = &mut reverse_sequences {
            let reverse_read = itr
                .next()
                .expect("Error -- read and reverse read files do not have matching lengths: ")
                .expect("Error -- could not parse reverse read. Input R2 data malformed.");
            let (score, reason) = pseudoalign(&reverse_read, index, &config);
            reverse_read_t = reverse_read.clone();
            rev_seq_score = Some(score);
            read_rev = Some(reverse_read);
            rev_filter_reason = reason;
        }

        let mut seq_score_names = Vec::new();
        let mut rev_seq_score_names = Vec::new();
        let mut seq_score_weighted = -1.0;
        let mut rev_seq_score_weighted = -1.0;
        let mut forward_filter_reason_unwrapped = None;
        let mut rev_filter_reason_unwrapped = None;

        if seq_score.is_some() {
            let t: &(Vec<u32>, f64, usize) = seq_score.as_ref().unwrap();
            seq_score_names = get_score_map_key(&t.0, reference_metadata, &config);
            seq_score_weighted = t.1;
        }

        if rev_seq_score.is_some() {
            if rev_seq_score.as_ref().unwrap().is_some() {
                let t: &(Vec<u32>, f64, usize) = rev_seq_score.as_ref().unwrap().as_ref().unwrap();
                rev_seq_score_names = get_score_map_key(&t.0, reference_metadata, &config);
                rev_seq_score_weighted = t.1;
            }
        }

        if forward_filter_reason.as_ref().is_some() {
            forward_filter_reason_unwrapped = Some(forward_filter_reason.as_ref().unwrap().0);
        }

        if rev_filter_reason.as_ref().is_some() {
            rev_filter_reason_unwrapped = Some(rev_filter_reason.as_ref().unwrap().0);
        }

        let feature_of_interest = "LILR";
        let mut paired_end_type = "";
        if seq_score_names.is_empty()
            && rev_seq_score_names.is_empty()
            && (forward_metadata.4.contains(feature_of_interest)
                || reverse_metadata.4.contains(feature_of_interest))
        //&& (forward_metadata.4 == feature_of_interest || reverse_metadata.4 == feature_of_interest)
        {
            //println!("{}", forward_metadata.4);
            //println!("{}", reverse_metadata.4);
            paired_end_type = "cellranger"
        } else if (!seq_score_names.is_empty() || !rev_seq_score_names.is_empty())
            && (forward_metadata.4.is_empty() && reverse_metadata.4.is_empty())
        {
            paired_end_type = "nimble"
        } else if (!seq_score_names.is_empty() || !rev_seq_score_names.is_empty()) //shared hits
            && (forward_metadata.4.contains(feature_of_interest) || reverse_metadata.4.contains(feature_of_interest))
        //&& (forward_metadata.4 == feature_of_interest || reverse_metadata.4 == feature_of_interest)
        {
            //println!("{}", forward_metadata.4);
            //println!("{}", reverse_metadata.4);
            paired_end_type = "both"
        }

        if !paired_end_type.is_empty() {
            //let file_path = Path::new("/home/hextraza/work/data/alignment_logs/alignment_dump.tsv");
            let file_path = Path::new(output_path);

            write_header_if_needed(&file_path)
                .expect("Header error -- Error writing header to TSV file");

            log_alignment_comparison(
                &file_path,
                &read.to_string(),
                &reverse_read_t.to_string(),
                &seq_score_names,
                seq_score_weighted,
                forward_filter_reason_unwrapped,
                &rev_seq_score_names,
                rev_seq_score_weighted,
                rev_filter_reason_unwrapped,
                &forward_metadata,
                &reverse_metadata,
                paired_end_type.to_owned(),
                &reference_orientation,
            )
            .expect("Log error -- Error writing to TSV file");
        }

        /*{
            println!("read: {:?}", read);
            println!("reverse read: {:?}", reverse_read_t);
            println!(
                "class, score, filter reason: {:?} | {:?} | {:?}",
                seq_score_names, seq_score_weighted, forward_filter_reason_unwrapped
            );
            println!(
                "reverse class, score, filter reason: {:?} | {:?} | {:?}",
                rev_seq_score_names, rev_seq_score_weighted, rev_filter_reason_unwrapped
            );
            println!("forward metadata: {:?}", forward_metadata);
            println!("reverse metadata: {:?}\n", reverse_metadata);
        }*/

        let (failed_score, failed_raw_score) = if reverse_sequences.is_some() {
            match (forward_filter_reason, rev_filter_reason) {
                (Some((fr, s, ns)), Some((rr, r, nr))) => {
                    if fr == rr {
                        debug_info.update(Some((fr, s, ns)))
                    } else if (fr == FilterReason::NoMatch
                        && rr == FilterReason::ScoreBelowThreshold)
                        || (rr == FilterReason::NoMatch && fr == FilterReason::ScoreBelowThreshold)
                    {
                        let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                        debug_info.update(Some((
                            FilterReason::NoMatchAndScoreBelowThreshold,
                            s,
                            ns,
                        )))
                    } else {
                        let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                        debug_info.update(Some((FilterReason::DifferentFilterReasons, s, ns)))
                    }
                }
                (None, Some((rr, r, nr))) => debug_info.update(Some((rr, r, nr))),
                (Some((fr, s, ns)), None) => debug_info.update(Some((fr, s, ns))),
                (None, None) => (0.0, 0),
            }
        } else {
            debug_info.update(forward_filter_reason)
        };

        // If there are no reverse sequences, ignore the require_valid_pair filter
        if reverse_sequences.is_some()
            && config.require_valid_pair
            && filter_pair(&seq_score, &rev_seq_score)
        {
            debug_info.update(Some((FilterReason::NotMatchingPair, 0.0, 0)));
            continue;
        }

        let (eqv_class, s, n_s) = if let Some((eqv_class, s, n_s)) = &seq_score {
            ((*eqv_class).clone(), *s, *n_s)
        } else {
            (Vec::new(), 0.0, 0)
        };

        let (r_eqv_class, r_s, r_n_s) = if let Some(Some((r_eqv_class, r_s, n_s))) = &rev_seq_score
        {
            ((*r_eqv_class).clone(), *r_s, *n_s)
        } else {
            (Vec::new(), 0.0, 0)
        };

        if !eqv_class.is_empty() || !r_eqv_class.is_empty() {
            let mut key = if !eqv_class.is_empty() {
                get_score_map_key(&eqv_class, reference_metadata, &config) // Process the equivalence class into a score key
            } else if !r_eqv_class.is_empty() {
                get_score_map_key(&r_eqv_class, reference_metadata, &config)
            } else {
                Vec::new()
            };

            key.string_sort_unstable(natural_lexical_cmp); // Sort for deterministic names

            let (v, s, r_s, n_s, r_n_s) = if !eqv_class.is_empty() && !r_eqv_class.is_empty() {
                (
                    (
                        PairState::Both,
                        Some((eqv_class, s)),
                        Some((r_eqv_class, r_s)),
                    ),
                    s,
                    r_s,
                    n_s,
                    r_n_s,
                )
            } else if !eqv_class.is_empty() {
                (
                    (PairState::First, Some((eqv_class, s)), None),
                    s,
                    0.0,
                    n_s,
                    0,
                )
            } else if !r_eqv_class.is_empty() {
                (
                    (PairState::Second, None, Some((r_eqv_class, r_s))),
                    0.0,
                    r_s,
                    0,
                    r_n_s,
                )
            } else {
                continue;
            };

            let read_key = match &read_rev {
                Some(rev) => read.to_string() + &rev.to_string(),
                None => read.to_string(),
            };

            match v.0 {
                PairState::First => {
                    read_matches.push((key.clone(), read.to_string(), s, n_s, read_key.clone()))
                }
                PairState::Second => match &read_rev {
                    Some(r) => read_matches.push((
                        key.clone(),
                        r.to_string(),
                        r_s,
                        r_n_s,
                        read_key.clone(),
                    )),
                    None => (),
                },
                PairState::Both => {
                    read_matches.push((key.clone(), read.to_string(), s, r_n_s, read_key.clone()))
                }
                PairState::None => (),
            };

            score_map.insert(read_key, v);
            debug_info.read_units_aligned += 1;
        } else {
            // If both equivalence classes are empty, the attempted alignment has failed, but we still report the failed alignment
            read_matches.push((
                Vec::new(),
                read.to_string(),
                failed_score,
                failed_raw_score,
                String::new(),
            ));
        }

        //let duration = start.elapsed();
        //println!("time to generate score for read: {:?}", duration);
    }

    //let duration = score_start.elapsed();
    //println!("time to generate score for UMI: {:?}", duration);
    (score_map, read_matches, debug_info)
}

// Determine whether a given pair of equivalence classes constitute a valid pair
fn filter_pair(
    seq_score: &Option<(Vec<u32>, f64, usize)>,
    rev_seq_score: &Option<Option<(Vec<u32>, f64, usize)>>,
) -> bool {
    // Unpack the data to check if the pairs have the same eq_class. If they both contain data, do the comparison
    if let (Some(Some((mut rev_eq_class, _, _))), Some((mut eq_class, _, _))) =
        ((rev_seq_score).clone(), (seq_score).clone())
    {
        // Sort the vectors and compare them by counting matching elements. If they don't match, don't modify the score for this read
        rev_eq_class.sort();
        eq_class.sort();
        let matching = eq_class
            .iter()
            .zip(&rev_eq_class)
            .filter(|&(classl, classr)| classl == classr)
            .count();

        if matching != eq_class.len() || matching != rev_eq_class.len() {
            return true;
        }
    } else {
        // Otherwise, require_valid_pair is on and rev_score != seq_score, or they're both None. In either case, don't modify score
        return true;
    }
    false
}

// Return matches that match in all of the given classes; if soft intersection is enabled, fall back to best read if one of the reads is empty
fn get_intersecting_reads(
    seq_score: &mut Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
    rev_seq_score: &mut Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
    fallback_on_intersect_fail: bool,
    debug_info: &mut AlignDebugInfo,
) -> Vec<u32> {
    let (seq_class, _) = match seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => {
            f.append(r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
            )
        }
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0.0),
    };

    let (r_seq_class, _) = match rev_seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => {
            f.append(r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
            )
        }
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0.0),
    };

    let class = seq_class.intersect(r_seq_class);

    if class.len() == 0 && fallback_on_intersect_fail {
        get_all_reads(seq_score, rev_seq_score)
    } else if class.len() != 0 {
        class
    } else {
        debug_info.update(Some((FilterReason::ForceIntersectFailure, 0.0, 0)));
        Vec::new()
    }
}

// Return best equivalence class out of the given classes
fn get_all_reads(
    seq_score: &mut Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
    rev_seq_score: &mut Option<(PairState, Option<(Vec<u32>, f64)>, Option<(Vec<u32>, f64)>)>,
) -> Vec<u32> {
    let (mut seq_class, _) = match seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => {
            f.append(r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
            )
        }
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0.0),
    };

    let (mut r_seq_class, _) = match rev_seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => {
            f.append(r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
            )
        }
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0.0),
    };

    seq_class.append(&mut r_seq_class);
    seq_class.unique()
}

/* Takes a equivalence class and returns a list of strings. If we're processing allele-level data, the strings will be
 * the nt_sequences of the relevant alleles. Otherwise, if we're doing a group_by, the equivalence class will be
 * filtered such that there is only one hit per group_by string (e.g. one hit per lineage) and the corresponding strings
 * (e.g. lineage name) will be returned. */
fn get_score_map_key(
    equiv_class: &Vec<u32>,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
) -> Vec<String> {
    if reference_metadata.headers[reference_metadata.group_on] == "nt_sequence" {
        equiv_class
            .into_iter()
            .map(|ref_idx| {
                reference_metadata.columns[reference_metadata.group_on][*ref_idx as usize].clone()
            })
            .collect()
    } else {
        let mut results = Vec::new();

        for ref_idx in equiv_class {
            let mut group =
                &reference_metadata.columns[reference_metadata.group_on][*ref_idx as usize];

            // If the group is an empty string, get the sequence name instead so we don't just return nothing
            if group.is_empty() {
                group = &reference_metadata.columns[reference_metadata.sequence_name_idx]
                    [*ref_idx as usize];
            }

            if !results.contains(group) {
                results.push(group.to_string());
            }
        }

        // Filter based on discard_multi_hits
        if config.discard_multi_hits > 0 && results.len() > config.discard_multi_hits {
            Vec::new()
        } else {
            results
        }
    }
}

// Align the given sequence against the given reference with a score threshold
fn pseudoalign(
    sequence: &DnaString,
    reference_index: &PseudoAligner,
    config: &AlignFilterConfig,
) -> (
    Option<(Vec<u32>, f64, usize)>,
    Option<(FilterReason, f64, usize)>,
) {
    // Filter short reads
    if sequence.to_string().len() < MIN_READ_LENGTH {
        return (None, Some((FilterReason::ShortRead, 0.0, 0)));
    };

    // Perform alignment
    match reference_index.map_read_with_mismatch(&sequence, config.num_mismatches) {
        Some((equiv_class, score, mismatches)) => {
            // Normalize score by read length
            let normalized_score = score as f64 / sequence.len() as f64;

            // Filter nonzero mismatch
            if config.discard_nonzero_mismatch && mismatches != 0 {
                return (None, Some((FilterReason::DiscardedNonzeroMismatch, 0.0, 0)));
            }

            // Filter by score and match threshold
            filter::align::filter_alignment_by_metrics(
                score,
                normalized_score,
                config.score_threshold,
                config.score_percent,
                equiv_class,
                config.discard_multiple_matches,
            )
        }
        None => (None, Some((FilterReason::NoMatch, 0.0, 0))),
    }
}

fn log_alignment_comparison(
    file_path: &Path,
    read: &str,
    reverse_read: &str,
    seq_score_names: &Vec<String>,
    seq_score_weighted: f64,
    forward_filter_reason_unwrapped: Option<FilterReason>,
    rev_seq_score_names: &Vec<String>,
    rev_seq_score_weighted: f64,
    rev_filter_reason_unwrapped: Option<FilterReason>,
    forward_metadata: &(
        &u8,
        &String,
        &String,
        &bool,
        &String,
        String,
        &String,
        &String,
        &String,
        &String,
    ),
    reverse_metadata: &(
        &u8,
        &String,
        &String,
        &bool,
        &String,
        String,
        &String,
        &String,
        &String,
        &String,
    ),
    alignment_type: String,
    reference_orientation: &String,
) -> std::io::Result<()> {
    let (mapq, orientation, pair, rev_comp, hit, qname_ascii, tx, umi, cb, an) = forward_metadata;
    let (mapq2, orientation2, pair2, rev_comp2, hit2, qname_ascii2, tx2, umi2, cb2, an2) =
        reverse_metadata;

    let mut forward_filter_reason = String::from("No Filter");
    let mut rev_filter_reason = String::from("No Filter");

    if forward_filter_reason_unwrapped.is_some() {
        forward_filter_reason = forward_filter_reason_unwrapped.unwrap().to_string();
    }

    if rev_filter_reason_unwrapped.is_some() {
        rev_filter_reason = rev_filter_reason_unwrapped.unwrap().to_string();
    }

    let mut file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_path)?;

    writeln!(
        file,
        "{}\t{}\t{:?}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{}\t{}\t{}\t{:?}\t{}\t{:?}\t{}\t{}\t{}\t{:?}\t{:?}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\n",
        read,
        reverse_read,
        seq_score_names,
        seq_score_weighted,
        forward_filter_reason,
        mapq,
        orientation,
        pair,
        rev_comp,
        hit,
        qname_ascii,
        tx,
        umi,
        cb,
        an,
        rev_seq_score_names,
        rev_seq_score_weighted,
        rev_filter_reason,
        mapq2,
        orientation2,
        pair2,
        rev_comp2,
        hit2,
        qname_ascii2,
        tx2,
        alignment_type,
        reference_orientation,
        umi2,
        cb2,
        an2,
    )?;

    Ok(())
}

fn write_header_if_needed(file_path: &Path) -> std::io::Result<()> {
    let mut needs_header = false;

    if !file_path.exists() {
        println!("needs a header/file creation");
        needs_header = true;
    } else {
        let metadata = std::fs::metadata(file_path)?;
        if metadata.len() == 0 {
            println!("needs a header/file creation");
            needs_header = true;
        }
    }

    if needs_header {
        let header = "read\treverse_read\tnimble_score_names\tnimble_score_weighted\tforward_filter_reason\tmapq\torientation\tpair\trev_comp\thit\tqname_ascii\ttx\tumi\tcb\tan\tnimble_rev_score_names\tnimble_rev_score_weighted\tnimble_rev_filter_reason\tmapq2\torientation2\tpair2\trev_comp2\thit2\tqname_ascii2\ttx2\talignment_type\treference_orientation\tumi\tcb\tan";

        println!("Opening tsv file for header writing");
        let mut file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(file_path)?;

        println!("Writing header to file");
        writeln!(file, "{}", header)?;
        println!("Finished writing header to file")
    }

    Ok(())
}
