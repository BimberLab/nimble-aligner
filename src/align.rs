use crate::filter;
use crate::reference_library;
use crate::utils::shannon_entropy;

use std::collections::HashMap;
use std::fmt;
use std::io::Error;

use std::default::Default;

use std::fmt::{Display, Formatter};

use array_tool::vec::Intersect;
use array_tool::vec::Uniq;
use debruijn::dna_string::DnaString;
use lexical_sort::{natural_lexical_cmp, StringSort};
use reference_library::ReferenceMetadata;

const MIN_READ_LENGTH: usize = 12;
const MIN_ENTROPY: f64 = 1.75;

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
    HighEntropy,
    SuccessfulMatch,
    StrandWasWrong,
    TriageEmptyEquivalenceClass,
    None,
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
            FilterReason::NotMatchingPair => write!(f, "Required Valid Pair Not Matching"),
            FilterReason::ForceIntersectFailure => write!(f, "Force Intersect Failure"),
            FilterReason::ShortRead => write!(f, "Short Read"),
            FilterReason::MaxHitsExceeded => write!(f, "Max Hits Exceeded"),
            FilterReason::HighEntropy => write!(f, "Low Entropy"),
            FilterReason::SuccessfulMatch => write!(f, "Successful Match"),
            FilterReason::StrandWasWrong => write!(f, "Strandedness Filtered"),
            FilterReason::TriageEmptyEquivalenceClass => write!(f, "Equivalence Class Empty After Filters"),
            FilterReason::None => write!(f, "None"),
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
    pub low_entropy: f64,
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
    None
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
            AlignmentDirection::None => write!(f, "None"),
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
        };

        dir
    }

    fn filter_hits(
        mut forward_hits: Option<(
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        )>,
        mut reverse_hits: Option<(
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        )>,
        results: &mut HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)>,
        reference_metadata: &ReferenceMetadata,
        config: &AlignFilterConfig,
        debug_info: &mut AlignDebugInfo,
        matched_sequences: &mut Vec<(Vec<String>, String, f64, usize, String)>,
        map_key: String,
        debug: bool,
        mut filtered_keys: &mut HashMap<String, (FilterReason, AlignmentDirection)>
    ) {
        if let Some((f_pair_state, _, _, _, _)) = forward_hits {
            if let Some((r_pair_state, _, _, _, _)) = reverse_hits {
                let dir =
                    AlignmentDirection::get_alignment_dir(f_pair_state, r_pair_state, debug_info);
                if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                    filtered_keys.insert(map_key.clone(), (FilterReason::StrandWasWrong, dir));
                    return;
                }
            }

            let dir =
                AlignmentDirection::get_alignment_dir(f_pair_state, PairState::None, debug_info);
            if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                filtered_keys.insert(map_key.clone(), (FilterReason::StrandWasWrong, dir));
                return;
            }
        } else if let Some((r_pair_state, _, _, _, _)) = reverse_hits {
            let dir =
                AlignmentDirection::get_alignment_dir(PairState::None, r_pair_state, debug_info);
            if AlignmentDirection::filter_read(dir, &config.strand_filter) {
                filtered_keys.insert(map_key.clone(), (FilterReason::StrandWasWrong, dir));
                return;
            }
        }

        // Take the "best" alignment. The specific behavior is determined by the intersect level set in the aligner config
        let (match_eqv_class, f_bam_data, r_bam_data) = match config.intersect_level {
            IntersectLevel::NoIntersect => get_all_reads(forward_hits, reverse_hits),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(forward_hits, reverse_hits, true, debug_info, map_key.clone(), &mut filtered_keys)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(forward_hits, reverse_hits, false, debug_info, map_key.clone(), &mut filtered_keys)
            }
        };

        let mut key = get_score_map_key(&match_eqv_class, reference_metadata, &config); // Process the equivalence class into a score key
        key.string_sort_unstable(natural_lexical_cmp); // Sort for deterministic names

        // Max hits filter
        if key.len() > config.max_hits_to_report {
            filtered_keys.insert(map_key.clone(), (FilterReason::MaxHitsExceeded, AlignmentDirection::None));
            return;
        }

        // Ensure we don't add empty keys, if any got to this point
        if key.len() == 0 {
            filtered_keys.insert(map_key.clone(), (FilterReason::TriageEmptyEquivalenceClass, AlignmentDirection::None));
            return;
        }

        let (accessor, f, r) =
            results
                .entry(key)
                .or_insert((0, Vec::new(), Vec::new()));
        *accessor += 1;
        *f = f_bam_data;
        *r = r_bam_data;
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
            AlignmentDirection::None => true,
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
            AlignmentDirection::None => true,
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
            AlignmentDirection::None => true,
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
    current_metadata_group: &Vec<Vec<String>>,
    index_pair: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, (FilterReason, FilterReason, FilterReason, FilterReason, FilterReason, AlignmentDirection)>
) {
    let mut filter_reasons: HashMap<String, (FilterReason, FilterReason, FilterReason, FilterReason, FilterReason, AlignmentDirection)> = HashMap::new();
    let mut filter_reasons_forward: HashMap<String, (FilterReason, FilterReason)> = HashMap::new();
    let mut filter_reasons_reverse: HashMap<String, (FilterReason, FilterReason)> = HashMap::new();
    let mut post_triaged_keys: HashMap<String, (FilterReason, AlignmentDirection)> = HashMap::new();

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
            &mut filter_reasons_forward
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
            &mut filter_reasons_reverse
        );

    forward_matched_sequences.append(&mut backward_matched_sequences);

    let mut results: HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)> = HashMap::new();

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
            &mut post_triaged_keys
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
            &mut post_triaged_keys
        );
    }

    for (key, value) in filter_reasons_forward.into_iter() {
        let reverse_value = filter_reasons_reverse.get(&key).unwrap();
        match post_triaged_keys.get(&key) {
            Some(triage) => filter_reasons.insert(key, (value.0, value.1, reverse_value.0, reverse_value.1, triage.0, triage.1)),
            None => filter_reasons.insert(key, (value.0, value.1, reverse_value.0, reverse_value.1, FilterReason::None, AlignmentDirection::None))
        };
    }

    let mut ret = Vec::new();
    for (key, value) in results.into_iter() {
        ret.push((key, value));
    }

    (ret, forward_matched_sequences, filter_reasons)
}

fn generate_score<'a>(
    sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    mut reverse_sequences: Option<Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>>,
    current_metadata_group: &Vec<Vec<String>>,
    index: &PseudoAligner,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
    reference_orientation: String,
    filter_reasons: &mut HashMap<String, (FilterReason, FilterReason)>
) -> (
    HashMap<
        String,
        (
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        ),
    >,
    Vec<(Vec<String>, String, f64, usize, String)>,
    AlignDebugInfo,
) {
    // HashMap of the alignment results. The keys are either strong hits or equivalence classes of hits
    let mut score_map: HashMap<
        String,
        (
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        ),
    > = HashMap::new();
    let mut debug_info: AlignDebugInfo = Default::default();
    let mut read_matches: Vec<(Vec<String>, String, f64, usize, String)> = Vec::new();

    // Iterate over every read/reverse read pair and align it, incrementing scores for the matching references/equivalence classes
    let mut metadata_iter = current_metadata_group.iter();
    for read in sequences {
        let forward_metadata = metadata_iter.next().unwrap_or(&Vec::new()).clone(); //r2
        let reverse_metadata = metadata_iter.next().unwrap_or(&Vec::new()).clone(); //r1

        let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
        let mut read_rev: Option<DnaString> = None;

        /* Generate score and equivalence class for this read by aligning the sequence against
         * the current reference, if there is a match.*/
        let (seq_score, forward_filter_reason) = pseudoalign(&read, index, &config);

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
        
        let read_key = match &read_rev {
            Some(rev) => read.to_string() + &rev.to_string(),
            None => read.to_string(),
        };

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

        let (failed_score, failed_raw_score) = if reverse_sequences.is_some() {
            match (forward_filter_reason, rev_filter_reason) {
                (Some((fr, s, ns)), Some((rr, r, nr))) => {
                    if fr == rr {
                        (s, ns)
                        //debug_info.update(Some((fr, s, ns)))
                    } else if (fr == FilterReason::NoMatch
                        && rr == FilterReason::ScoreBelowThreshold)
                        || (rr == FilterReason::NoMatch && fr == FilterReason::ScoreBelowThreshold)
                    {
                        let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                        /*debug_info.update(Some((
                            FilterReason::NoMatchAndScoreBelowThreshold,
                            s,
                            ns,
                        )))*/
                        (s, ns)
                    } else {
                        let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                        (s, ns)
                        //debug_info.update(Some((FilterReason::DifferentFilterReasons, s, ns)))
                    }
                }
                (None, Some((rr, r, nr))) => (r, nr),//debug_info.update(Some((rr, r, nr))),
                (Some((fr, s, ns)), None) => (s, ns),//debug_info.update(Some((fr, s, ns))),
                (None, None) => (0.0, 0),
            }
        } else {
            //debug_info.update(forward_filter_reason)
            match forward_filter_reason {
                Some(r) => (r.1, r.2),
                None => (0.0, 0)
            } 
        };

        // If there are no reverse sequences, ignore the require_valid_pair filter
        if reverse_sequences.is_some()
            && config.require_valid_pair
            && filter_pair(&seq_score, &rev_seq_score)
        {
            filter_reasons.insert(read_key.clone(), (FilterReason::NotMatchingPair, FilterReason::NotMatchingPair));
            continue;
        } else {
            filter_reasons.insert(read_key.clone(), (
                match forward_filter_reason {
                    Some((reason, _, _)) => reason,
                    None => FilterReason::SuccessfulMatch
                },
                match rev_filter_reason {
                    Some((reason, _, _)) => reason,
                    None => FilterReason::SuccessfulMatch
                }
            ));
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
                        forward_metadata,
                        reverse_metadata,
                    ),
                    s,
                    r_s,
                    n_s,
                    r_n_s,
                )
            } else if !eqv_class.is_empty() {
                (
                    (
                        PairState::First,
                        Some((eqv_class, s)),
                        None,
                        forward_metadata,
                        reverse_metadata,
                    ),
                    s,
                    0.0,
                    n_s,
                    0,
                )
            } else if !r_eqv_class.is_empty() {
                (
                    (
                        PairState::Second,
                        None,
                        Some((r_eqv_class, r_s)),
                        forward_metadata,
                        reverse_metadata,
                    ),
                    0.0,
                    r_s,
                    0,
                    r_n_s,
                )
            } else {
                continue;
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
    }

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
    seq_score: Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    rev_seq_score: Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    fallback_on_intersect_fail: bool,
    debug_info: &mut AlignDebugInfo,
    map_key: String,
    filtered_keys: &mut HashMap<String, (FilterReason, AlignmentDirection)>
) -> (Vec<u32>, Vec<String>, Vec<String>) {
    let (seq_class, _, seq_f_data, seq_r_data) = match seq_score {
        Some((
            _,
            Some((ref f_unowned, ref fs)),
            Some((ref r_unowned, ref rs)),
            ref f_data,
            ref r_data,
        )) => {
            let mut f = f_unowned.clone();
            let mut r = r_unowned.clone();
            f.append(&mut r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
                f_data.clone(),
                r_data.clone(),
            )
        }
        Some((_, None, Some((ref r, ref rs)), ref f_data, ref r_data)) => {
            (r.to_owned(), rs.to_owned(), f_data.clone(), r_data.clone())
        }
        Some((_, Some((ref f, ref fs)), None, ref f_data, ref r_data)) => {
            (f.to_owned(), fs.to_owned(), f_data.clone(), r_data.clone())
        }
        _ => (Vec::new(), 0.0, Vec::new(), Vec::new()),
    };

    let (r_seq_class, _, rev_seq_f_data, rev_seq_r_data) = match rev_seq_score {
        Some((_, Some((ref f_unowned, fs)), Some((ref r_unowned, rs)), ref f_data, ref r_data)) => {
            let mut f = f_unowned.clone();
            let mut r = r_unowned.clone();
            f.append(&mut r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
                f_data.clone(),
                r_data.clone(),
            )
        }
        Some((_, None, Some((ref r, ref rs)), ref f_data, ref r_data)) => {
            (r.to_owned(), rs.to_owned(), f_data.clone(), r_data.clone())
        }
        Some((_, Some((ref f, ref fs)), None, ref f_data, ref r_data)) => {
            (f.to_owned(), fs.to_owned(), f_data.clone(), r_data.clone())
        }
        _ => (Vec::new(), 0.0, Vec::new(), Vec::new()),
    };

    let class = seq_class.intersect(r_seq_class);

    if class.len() == 0 && fallback_on_intersect_fail {
        get_all_reads(seq_score, rev_seq_score)
    } else if class.len() != 0 {
        if seq_f_data.is_empty() && seq_r_data.is_empty() {
            (class, rev_seq_f_data, rev_seq_r_data)
        } else {
            (class, seq_f_data, seq_r_data)
        }
    } else {
        filtered_keys.insert(map_key, (FilterReason::ForceIntersectFailure, AlignmentDirection::None));
        (Vec::new(), seq_f_data, seq_r_data)
    }
}

// Return best equivalence class out of the given classes
fn get_all_reads(
    seq_score: Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    rev_seq_score: Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
) -> (Vec<u32>, Vec<String>, Vec<String>) {
    let (mut seq_class, _, seq_f_data, seq_r_data) = match seq_score {
        Some((_, Some((ref f_unowned, fs)), Some((ref r_unowned, rs)), f_data, r_data)) => {
            let mut f = f_unowned.clone();
            let mut r = r_unowned.clone();
            f.append(&mut r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
                f_data,
                r_data,
            )
        }
        Some((_, None, Some((r, rs)), f_data, r_data)) => {
            (r.to_owned(), rs.to_owned(), f_data, r_data)
        }
        Some((_, Some((f, fs)), None, f_data, r_data)) => {
            (f.to_owned(), fs.to_owned(), f_data, r_data)
        }
        Some((_, None, None, f_data, r_data)) => (Vec::new(), 0.0, f_data, r_data),
        _ => (Vec::new(), 0.0, Vec::new(), Vec::new()),
    };

    let (mut r_seq_class, _, rev_seq_f_data, rev_seq_r_data) = match rev_seq_score {
        Some((_, Some((ref f_unowned, fs)), Some((ref r_unowned, rs)), f_data, r_data)) => {
            let mut f = f_unowned.clone();
            let mut r = r_unowned.clone();
            f.append(&mut r);
            (
                f.to_owned(),
                if fs > rs {
                    fs.to_owned()
                } else {
                    rs.to_owned()
                },
                f_data,
                r_data,
            )
        }
        Some((_, None, Some((r, rs)), f_data, r_data)) => {
            (r.to_owned(), rs.to_owned(), f_data, r_data)
        }
        Some((_, Some((f, fs)), None, f_data, r_data)) => {
            (f.to_owned(), fs.to_owned(), f_data, r_data)
        }
        Some((_, None, None, f_data, r_data)) => (Vec::new(), 0.0, f_data, r_data),
        _ => (Vec::new(), 0.0, Vec::new(), Vec::new()),
    };

    seq_class.append(&mut r_seq_class);
    seq_class.unique();

    if seq_f_data.is_empty() && seq_r_data.is_empty() {
        (seq_class, rev_seq_f_data, rev_seq_r_data)
    } else {
        (seq_class, seq_f_data, seq_r_data)
    }
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
    config: &AlignFilterConfig
) -> (
    Option<(Vec<u32>, f64, usize)>,
    Option<(FilterReason, f64, usize)>,
) {
    // Filter short reads
    if sequence.to_string().len() < MIN_READ_LENGTH {
        return (None, Some((FilterReason::ShortRead, 0.0, 0)));
    };

    // Filter high-entropy
    if shannon_entropy(&sequence.to_string()) < MIN_ENTROPY {
        return (None, Some((FilterReason::HighEntropy, 0.0, 0)))
    }

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
