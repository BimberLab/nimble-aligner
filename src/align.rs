use crate::filter;
use crate::reference_library;

use std::collections::HashMap;
use std::io::Error;

use array_tool::vec::Intersect;
use debruijn::dna_string::DnaString;
use reference_library::ReferenceMetadata;
use lexical_sort::{StringSort, natural_lexical_cmp};


const MIN_READ_LENGTH: usize = 12;
const SCORE_BASELINE: usize = 20;


pub type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<
    debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>,
>;

pub enum IntersectLevel {
    NoIntersect,
    IntersectWithFallback,
    ForceIntersect,
}

#[derive(Debug, PartialEq)]
pub enum FilterReason {
    ScoreBelowThreshold,
    DiscardedMultipleMatch,
    DiscardedNonzeroMismatch,
    NoMatch,
    NoMatchAndScoreBelowThreshold,
    DifferentFilterReasons,
    NotMatchingPair,
    ForceIntersectFailure,
    DisjointPairIntersection,
    BestClassEmpty,
    ShortRead
}

pub struct AlignFilterConfig {
    pub reference_genome_size: usize,
    pub score_threshold: usize,
    pub num_mismatches: usize,
    pub discard_nonzero_mismatch: bool,
    pub discard_multiple_matches: bool,
    pub score_filter: i32,
    pub intersect_level: IntersectLevel,
    pub require_valid_pair: bool,
    pub discard_multi_hits: usize,
    pub max_hits_to_report: usize
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
    pub disjoint_pair_intersection: usize,
    pub best_class_empty: usize,
    pub forward_runs_discarded: usize,
    pub backward_runs_discarded: usize,
    pub reverse_read_sets_discarded_noneven: usize,
    pub short_read: usize
}

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
    I
}

pub enum LibraryType {
   Unstranded,
   FivePrime,
   ThreePrime 
}

impl AlignmentDirection {
    fn get_alignment_dir(forward_pair_state: PairState, reverse_pair_state: PairState) -> AlignmentDirection {
        match (forward_pair_state, reverse_pair_state) {
            (PairState::First, PairState::First) => AlignmentDirection::FU,
            (PairState::First, PairState::Second) => AlignmentDirection::FR,
            (PairState::Second, PairState::First) => AlignmentDirection::RF,
            (PairState::Second, PairState::Second) => AlignmentDirection::RU,
            (PairState::First, PairState::None) => AlignmentDirection::FU,
            (PairState::Second, PairState::None) => AlignmentDirection::RU,
            (PairState::None, PairState::First) => AlignmentDirection::FU,
            (PairState::None, PairState::Second) => AlignmentDirection::RU,
            (PairState::Intersect, _) => AlignmentDirection::I,
            (_, PairState::Intersect) => AlignmentDirection::I,
            (PairState::None, PairState::None) => AlignmentDirection::UU
        }
    }

    fn filter_read(dir: AlignmentDirection, lib_type: LibraryType) -> bool {
        match lib_type {
            LibraryType::Unstranded => AlignmentDirection::filter_unstranded(dir),
            LibraryType::FivePrime => AlignmentDirection::filter_fiveprime(dir),
            LibraryType::ThreePrime => AlignmentDirection::filter_threeprime(dir)
        }
    }

    fn filter_unstranded(dir: AlignmentDirection) -> bool {
        match dir {
            FF => false,
            RR => true,
            UU => true,
            FR => false,
            FU => false,
            RF => false,
            RU => false,
            UF => false,
            UR => false,
            I => false,
        }
    }

    fn filter_fiveprime(dir: AlignmentDirection) -> bool {
        match dir {
            FF => false,
            RR => true,
            UU => true,
            FR => false,
            FU => false,
            RF => false,
            RU => false,
            UF => false,
            UR => false,
            I => false,
        }
    }

    fn filter_threeprime(dir: AlignmentDirection) -> bool {
        match dir {
            FF => false,
            RR => true,
            UU => true,
            FR => false,
            FU => false,
            RF => false,
            RU => false,
            UF => false,
            UR => false,
            I => false,
        }
    }
}

pub enum PairState {
    First,
    Second,
    Intersect,
    None
}

impl AlignDebugInfo {
    fn update(&mut self, reason: Option<FilterReason>) {
        match reason {
            Some(FilterReason::ScoreBelowThreshold) => self.score_below_threshold += 1,
            Some(FilterReason::DiscardedMultipleMatch) => self.discarded_multiple_match += 1,
            Some(FilterReason::DiscardedNonzeroMismatch) => self.discarded_nonzero_mismatch += 1,
            Some(FilterReason::NoMatch) => self.no_match += 1,
            Some(FilterReason::NoMatchAndScoreBelowThreshold) => self.no_match_and_score_below_threshold += 1,
            Some(FilterReason::DifferentFilterReasons) => self.different_filter_reasons += 1,
            Some(FilterReason::NotMatchingPair) => self.not_matching_pair += 1,
            Some(FilterReason::ForceIntersectFailure) => self.force_intersect_failure += 1,
            Some(FilterReason::DisjointPairIntersection) => self.disjoint_pair_intersection += 1,
            Some(FilterReason::BestClassEmpty) => self.best_class_empty += 1,
            Some(FilterReason::ShortRead) => self.short_read += 1,
            None => (),
        }
    }

    fn merge(&mut self, info: AlignDebugInfo) {
        self.read_units_aligned += info.read_units_aligned;
        self.read_units_aligned += info.read_units_aligned;
        self.score_below_threshold += info.score_below_threshold;
        self.discarded_multiple_match += info.discarded_multiple_match;
        self.discarded_nonzero_mismatch += info.discarded_nonzero_mismatch;
        self.no_match += info.no_match;
        self.no_match_and_score_below_threshold += info.no_match_and_score_below_threshold;
        self.different_filter_reasons += info.different_filter_reasons;
        self.not_matching_pair += info.not_matching_pair;
        self.force_intersect_failure += info.force_intersect_failure;
        self.disjoint_pair_intersection += info.disjoint_pair_intersection;
        self.best_class_empty += info.best_class_empty;
        self.short_read += info.short_read;
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
    index_pair: &(PseudoAligner, PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
    debug_info: Option<&mut AlignDebugInfo>
) -> (Vec<(Vec<String>, i32)>, Vec<(Vec<String>, String, usize)>) {
    let (index_forward, index_backward) = index_pair;
    let (sequences, sequences_2) = sequence_iter_pair;
    let (reverse_sequences, reverse_sequences_2) = match reverse_sequence_iter_pair {
        Some((l, r)) => (Some(l), Some(r)),
        None => (None, None),
    };

    let (forward_score, forward_matched_sequences, forward_align_debug_info, forward_pair_states) = generate_score(
        sequences,
        reverse_sequences,
        index_forward,
        reference_metadata,
        config,
    );
    let (backward_score, backward_matched_sequences, backward_align_debug_info, reverse_pair_states) = generate_score(
        sequences_2,
        reverse_sequences_2,
        index_backward,
        reference_metadata,
        config,
    );

    let filter_list = Vec::new();
    for (match_names, sequence, _) in forward_matched_sequences.iter() {
        let forward_pair_state = match forward_pair_states.get(sequence) {
            Some(state) => state,
            None => &PairState::None 
        };
        
       filter_list.append(AlignmentDirection::get_filter_list((match_names, forward_pair_state));
    };

    for (match_names, sequence, _) in backward_matched_sequences.iter() {
        let reverse_pair_state = match reverse_pair_states.get(sequence) {
            Some(state) => state,
            None => &PairState::None
        };

       filter_list.append(AlignmentDirection::get_filter_list((match_names, forward_pair_state));
    }

    if forward_score.len() > backward_score.len() {
        if let Some(debug_info) = debug_info {
            debug_info.merge(forward_align_debug_info);
            debug_info.backward_runs_discarded += 1;
        }

        (forward_score, forward_matched_sequences)
    } else {
        if let Some(debug_info) = debug_info {
            debug_info.merge(backward_align_debug_info);
            debug_info.forward_runs_discarded += 1;
        }

        (backward_score, backward_matched_sequences)
    }
}

fn generate_score<'a>(
    sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    mut reverse_sequences: Option<Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>>,
    index: &PseudoAligner,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
) -> (Vec<(Vec<String>, i32)>, Vec<(Vec<String>, String, usize)>, AlignDebugInfo, HashMap<String, PairState>) {
    // HashMap of the alignment results. The keys are either strong hits or equivalence classes of hits
    let mut score_map: HashMap<Vec<String>, i32> = HashMap::new();
    let mut debug_info: AlignDebugInfo = Default::default();
    let mut read_matches: Vec<(Vec<String>, String, usize)> = Vec::new();
    let mut pair_states: HashMap<String, PairState> = HashMap::new();

    // Iterate over every read/reverse read pair and align it, incrementing scores for the matching references/equivalence classes
    for read in sequences {
        let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
        
        /* Generate score and equivalence class for this read by aligning the sequence against
         * the current reference, if there is a match.*/
        let (seq_score, forward_filter_reason) = pseudoalign(&read, index, &config);

        // If there's a reversed sequence, do the paired-end alignment
        let mut rev_seq_score = None;
        let mut rev_filter_reason: Option<FilterReason> = None;
        if let Some(itr) = &mut reverse_sequences {
            let reverse_read = itr
                .next()
                .expect("Error -- read and reverse read files do not have matching lengths: ")
                .expect("Error -- could not parse reverse read. Input R2 data malformed.");
            let (score, reason) = pseudoalign(&reverse_read, index, &config);
            rev_seq_score = Some(score);
            rev_filter_reason = reason;
        }

        if reverse_sequences.is_some() {
            match (forward_filter_reason, rev_filter_reason) {
                (Some(fr), Some(rr)) => {
                    if fr == rr {
                        debug_info.update(Some(fr));
                    } else if (fr == FilterReason::NoMatch && rr == FilterReason::ScoreBelowThreshold) ||
                              (rr == FilterReason::NoMatch && fr == FilterReason::ScoreBelowThreshold) { 

                        debug_info.update(Some(FilterReason::NoMatchAndScoreBelowThreshold));
                    } else {
                        debug_info.update(Some(FilterReason::DifferentFilterReasons));
                    }
                },
                (None, Some(rr)) => debug_info.update(Some(rr)),
                (Some(fr), None) => debug_info.update(Some(fr)),
                (None, None) => (),
            };
        } else {
            debug_info.update(forward_filter_reason);
        }

        // If there are no reverse sequences, ignore the require_valid_pair filter
        if reverse_sequences.is_some()
            && config.require_valid_pair
            && filter_pair(&seq_score, &rev_seq_score)
        {
            debug_info.update(Some(FilterReason::NotMatchingPair));
            continue;
        }

        // Take the "best" alignment. The specific behavior is determined by the intersect level set in the aligner config
        let (match_eqv_class, pseudoaligner_score, pair_state) = match config.intersect_level {
            IntersectLevel::NoIntersect => get_best_reads(&seq_score, &rev_seq_score),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(&seq_score, &rev_seq_score, true, &mut debug_info)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(&seq_score, &rev_seq_score, false, &mut debug_info)
            }
        };

        if !match_eqv_class.is_empty() {
            let mut key = get_score_map_key(match_eqv_class, reference_metadata, &config); // Process the equivalence class into a score key
            key.string_sort_unstable(natural_lexical_cmp); // Sort for deterministic names 

            // Max hits filter
            if key.len() > config.max_hits_to_report {
                continue;
            }

            // Ensure we don't add empty keys, if any got to this point
            if key.len() == 0 {
                continue;
            }

            read_matches.push((key.clone(), read.to_string(), pseudoaligner_score));
            pair_states.insert(read.to_string(), pair_state);

            // Add the key to the score map and increment the score
            let accessor = score_map.entry(key).or_insert(0);
            *accessor += 1;
            debug_info.read_units_aligned += 1;
        }
    }

    // Update the results map
    let mut results = Vec::new();
    for (key, value) in score_map.into_iter() {
        results.push((key, value));
    }

    (results, read_matches, debug_info, pair_states)
}

// Determine whether a given pair of equivalence classes constitute a valid pair
fn filter_pair(
    seq_score: &Option<(Vec<u32>, usize)>,
    rev_seq_score: &Option<Option<(Vec<u32>, usize)>>,
) -> bool {
    // Unpack the data to check if the pairs have the same eq_class. If they both contain data, do the comparison
    if let (Some(Some((mut rev_eq_class, _))), Some((mut eq_class, _))) =
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

// Return matches that match in both seq_score and rev_seq_score; if soft intersection is enabled, fall back to best read if one of the reads is empty
fn get_intersecting_reads(
    seq_score: &Option<(Vec<u32>, usize)>,
    rev_seq_score: &Option<Option<(Vec<u32>, usize)>>,
    fallback_on_intersect_fail: bool,
    debug_info: &mut AlignDebugInfo
) -> (Vec<u32>, usize, PairState) {
    if let (Some((eqv_class_seq, score)), Some(Some((eqv_class_rev_seq, _rev_score)))) =
        (&seq_score, &rev_seq_score)
    {
        let class = eqv_class_seq.intersect(eqv_class_rev_seq.to_vec());

        if class.len() == 0 {
            debug_info.update(Some(FilterReason::DisjointPairIntersection))
        }

        (class, *score, PairState::Intersect)
    } else if fallback_on_intersect_fail {
        let (class, score, pair_state) = get_best_reads(seq_score, rev_seq_score);
        
        if class.len() == 0 {
            debug_info.update(Some(FilterReason::BestClassEmpty));
        }
        (class, score, pair_state)
    } else {
        debug_info.update(Some(FilterReason::ForceIntersectFailure));
        (Vec::new(), 0, PairState::None)
    }
}

// Return matches from seq_score -- otherwise, return matches from rev_seq_score
fn get_best_reads(
    seq_score: &Option<(Vec<u32>, usize)>,
    rev_seq_score: &Option<Option<(Vec<u32>, usize)>>,
) -> (Vec<u32>, usize, PairState) {
    if let Some((eqv_class, s)) = &seq_score {
        if let Some(Some((r_eqv_class, r_s))) = &rev_seq_score {
           if s > r_s { return ((*eqv_class).clone(), *s, PairState::First) } else { return ((*r_eqv_class).clone(), *r_s, PairState::Second) }
        }
        return ((*eqv_class).clone(), *s, PairState::First)
    } else if let Some(Some((r_eqv_class, r_s))) = &rev_seq_score {
        return ((*r_eqv_class).clone(), *r_s, PairState::Second)
    }
    (Vec::new(), 0, PairState::None)
}

/* Takes a equivalence class and returns a list of strings. If we're processing allele-level data, the strings will be
 * the nt_sequences of the relevant alleles. Otherwise, if we're doing a group_by, the equivalence class will be
 * filtered such that there is only one hit per group_by string (e.g. one hit per lineage) and the corresponding strings
 * (e.g. lineage name) will be returned. */
fn get_score_map_key(
    equiv_class: Vec<u32>,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
) -> Vec<String> {
    if reference_metadata.headers[reference_metadata.group_on] == "nt_sequence" {
        equiv_class
            .into_iter()
            .map(|ref_idx| {
                reference_metadata.columns[reference_metadata.group_on][ref_idx as usize].clone()
            })
            .collect()
    } else {
        let mut results = Vec::new();

        for ref_idx in equiv_class {
            let mut group = &reference_metadata.columns[reference_metadata.group_on][ref_idx as usize];

            // If the group is an empty string, get the sequence name instead so we don't just return nothing
            if group.is_empty() {
                group = &reference_metadata.columns[reference_metadata.sequence_name_idx][ref_idx as usize];
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
) -> (Option<(Vec<u32>, usize)>, Option<FilterReason>){
    // Filter short reads
    if sequence.to_string().len() < MIN_READ_LENGTH {
        return (None, Some(FilterReason::ShortRead))
    };
        
    // Perform alignment
    match reference_index.map_read_with_mismatch(sequence, config.num_mismatches) {
        Some((equiv_class, score, mismatches)) => {
            let score = score - SCORE_BASELINE;
            // Filter nonzero mismatch
            if config.discard_nonzero_mismatch && mismatches != 0 {
                return (None, Some(FilterReason::DiscardedNonzeroMismatch));
            }

            // Filter by score and match threshold
            filter::align::filter_alignment_by_metrics(
                score,
                equiv_class,
                config.score_threshold,
                config.discard_multiple_matches,
            )
        }
        None => (None, Some(FilterReason::NoMatch)),
    }
}
