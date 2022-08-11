use crate::filter;
use crate::reference_library;

use std::collections::HashMap;
use std::io::Error;

use array_tool::vec::Intersect;
use debruijn::dna_string::DnaString;
use reference_library::ReferenceMetadata;
use lexical_sort::{StringSort, natural_lexical_cmp};
use crate::utils::revcomp;

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
    pub max_hits_to_report: usize,
    pub strand_filter: StrandFilter
}

pub enum StrandFilter {
   Unstranded,
   FivePrime,
   ThreePrime,
   None
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

#[derive(Debug)]
pub enum AlignmentDirection {
    FF,
    RR,
    UU,
    FR,
    FU,
    RF,
    RU,
    UF,
    UR
}

impl AlignmentDirection {
    fn get_alignment_dir(forward_pair_state: PairState, reverse_pair_state: PairState) -> AlignmentDirection {
        match (forward_pair_state, reverse_pair_state) {
            (PairState::First, PairState::First) => AlignmentDirection::FF,
            (PairState::First, PairState::Second) => AlignmentDirection::FR,
            (PairState::Second, PairState::First) => AlignmentDirection::RF,
            (PairState::Second, PairState::Second) => AlignmentDirection::RR,
            (PairState::First, PairState::None) => AlignmentDirection::FU,
            (PairState::Second, PairState::None) => AlignmentDirection::RU,
            (PairState::None, PairState::First) => AlignmentDirection::UF,
            (PairState::None, PairState::Second) => AlignmentDirection::UR,
            (PairState::Both, PairState::First) => AlignmentDirection::FF,
            (PairState::First, PairState::Both) => AlignmentDirection::FF,
            (PairState::Both, PairState::Second) => AlignmentDirection::RR,
            (PairState::Second, PairState::Both) => AlignmentDirection::RR,
            (PairState::Both, PairState::Both) => AlignmentDirection::FR, // TODO: how to deal with these remaining cases across library types
            (PairState::Both, PairState::None) => AlignmentDirection::FU,
            (PairState::None, PairState::Both) => AlignmentDirection::UF,
            (PairState::None, PairState::None) => AlignmentDirection::UU
        }
    }
    
    fn filter_hits(mut forward_hits: Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>, mut reverse_hits: Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>,
                    results: &mut HashMap<Vec<String>, i32>, reference_metadata: &ReferenceMetadata, config: &AlignFilterConfig) {
        if let Some((f_pair_state, _, _)) = forward_hits {
            if let Some((r_pair_state, _, _)) = reverse_hits {
                if AlignmentDirection::filter_read(AlignmentDirection::get_alignment_dir(f_pair_state, r_pair_state), &config.strand_filter) {
                    return
                }
            }

           if AlignmentDirection::filter_read(AlignmentDirection::get_alignment_dir(f_pair_state, PairState::None), &config.strand_filter) {
                return
           }
        } else if let Some((r_pair_state, _, _)) = reverse_hits {
           if AlignmentDirection::filter_read(AlignmentDirection::get_alignment_dir(PairState::None, r_pair_state), &config.strand_filter) {
                return
           }
        }

        // Take the "best" alignment. The specific behavior is determined by the intersect level set in the aligner config
        let match_eqv_class = match config.intersect_level {
            IntersectLevel::NoIntersect => get_best_reads(&mut forward_hits, &mut reverse_hits),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(&mut forward_hits, &mut reverse_hits, true)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(&mut forward_hits, &mut reverse_hits, false)
            }
        };

        let mut key = get_score_map_key(&match_eqv_class, reference_metadata, &config); // Process the equivalence class into a score key
        key.string_sort_unstable(natural_lexical_cmp); // Sort for deterministic names


        // Max hits filter
        if key.len() > config.max_hits_to_report {
           return 
        }

        // Ensure we don't add empty keys, if any got to this point
        if key.len() == 0 {
           return 
        }
           
        let accessor = results.entry(key).or_insert(0);
       *accessor += 1;
    }

    fn filter_read(dir: AlignmentDirection, lib_type: &StrandFilter) -> bool {
        match lib_type {
            StrandFilter::Unstranded => AlignmentDirection::filter_unstranded(dir),
            StrandFilter::FivePrime => AlignmentDirection::filter_fiveprime(dir),
            StrandFilter::ThreePrime => AlignmentDirection::filter_threeprime(dir),
            StrandFilter::None => false
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
            AlignmentDirection::UR => false
        }
    }

    fn filter_fiveprime(dir: AlignmentDirection) -> bool {
        match dir {
            AlignmentDirection::FF => true,
            AlignmentDirection::RR => true,
            AlignmentDirection::UU => true,
            AlignmentDirection::FR => false,
            AlignmentDirection::FU => false,
            AlignmentDirection::RF => true,
            AlignmentDirection::RU => true,
            AlignmentDirection::UF => false,
            AlignmentDirection::UR => false
        }
    }

    fn filter_threeprime(dir: AlignmentDirection) -> bool {
        match dir {
            AlignmentDirection::FF => true,
            AlignmentDirection::RR => true,
            AlignmentDirection::UU => true,
            AlignmentDirection::FR => true,
            AlignmentDirection::FU => false,
            AlignmentDirection::RF => false,
            AlignmentDirection::RU => false,
            AlignmentDirection::UF => false,
            AlignmentDirection::UR => true
        }
    }
}

#[derive(PartialEq, Debug, Clone, Copy)]
pub enum PairState {
    First,
    Second,
    Both,
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

    let (forward_score, forward_matched_sequences, forward_align_debug_info) = generate_score(
        sequences,
        reverse_sequences,
        index_forward,
        reference_metadata,
        config,
    );
    let (mut backward_score, backward_matched_sequences, backward_align_debug_info) = generate_score(
        sequences_2,
        reverse_sequences_2,
        index_backward,
        reference_metadata,
        config,
    );

    let mut results: HashMap<Vec<String>, i32> = HashMap::new();

    for (key, f) in forward_score.into_iter() {
        let r = match backward_score.get(&key) {
            Some(_) => backward_score.remove(&key),
            None => None
        };

        AlignmentDirection::filter_hits(Some(f), r, &mut results, reference_metadata, config);
    };

    for (_, r) in backward_score.into_iter() {
        AlignmentDirection::filter_hits(None, Some(r), &mut results, reference_metadata, config);
    }

    // Update the results map
    let mut ret = Vec::new();
    for (key, value) in results.into_iter() {
        ret.push((key, value));
    }

    (ret, forward_matched_sequences)
}

fn generate_score<'a>(
    sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    mut reverse_sequences: Option<Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>>,
    index: &PseudoAligner,
    reference_metadata: &ReferenceMetadata,
    config: &AlignFilterConfig,
) -> (HashMap<String, (PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>, Vec<(Vec<String>, String, usize)>, AlignDebugInfo) {
    // HashMap of the alignment results. The keys are either strong hits or equivalence classes of hits
    let mut score_map: HashMap<String, (PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)> = HashMap::new();
    let mut debug_info: AlignDebugInfo = Default::default();
    let mut read_matches: Vec<(Vec<String>, String, usize)> = Vec::new();

    // Iterate over every read/reverse read pair and align it, incrementing scores for the matching references/equivalence classes
    for read in sequences {
        let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
        let mut read_rev = None;
        
        /* Generate score and equivalence class for this read by aligning the sequence against
         * the current reference, if there is a match.*/
        let (seq_score, forward_filter_reason) = pseudoalign(&read, index, &config, false);

        // If there's a reversed sequence, do the paired-end alignment
        let mut rev_seq_score = None;
        let mut rev_filter_reason: Option<FilterReason> = None;
        if let Some(itr) = &mut reverse_sequences {
            let reverse_read = itr
                .next()
                .expect("Error -- read and reverse read files do not have matching lengths: ")
                .expect("Error -- could not parse reverse read. Input R2 data malformed.");
            let (score, reason) = pseudoalign(&reverse_read, index, &config, true);
            rev_seq_score = Some(score);
            read_rev = Some(reverse_read);
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

        let (eqv_class, s) = if let Some((eqv_class, s)) = &seq_score {
            ((*eqv_class).clone(), *s)
        } else {
            (Vec::new(), 0)
        };

        let (r_eqv_class, r_s) = if let Some(Some((r_eqv_class, r_s))) = &rev_seq_score {
            ((*r_eqv_class).clone(), *r_s)
        } else {
            (Vec::new(), 0)
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

            let (v, s, r_s) = if !eqv_class.is_empty() && !r_eqv_class.is_empty() {
                ((PairState::Both, Some((eqv_class, s)), Some((r_eqv_class, r_s))), s, r_s)
            } else if !eqv_class.is_empty() {
                ((PairState::First, Some((eqv_class, s)), None), s, 0)
            } else if !r_eqv_class.is_empty() {
                ((PairState::Second, None, Some((r_eqv_class, r_s))), 0, r_s)
            } else {
                continue
            };

            match v.0 {
                PairState::First => read_matches.push((key.clone(), read.to_string(), s)),
                PairState::Second => match &read_rev {
                    Some(r) => read_matches.push((key.clone(), r.to_string(), r_s)),
                    None => ()
                },
                PairState::Both => read_matches.push((key.clone(), read.to_string(), s)),
                PairState::None => ()
            };

            let read_key = match &read_rev {
                Some(rev) => read.to_string() + &rev.to_string(),
                None => read.to_string()
            };

            score_map.insert(read_key, v);
            debug_info.read_units_aligned += 1;
        }
    }

    (score_map, read_matches, debug_info)
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

// Return matches that match in all of the given classes; if soft intersection is enabled, fall back to best read if one of the reads is empty
fn get_intersecting_reads(
    seq_score: &mut Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>,
    rev_seq_score: &mut Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>,
    fallback_on_intersect_fail: bool
) -> Vec<u32> {

    let (seq_class, _) = match seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => { f.append(r); (f.to_owned(), if fs > rs { fs.to_owned() } else { rs.to_owned() })},
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0),
    };

    let (r_seq_class, _) = match rev_seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => { f.append(r); (f.to_owned(), if fs > rs { fs.to_owned() } else { rs.to_owned() })},
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0),
    };

    let class = seq_class.intersect(r_seq_class);

    if class.len() == 0 && fallback_on_intersect_fail {
        get_best_reads(seq_score, rev_seq_score)
    } else if class.len() != 0 {
        class
    } else {
        Vec::new()
    }
}

// Return best equivalence class out of the given classes
fn get_best_reads(
    seq_score: &mut Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>,
    rev_seq_score: &mut Option<(PairState, Option<(Vec<u32>, usize)>, Option<(Vec<u32>, usize)>)>,
) -> Vec<u32> {
    let (seq_class, seq_score) = match seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => { f.append(r); (f.to_owned(), if fs > rs { fs.to_owned() } else { rs.to_owned() })},
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0),
    };

    let (r_seq_class, r_seq_score) = match rev_seq_score {
        Some((_, Some((ref mut f, fs)), Some((ref mut r, rs)))) => { f.append(r); (f.to_owned(), if fs > rs { fs.to_owned() } else { rs.to_owned() })},
        Some((_, None, Some((r, rs)))) => (r.to_owned(), rs.to_owned()),
        Some((_, Some((f, fs)), None)) => (f.to_owned(), fs.to_owned()),
        _ => (Vec::new(), 0),
    };

    if seq_score > r_seq_score { seq_class.to_owned() } else { r_seq_class.to_owned() } 
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
            let mut group = &reference_metadata.columns[reference_metadata.group_on][*ref_idx as usize];

            // If the group is an empty string, get the sequence name instead so we don't just return nothing
            if group.is_empty() {
                group = &reference_metadata.columns[reference_metadata.sequence_name_idx][*ref_idx as usize];
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
    reverse_comp: bool
) -> (Option<(Vec<u32>, usize)>, Option<FilterReason>){
    // Filter short reads
    if sequence.to_string().len() < MIN_READ_LENGTH {
        return (None, Some(FilterReason::ShortRead))
    };

    // TODO figure out reverse_comp of second in pair
    /*let sequence = if reverse_comp {
        DnaString::from_dna_string(&revcomp(&sequence.to_string())).to_owned()
    } else {
        sequence.to_owned()
    };*/

    // Perform alignment
    match reference_index.map_read_with_mismatch(&sequence, config.num_mismatches) {
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
