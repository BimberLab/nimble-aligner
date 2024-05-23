use crate::filter;
use crate::reference_library;
use crate::utils::shannon_entropy;

use std::collections::HashMap;
use std::fmt;
use std::io::Error;

use std::fmt::{Display, Formatter};

use array_tool::vec::Intersect;
use array_tool::vec::Uniq;
use debruijn::dna_string::DnaString;
use lexical_sort::{natural_lexical_cmp, StringSort};
use reference_library::Reference;

const MIN_READ_LENGTH: usize = 12;
const MIN_ENTROPY_SCORE: f64 = 1.75; // Higher score = lower entropy

pub type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::Kmer30>;
pub type AlignmentScore = Option<(Vec<u32>, f64, usize)>;
pub type Filter = Option<(FilterReason, f64, usize)>;

#[derive(Debug, PartialEq)]
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
    pub strand_filter: LibraryChemistry,
}

#[derive(Debug, Copy, Clone)]
pub enum LibraryChemistry {
    Unstranded,
    FivePrime,
    ThreePrime,
    None,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum AlignmentOrientation {
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

impl fmt::Display for AlignmentOrientation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentOrientation::FF => write!(f, "FF"),
            AlignmentOrientation::RR => write!(f, "RR"),
            AlignmentOrientation::UU => write!(f, "UU"),
            AlignmentOrientation::FR => write!(f, "FR"),
            AlignmentOrientation::FU => write!(f, "FU"),
            AlignmentOrientation::RF => write!(f, "RF"),
            AlignmentOrientation::RU => write!(f, "RU"),
            AlignmentOrientation::UF => write!(f, "UF"),
            AlignmentOrientation::UR => write!(f, "UR"),
            AlignmentOrientation::None => write!(f, "None"),
        }
    }
}

/* Get the orientation of a set of alignments for a single read-pair. If there's a call of the same feature for both mates,
 * but with opposite orientations, we return FR or RF -- the strongest evidence of a successful mapping, given that the read-mates are
 * sequenced from opposite ends of the input molecule.
 * Single hits with one mate get FU, UF, RU, or UR respectively. A feature call of both orientations for a single mate return FF or RR,
 * which is generally a sign of an errant mapping.
 * Finally, if both mates hit a single orientation, with one mate hitting the other orientation, or both mates hit both orientations,
 * we return FF or RR, which generally get filtered, or UU which is considered invalid. */
impl AlignmentOrientation {
    fn get_alignment_orientation(
        forward_pair_state: PairState,
        reverse_pair_state: PairState,
    ) -> AlignmentOrientation {
        let dir = match (forward_pair_state, reverse_pair_state) {
            (PairState::First, PairState::First) => AlignmentOrientation::FF,
            (PairState::First, PairState::Second) => AlignmentOrientation::FR,
            (PairState::Second, PairState::First) => AlignmentOrientation::RF,
            (PairState::Second, PairState::Second) => AlignmentOrientation::RR,
            (PairState::First, PairState::None) => AlignmentOrientation::FU,
            (PairState::Second, PairState::None) => AlignmentOrientation::RU,
            (PairState::None, PairState::First) => AlignmentOrientation::UF,
            (PairState::None, PairState::Second) => AlignmentOrientation::UR,

            (PairState::Both, PairState::First) => AlignmentOrientation::FF,
            (PairState::First, PairState::Both) => AlignmentOrientation::FF,
            (PairState::Both, PairState::Second) => AlignmentOrientation::RR,
            (PairState::Second, PairState::Both) => AlignmentOrientation::RR,
            (PairState::Both, PairState::Both) => AlignmentOrientation::UU,
            (PairState::Both, PairState::None) => AlignmentOrientation::UU,
            (PairState::None, PairState::Both) => AlignmentOrientation::UU,
            (PairState::None, PairState::None) => AlignmentOrientation::UU,
        };

        dir
    }

    
    /* The value in sequence_call contains a pair_state that reports whether one, both, or neither mates in the given read-pair aligned
     * to the given feature. sequence_call_revcomp_library also contains this information, but for the revcomp version of the feature.
     * As such, we can compute a "combined orientation" for all 4 alignments (both read-mates for both versions of the feature),
     * which we use to filter out calls generated by alignments in a "suspect" combined orientation. The specific set of orientations
     * considered suspect -- i.e., should be filtered, are determined based on the library chemistry flag in the aligner configuration. */
    fn filter_and_coerce_sequence_call_orientations(
        sequence_call: Option<(
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        )>,
        sequence_call_revcomp_library: Option<(
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        )>,
        results: &mut HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)>,
        reference_metadata: &Reference,
        config: &AlignFilterConfig,
        read_key: String,
        mut filtered_keys: &mut HashMap<String, (FilterReason, AlignmentOrientation)>
    ) {
        if let Some((pair_state, _, _, _, _)) = sequence_call {
            // If we passed both sequence_calls and sequence_calls_revcomp_library, determine whether or not to filter the call based on both pair_states
            if let Some((pair_state_revcomp_library, _, _, _, _)) = sequence_call_revcomp_library {
                let alignment_orientation =
                    AlignmentOrientation::get_alignment_orientation(pair_state, pair_state_revcomp_library);

                if AlignmentOrientation::filter_orientation_on_library_chemistry(alignment_orientation, &config.strand_filter) {
                    filtered_keys.insert(read_key.clone(), (FilterReason::StrandWasWrong, alignment_orientation));
                    return;
                }

            // Otherwise, filter the call based on only the sequence_calls pair state
            } else {
                let alignment_orientation =
                    AlignmentOrientation::get_alignment_orientation(pair_state, PairState::None);

                if AlignmentOrientation::filter_orientation_on_library_chemistry(alignment_orientation, &config.strand_filter) {
                    filtered_keys.insert(read_key.clone(), (FilterReason::StrandWasWrong, alignment_orientation));
                    return;
                }
            }
        // If there are no sequence_calls, filter based only on the pair state in sequence_call_revcomp_library
        } else if let Some((pair_state_revcomp_library, _, _, _, _)) = sequence_call_revcomp_library {
            let alignment_orientation =
                AlignmentOrientation::get_alignment_orientation(PairState::None, pair_state_revcomp_library);

            if AlignmentOrientation::filter_orientation_on_library_chemistry(alignment_orientation, &config.strand_filter) {
                filtered_keys.insert(read_key.clone(), (FilterReason::StrandWasWrong, alignment_orientation));
                return;
            }
        }

        /* We have up to four different potential calls for a given read-pair, which are coerced into a single callset here.
         * There are a few different options:
         *  - No Intersect: merge all the calls together across normal and revcomped versions of the reference features
         *  - Intersect With Fallback: intersect calls between the orientations, but fall back to No Intersect if the equivalence class is destroyed
         *  - Force Intersect: intersect calls between the orientations, discarding the read upon intersection failure */
        let (final_callset, sequence_metadata, mate_sequence_metadata) = match config.intersect_level {
            IntersectLevel::NoIntersect => get_all_calls(&sequence_call, &sequence_call_revcomp_library),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(&sequence_call, &sequence_call_revcomp_library, true, read_key.clone(), &mut filtered_keys)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(&sequence_call, &sequence_call_revcomp_library, false, read_key.clone(), &mut filtered_keys)
            }
        };

        let feature_callset = process_equivalence_class_to_feature_list(&final_callset, reference_metadata, &config);

        // Max hits filter. group_by rollup reduces the number of called features, so we do this after process_equivalence_class_to_feature_list
        if feature_callset.len() > config.max_hits_to_report {
            filtered_keys.insert(read_key.clone(), (FilterReason::MaxHitsExceeded, AlignmentOrientation::None));
            return;
        }

        // Ensure we don't add empty keys, if any got to this point
        if feature_callset.len() == 0 {
            filtered_keys.insert(read_key.clone(), (FilterReason::TriageEmptyEquivalenceClass, AlignmentOrientation::None));
            return;
        }

        // Add the call for this read-pair to the final result
        let (count_accessor, sequence_metadata_accessor, mate_sequence_metadata_accessor) =
            results
                .entry(feature_callset)
                .or_insert((0, Vec::new(), Vec::new()));
        *count_accessor += 1;
        *sequence_metadata_accessor = sequence_metadata;
        *mate_sequence_metadata_accessor = mate_sequence_metadata;
    }

    // Given the alignment orientation derived from the set of alignments for a read-pair, determine whether to filter it or not based on the library chemistry
    fn filter_orientation_on_library_chemistry(dir: AlignmentOrientation, lib_type: &LibraryChemistry) -> bool {
        match lib_type {
            LibraryChemistry::Unstranded => AlignmentOrientation::filter_unstranded(dir),
            LibraryChemistry::FivePrime => AlignmentOrientation::filter_fiveprime(dir),
            LibraryChemistry::ThreePrime => AlignmentOrientation::filter_threeprime(dir),
            LibraryChemistry::None => false,
        }
    }

    fn filter_unstranded(dir: AlignmentOrientation) -> bool {
        match dir {
            AlignmentOrientation::FF => true,
            AlignmentOrientation::RR => true,
            AlignmentOrientation::UU => true,
            AlignmentOrientation::FR => false,
            AlignmentOrientation::FU => false,
            AlignmentOrientation::RF => false,
            AlignmentOrientation::RU => false,
            AlignmentOrientation::UF => false,
            AlignmentOrientation::UR => false,
            AlignmentOrientation::None => true,
        }
    }

    fn filter_fiveprime(dir: AlignmentOrientation) -> bool {
        match dir {
            AlignmentOrientation::FF => true,
            AlignmentOrientation::RR => true,
            AlignmentOrientation::UU => true,
            AlignmentOrientation::FR => true,
            AlignmentOrientation::FU => true,
            AlignmentOrientation::RF => false,
            AlignmentOrientation::RU => false,
            AlignmentOrientation::UF => false,
            AlignmentOrientation::UR => true,
            AlignmentOrientation::None => true,
        }
    }

    fn filter_threeprime(dir: AlignmentOrientation) -> bool {
        match dir {
            AlignmentOrientation::FF => true,
            AlignmentOrientation::RR => true,
            AlignmentOrientation::UU => true,
            AlignmentOrientation::FR => false,
            AlignmentOrientation::FU => false,
            AlignmentOrientation::RF => true,
            AlignmentOrientation::RU => true,
            AlignmentOrientation::UF => true,
            AlignmentOrientation::UR => false,
            AlignmentOrientation::None => true,
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

/* Takes a set of sequences and optionally, mate sequences, associated sequence metadata,
 * debrujin indices of the reference library of interest,
 * the associated reference library metadata, and the aligner configuration, and performs a
 * debrujin-graph based pseduoalignment, returning a set of feature calls for each read-mate
 * that passes alignment and filtration. 
 * The score set is filtered based on the provided aligner configuration. */
pub fn get_calls<'a>(
    sequence_iterators: (
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    ),
    mate_sequence_iterators: Option<(
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    )>,
    sequence_metadata: &Vec<Vec<String>>,
    indices: &(PseudoAligner, PseudoAligner),
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)>
) {
    let mut filter_reasons: HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)> = HashMap::new();
    let mut filter_reasons_forward: HashMap<String, ((FilterReason, usize), (FilterReason, usize))> = HashMap::new();
    let mut filter_reasons_reverse: HashMap<String, ((FilterReason, usize), (FilterReason, usize))> = HashMap::new();
    let mut post_triaged_keys: HashMap<String, (FilterReason, AlignmentOrientation)> = HashMap::new();

    // Unpack the index and sequence pairs
    let (index, index_revcomp) = indices;
    let (sequences, sequences_clone) = sequence_iterators;
    let (mate_sequences, mate_sequences_clone) = match mate_sequence_iterators {
        Some((l, r)) => (Some(l), Some(r)),
        None => (None, None),
    };

    // Generate a set of passing scores for the sequences (mate sequences optional), for both the regular and reverse complemented versions of the reference
    let (sequence_scores, mut matched_sequences) =
        score_sequences(
            sequences,
            mate_sequences,
            sequence_metadata,
            index,
            reference,
            aligner_config,
            &mut filter_reasons_forward
        );
    let (mut sequence_scores_revcomp_library, mut matched_sequences_revcomp_library) =
        score_sequences(
            sequences_clone,
            mate_sequences_clone,
            sequence_metadata,
            index_revcomp,
            reference,
            aligner_config,
            &mut filter_reasons_reverse
        );

    matched_sequences.append(&mut matched_sequences_revcomp_library);

    // Final results hashmap, which will contain all calls that remain after orientation filtering
    let mut results: HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)> = HashMap::new();

    /* Iterate all of the calls generated from aligning the read-pairs to the reference library in the normal orientation.
     * Send those calls to the alignment orientation filtering pipeline.
     * If the same read-pair got a match to the revcomped version of the reference genome, grab that call and also send it to the
     * orientation filter. Once this process is complete, the results hashmap should contain one value per unique feature, and a count of the
     * number of read-pairs that returned a call for that feature. */
    for (read_pair_key, call) in sequence_scores.into_iter() {
        let call_revcomp_library = match sequence_scores_revcomp_library.get(&read_pair_key) {
            Some(_) => sequence_scores_revcomp_library.remove(&read_pair_key),
            None => None,
        };

        AlignmentOrientation::filter_and_coerce_sequence_call_orientations(
            Some(call),
            call_revcomp_library,
            &mut results,
            reference,
            aligner_config,
            read_pair_key,
            &mut post_triaged_keys
        );
    }

    // All sequence_scores have been processed -- process the remaining calls unique to sequence_scores_revcomp_library
    for (key, r) in sequence_scores_revcomp_library.into_iter() {
        AlignmentOrientation::filter_and_coerce_sequence_call_orientations(
            None,
            Some(r),
            &mut results,
            reference,
            aligner_config,
            key,
            &mut post_triaged_keys
        );
    }

    /*  Merge all the read-filtration information into a single collection, including normal/revcomp library
     *   filter information, and the alignment orientation filter pipeline's returned info */
    for (key, value) in filter_reasons_forward.into_iter() {
        let reverse_value = filter_reasons_reverse.get(&key).unwrap();
        match post_triaged_keys.get(&key) {
            Some(triage) => filter_reasons.insert(key, (value.0, value.1, reverse_value.0, reverse_value.1, triage.0, triage.1)),
            None => filter_reasons.insert(key, (value.0, value.1, reverse_value.0, reverse_value.1, FilterReason::None, AlignmentOrientation::None))
        };
    }

    // Iterate the keys and values in the results hashmap, pushing them to the results vector to return
    let mut ret = Vec::new();
    for (key, value) in results.into_iter() {
        ret.push((key, value));
    }

    (ret, matched_sequences, filter_reasons)
}


/* 
  * The core alignment function.
  * Takes a list of sequences, and optionally reverse sequences, and produces a map of those reads/read-pairs and their corresponding feature calls
  * based on a given reference library and aligner configuration
  */
fn score_sequences<'a>(
    sequences: Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    mut mate_sequences: Option<Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>>,
    sequence_metadata: &Vec<Vec<String>>,
    index: &PseudoAligner,
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
    filter_reasons: &mut HashMap<String, ((FilterReason, usize), (FilterReason, usize))>
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
    Vec<(Vec<String>, String, f64, usize, String)>
) {
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
    let mut read_matches: Vec<(Vec<String>, String, f64, usize, String)> = Vec::new();
    let mut metadata_iter = sequence_metadata.iter();

    // Iterate over every read/read pair and align it, adding the results to a hashmap of feature calls
    for read in sequences {

        // For each sequence, get the according metadata object
        // The fastq pipeline cannot pass a metadata object, so we assume paired-end organization for the metadata iterator
        let sequence_metadata = metadata_iter.next().unwrap_or(&Vec::new()).clone();
        let mate_sequence_metadata = metadata_iter.next().unwrap_or(&Vec::new()).clone();

        let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
        let mut read_rev: Option<DnaString> = None;

        // Generate score and equivalence class for this read by aligning the sequence against the reference.
        let (sequence_alignment, sequence_filter_reason) = pseudoalign(&read, index, &aligner_config);

        // If there's a mate sequence, also perform the alignment for it
        let mut mate_sequence_alignment: Option<AlignmentScore> = None;
        let mut mate_sequence_filter_reason: Filter = None;

        if let Some(itr) = &mut mate_sequences {
            let reverse_read = itr
                .next()
                .expect("Error -- read and reverse read files do not have matching lengths: ")
                .expect("Error -- could not parse reverse read. Input R2 data malformed.");

            let (score, filter_reason) = pseudoalign(&reverse_read, index, &aligner_config);

            mate_sequence_alignment = Some(score);
            read_rev = Some(reverse_read);
            mate_sequence_filter_reason = filter_reason;
        }

        // Unpack the alignments for each sequence into separate variables
        let (sequence_equivalence_class, normalized_sequence_score, sequence_score) = if let Some((eqv_class, n_s, s)) = &sequence_alignment {
            ((*eqv_class).clone(), *n_s, *s)
        } else {
            (Vec::new(), 0.0, 0)
        };

        let (mate_sequence_equivalence_class, normalized_mate_sequence_score, mate_sequence_score) = if let Some(Some((eqv_class, n_s, s))) = &mate_sequence_alignment
        {
            ((*eqv_class).clone(), *n_s, *s)
        } else {
            (Vec::new(), 0.0, 0)
        };

        // Create a key for the read-pair using the sequence data, for registration in the output score hashmap
        let read_key = match &read_rev {
            Some(rev) => read.to_string() + &rev.to_string(),
            None => read.to_string(),
        };

        // Check if there are mate sequences. If there are, and require_valid_pair is enabled, we filter the pair if both alignments are not identical
        if mate_sequences.is_some()
            && aligner_config.require_valid_pair
            && filter_pair(&sequence_equivalence_class, &mate_sequence_equivalence_class)
        {
            filter_reasons.insert(read_key.clone(), ((FilterReason::NotMatchingPair, sequence_score), (FilterReason::NotMatchingPair, mate_sequence_score)));

            continue;
        } else {
            // Otherwise, this was either a successful match or it was filtered at a pseudoalignment level. Either way, we add it to the filter report
            filter_reasons.insert(read_key.clone(), (
                match sequence_filter_reason {
                    Some((reason, _, _)) => (reason, sequence_score),
                    None => (FilterReason::SuccessfulMatch, sequence_score)
                },
                match mate_sequence_filter_reason {
                    Some((reason, _, _)) => (reason, mate_sequence_score),
                    None => (FilterReason::SuccessfulMatch, mate_sequence_score)
                }
            ));
        }

        // At this point, if there are unfiltered alignments, we convert the equivalence classes to feature lists by looking them up in the reference
        if !sequence_equivalence_class.is_empty() || !mate_sequence_equivalence_class.is_empty() {
            let feature_list = if !sequence_equivalence_class.is_empty() {
                process_equivalence_class_to_feature_list(&sequence_equivalence_class, reference, &aligner_config)
            } else if !mate_sequence_equivalence_class.is_empty() {
                process_equivalence_class_to_feature_list(&mate_sequence_equivalence_class, reference, &aligner_config)
            } else {
                Vec::new()
            };

            // Package the equivalence classes and metadata into a result value for the score map
            let (pair_score,
                normalized_sequence_score,
                normalized_mate_sequence_score,
                sequence_score,
                mate_sequence_score) = if !sequence_equivalence_class.is_empty() && !mate_sequence_equivalence_class.is_empty() {
                (
                    (
                        PairState::Both,
                        Some((sequence_equivalence_class, normalized_sequence_score)),
                        Some((mate_sequence_equivalence_class, normalized_mate_sequence_score)),
                        sequence_metadata,
                        mate_sequence_metadata,
                    ),
                    normalized_sequence_score,
                    normalized_mate_sequence_score,
                    sequence_score,
                    mate_sequence_score,
                )
            } else if !sequence_equivalence_class.is_empty() {
                (
                    (
                        PairState::First,
                        Some((sequence_equivalence_class, normalized_sequence_score)),
                        None,
                        sequence_metadata,
                        mate_sequence_metadata,
                    ),
                    normalized_sequence_score,
                    0.0,
                    sequence_score,
                    0,
                )
            } else if !mate_sequence_equivalence_class.is_empty() {
                (
                    (
                        PairState::Second,
                        None,
                        Some((mate_sequence_equivalence_class, normalized_mate_sequence_score)),
                        sequence_metadata,
                        mate_sequence_metadata,
                    ),
                    0.0,
                    normalized_mate_sequence_score,
                    0,
                    mate_sequence_score,
                )
            } else {
                continue;
            };

            // At this point, we have a match and add it to read_matches and score_map, which are this function's exports
            match pair_score.0 {
                PairState::First => {
                    read_matches.push((feature_list.clone(), read.to_string(), normalized_sequence_score, sequence_score, read_key.clone()))
                }
                PairState::Second => match &read_rev {
                    Some(r) => read_matches.push((
                        feature_list.clone(),
                        r.to_string(),
                        normalized_mate_sequence_score,
                        mate_sequence_score,
                        read_key.clone(),
                    )),
                    None => (),
                },
                PairState::Both => {
                    read_matches.push((feature_list.clone(), read.to_string(), normalized_sequence_score, sequence_score, read_key.clone()))
                }
                PairState::None => (),
            };
            
            score_map.insert(read_key, pair_score);
        } else {
            // If both equivalence classes are empty, the attempted alignment has failed, but we still report the failed alignment using read_matches
            let (failed_score, failed_raw_score) = if mate_sequences.is_some() {
                match (sequence_filter_reason, mate_sequence_filter_reason) {
                    (Some((filter_reason, s, ns)), Some((mate_filter_reason, r, nr))) => {
                        if filter_reason == mate_filter_reason {
                            (s, ns)
                        } else if (filter_reason == FilterReason::NoMatch
                            && mate_filter_reason == FilterReason::ScoreBelowThreshold)
                            || (mate_filter_reason == FilterReason::NoMatch && filter_reason == FilterReason::ScoreBelowThreshold)
                        {
                            let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                            (s, ns)
                        } else {
                            let (s, ns) = if s > r { (s, ns) } else { (r, nr) };

                            (s, ns)
                        }
                    }
                    (None, Some((_rr, r, nr))) => (r, nr),
                    (Some((_fr, s, ns)), None) => (s, ns),
                    (None, None) => (0.0, 0),
                }
            } else {
                // If there was a match to the first sequence but it got filtered, return the scores. Otherwise return 0 (no match)
                match sequence_filter_reason {
                    Some(r) => (r.1, r.2),
                    None => (0.0, 0)
                } 
            };

            read_matches.push((
                Vec::new(),
                read.to_string(),
                failed_score,
                failed_raw_score,
                String::new(),
            ));
        }
    }

    (score_map, read_matches)
}

// Determine whether a given pair of equivalence classes constitute a valid alignment for a read pair
fn filter_pair(
    sequence_equivalence_class: &Vec<u32>,
    mate_sequence_equivalence_class: &Vec<u32>,
) -> bool {
    let mut sequence_equivalence_class = sequence_equivalence_class.clone();
    let mut mate_sequence_equivalence_class = mate_sequence_equivalence_class.clone();

    // Check if the pairs have identical equivalence classes. If they both contain data, do the comparison
    if !sequence_equivalence_class.is_empty() && !mate_sequence_equivalence_class.is_empty() 
    {
        // Sort the vectors and compare them by counting matching elements
        mate_sequence_equivalence_class.sort();
        sequence_equivalence_class.sort();
        let matching = sequence_equivalence_class
            .iter()
            .zip(&mate_sequence_equivalence_class)
            .filter(|&(classl, classr)| classl == classr)
            .count();

        if matching != sequence_equivalence_class.len() || matching != mate_sequence_equivalence_class.len() {
            return true;
        }
    } else {
        // Since this function is only called when require_valid_pair is on, if one doesn't contain data, the equivalence classes are not identical
        // and this pair should be filtered. Or, if they both don't contain data, it should be filtered regardless.
        return true;
    }
    false
}

// For two given sets of calls, intersect them to get a final callset for the read-pair. Optionally, fall back to a more permissive strategy if intersection fails
fn get_intersecting_reads(
    sequence_call: &Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    sequence_call_revcomp_library: &Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    fallback_on_intersect_fail: bool,
    read_key: String,
    filtered_keys: &mut HashMap<String, (FilterReason, AlignmentOrientation)>
) -> (Vec<u32>, Vec<String>, Vec<String>) {

    // Coerce the callsets for both sets of alignments
    let (sequence_equivalence_class, _, sequence_metadata, mate_sequence_metadata) = coerce_sequence_callset(sequence_call);
    let (sequence_equivalence_class_revcomp_library, _, sequence_metadata_revcomp_library, mate_sequence_metadata_revcomp_library) = coerce_sequence_callset(sequence_call_revcomp_library);

    // Get only the calls shared between normal and revcomped alignments
    let class = sequence_equivalence_class.intersect(sequence_equivalence_class_revcomp_library);

    // If that destroyed the alignment and fallback is on, use get_call_calls instead
    if class.len() == 0 && fallback_on_intersect_fail {
        get_all_calls(sequence_call, sequence_call_revcomp_library)
    
    // Otherwise, if the intersection succeeded, return the class with any existing metadata
    } else if class.len() != 0 {
        if sequence_metadata.is_empty() && mate_sequence_metadata.is_empty() {
            (class, sequence_metadata_revcomp_library, mate_sequence_metadata_revcomp_library)
        } else {
            (class, sequence_metadata, mate_sequence_metadata)
        }

    // If it failed and fallback is off, add this read-pair to the filtered list
    } else {
        filtered_keys.insert(read_key, (FilterReason::ForceIntersectFailure, AlignmentOrientation::None));
        (Vec::new(), sequence_metadata, mate_sequence_metadata)
    }
}

// Permissively append all the calls for a given set of sequence calls, creating one large combined equivalence class
fn get_all_calls(
    sequence_call: &Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
    sequence_call_revcomp_library: &Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
) -> (Vec<u32>, Vec<String>, Vec<String>) {
    let (mut sequence_equivalence_class, _, sequence_metadata, mate_sequence_metadata) = coerce_sequence_callset(sequence_call);
    let (mut sequence_equivalence_class_revcomp_library, _, sequence_metadata_revcomp_library, mate_sequence_metadata_revcomp_library) = coerce_sequence_callset(sequence_call_revcomp_library);

    // Take the two coerced equivalence classes and merge them
    sequence_equivalence_class.append(&mut sequence_equivalence_class_revcomp_library);
    sequence_equivalence_class.unique();

    // Return them and any existing metadata objects
    if sequence_metadata.is_empty() && mate_sequence_metadata.is_empty() {
        (sequence_equivalence_class, sequence_metadata_revcomp_library, mate_sequence_metadata_revcomp_library)
    } else {
        (sequence_equivalence_class, sequence_metadata, mate_sequence_metadata)
    }
}

// Given a set of calls for a read-pair, coerce them down to a single equivalence class by appending and return that new class + the read metadata
fn coerce_sequence_callset(
    sequence_call: &Option<(
        PairState,
        Option<(Vec<u32>, f64)>,
        Option<(Vec<u32>, f64)>,
        Vec<String>,
        Vec<String>,
    )>,
) -> (Vec<u32>, f64, Vec<String>, Vec<String>) {
    // Take a given sequence call and match on the calls
    match sequence_call {
        // If there are calls for both read-mates, we need to coerce them down to a single equivalence class and score
        Some((_, Some((ref sequence_equivalence_class, normalized_score)),
                 Some((ref mate_sequence_equivalence_class, mate_normalized_score)),
                 sequence_metadata, mate_sequence_metadata)) => {
            
            // Append the equivalence classes, and return whichever normalized score is higher between the two
            let mut sequence_equivalence_class = sequence_equivalence_class.clone();
            let mut mate_sequence_equivalence_class = mate_sequence_equivalence_class.clone();
            sequence_equivalence_class.append(&mut mate_sequence_equivalence_class);
            (
                sequence_equivalence_class.to_owned(),
                if normalized_score > mate_normalized_score {
                    normalized_score.to_owned()
                } else {
                    mate_normalized_score.to_owned()
                },
                sequence_metadata.clone(),
                mate_sequence_metadata.clone(),
            )
        }
        // Otherwise, if there is only one score, return that score. Or if there are none, return none.
        Some((_, None, Some((mate_sequence_equivalence_class, mate_normalized_score)), sequence_metadata, mate_sequence_metadata)) => {
            (mate_sequence_equivalence_class.to_owned(), mate_normalized_score.to_owned(), sequence_metadata.clone(), mate_sequence_metadata.clone())
        }
        Some((_, Some((sequence_equivalence_class, sequence_normalized_score)), None, sequence_metadata, mate_sequence_metadata)) => {
            (sequence_equivalence_class.to_owned(), sequence_normalized_score.to_owned(), sequence_metadata.clone(), mate_sequence_metadata.clone())
        }
        Some((_, None, None, sequence_metadata, mate_sequence_metadata)) => (Vec::new(), 0.0, sequence_metadata.clone(), mate_sequence_metadata.clone()),
        _ => (Vec::new(), 0.0, Vec::new(), Vec::new()),
    }
}

/* Takes an equivalence class and returns a list of strings. If we're processing allele-level data, the strings will be
 * the sequence name of the relevant feature(s). Otherwise, if we're doing a group_by, the equivalence class will be
 * filtered such that there is only one hit per group_by string (e.g. one hit per lineage) and the corresponding strings
 * (e.g. lineage name) will be returned instead of individual sequence names. */
fn process_equivalence_class_to_feature_list(
    equivalence_class: &Vec<u32>,
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
) -> Vec<String> {
    /* If the group_on column is the default, nt_sequence, no feature roll-up occurs. We map over the equivalence class and translate one-to-one from
     * values in that column. */
    let mut results = if reference.headers[reference.group_on] == "nt_sequence" {
        equivalence_class
            .into_iter()
            .map(|ref_idx| {
                reference.columns[reference.group_on][*ref_idx as usize].clone()
            })
            .collect()
    } else {
        let mut grouped_feature_list = Vec::new();

        // Otherwise, we iterate the equivalence class, and grab the value from the group_on column instead
        for ref_idx in equivalence_class {
            let mut group =
                &reference.columns[reference.group_on][*ref_idx as usize];

            // If the group_on value is empty, we fall back to the feature name column, as if group_on == 'nt_sequence'
            if group.is_empty() {
                group = &reference.columns[reference.sequence_name_idx]
                    [*ref_idx as usize];
            }

            /* Then, we push to the results array, and only push unique values. Therefore, features with idential group_on values will only
             * contribute to the feature list once. */
            if !grouped_feature_list.contains(group) {
                grouped_feature_list.push(group.to_string());
            }
        }

        grouped_feature_list 
    };

    // Finally, if the library has set a discard_multi_hits filter and the resultant feature list is larger, drop this alignment
    if aligner_config.discard_multi_hits > 0 && results.len() > aligner_config.discard_multi_hits {
        Vec::new()
    } else {
        // Sort for deterministic feature lists
        results.string_sort_unstable(natural_lexical_cmp);
        results
    }
}

// Align a sequence against the reference library, with settings specified by the aligner configuration
fn pseudoalign(
    sequence: &DnaString,
    reference_index: &PseudoAligner,
    config: &AlignFilterConfig
) -> (
    AlignmentScore,
    Filter
) {
    // Ensure all reads are above a minimum length
    if sequence.len() < MIN_READ_LENGTH {
        return (None, Some((FilterReason::ShortRead, 0.0, 0)));
    };

    // Remove reads that do not pass an entropy threshold 
    if shannon_entropy(&sequence.to_string()) < MIN_ENTROPY_SCORE {
        return (None, Some((FilterReason::HighEntropy, 0.0, 0)))
    }

    // Perform the alignment of the sequence to the reference debruijn graph. Pass the number of allowed mismatches
    match reference_index.map_read_with_mismatch(&sequence, config.num_mismatches) {
        Some((equivalence_class, score, mismatches)) => {

            // Normalize score by read length
            let normalized_score = score as f64 / sequence.len() as f64;

            // Filter nonzero mismatch
            if config.discard_nonzero_mismatch && mismatches != 0 {
                return (None, Some((FilterReason::DiscardedNonzeroMismatch, 0.0, 0)));
            }

            // Filter by score thresholds 
            filter::align::filter_alignment_by_metrics(
                equivalence_class,
                score,
                normalized_score,
                config.score_threshold,
                config.score_percent,
                config.discard_multiple_matches,
            )
        }
        None => (None, Some((FilterReason::NoMatch, 0.0, 0))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use debruijn::kmer::Kmer30;

    fn setup_pseudoaligner() -> PseudoAligner {
        let reference_sequences = vec![
            DnaString::from_dna_string("ACGTACGTACGTACGTACGTACGTACGTACGT"),
            DnaString::from_dna_string("TGCATGCATGCATGCATGCATGCATGCATGCA"),
        ];
        let reference_feature_names = vec![
            String::from("Gene1"),
            String::from("Gene2"),
        ];

        debruijn_mapping::build_index::build_index::<Kmer30>(
            &reference_sequences,
            &reference_feature_names,
            &HashMap::new(),
            1
        ).expect("Failed to build index")
    }

    fn setup_reference() -> Reference {
        Reference {
            group_on: 0,
            headers: vec!["nt_sequence".to_string(), "gene".to_string()],
            columns: vec![
                vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()],
                vec!["geneA".to_string(), "geneB".to_string(), "geneA".to_string()]
            ],
            sequence_name_idx: 0,
            sequence_idx: 0,
        }
    }

    fn setup_config() -> AlignFilterConfig {
        AlignFilterConfig {
            reference_genome_size: 1000,
            score_percent: 0.1,
            score_threshold: 50,
            num_mismatches: 3,
            discard_nonzero_mismatch: false,
            discard_multiple_matches: false,
            score_filter: 10,
            intersect_level: IntersectLevel::IntersectWithFallback,
            require_valid_pair: false,
            discard_multi_hits: 0,
            max_hits_to_report: 5,
            strand_filter: LibraryChemistry::FivePrime,
        }
    }

    #[test]
    fn test_short_read() {
        let sequence = DnaString::from_dna_string("ACG");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config);
        assert_eq!(result.1, Some((FilterReason::ShortRead, 0.0, 0)));
    }

    #[test]
    fn test_high_entropy_read() {
        let sequence = DnaString::from_dna_string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config);
        assert_eq!(result.1, Some((FilterReason::HighEntropy, 0.0, 0)));
    }

    #[test]
    fn test_no_alignment_match() {
        let sequence = DnaString::from_dna_string("CCTGAGATTTCGAGCTCGTAACGTGACCTACGGACAC");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config);
        assert_eq!(result.1, Some((FilterReason::NoMatch, 0.0, 0)));
    }

    #[test]
    fn test_valid_alignment() {
        let sequence = DnaString::from_dna_string("TGCATGCATGCATGCATGCATGCATGCATGCA");
        let mut config = setup_config();
        let reference_index = setup_pseudoaligner();
        config.score_threshold = 32;
        let result = pseudoalign(&sequence, &reference_index, &config);
        assert_eq!(result.0, Some((vec![1], 1.0, 32)));
        assert_eq!(result.1, None);
    }

    #[test]
    fn test_score_and_match_threshold_filtering() {
        let sequence = DnaString::from_dna_string("TGCATGCATGCATGCATGCATGCATGCATGCA");
        let mut config = setup_config();
        config.score_threshold = 1000;
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config);
        assert_eq!(result.1, Some((FilterReason::ScoreBelowThreshold, 1.0, 32)));
    }

    #[test]
    fn test_empty_equivalence_classes() {
        let seq_classes = vec![];
        let mate_classes = vec![];
        assert_eq!(filter_pair(&seq_classes, &mate_classes), true);
    }

    #[test]
    fn test_one_empty_one_non_empty() {
        let non_empty_classes = vec![1, 2, 3];
        let empty_classes = vec![];
        assert_eq!(filter_pair(&non_empty_classes, &empty_classes), true);
        assert_eq!(filter_pair(&empty_classes, &non_empty_classes), true);
    }

    #[test]
    fn test_different_non_empty_classes() {
        let seq_classes = vec![1, 2, 3];
        let mate_classes = vec![4, 5, 6];
        assert_eq!(filter_pair(&seq_classes, &mate_classes), true);
    }

    #[test]
    fn test_identical_non_empty_classes() {
        let seq_classes = vec![1, 2, 3];
        let mate_classes = vec![1, 2, 3];
        assert_eq!(filter_pair(&seq_classes, &mate_classes), false);
    }

    #[test]
    fn test_different_lengths_some_matching_elements() {
        let seq_classes = vec![1, 2, 3, 4];
        let mate_classes = vec![1, 2, 3];
        assert_eq!(filter_pair(&seq_classes, &mate_classes), true);
    }

    #[test]
    fn test_group_by_nt_sequence() {
        let ref_data = setup_reference();
        let config = setup_config();
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config),
            vec!["seq1", "seq2", "seq3"]
        );
    }

    #[test]
    fn test_group_by_gene() {
        let mut ref_data = setup_reference();
        ref_data.group_on = 1;
        let config = setup_config();
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config),
            vec!["geneA", "geneB"]
        );
    }

    #[test]
    fn test_fallback_to_feature_name() {
        let mut ref_data = setup_reference();
        ref_data.columns[1] = vec!["geneA".to_string(), "".to_string(), "geneA".to_string()];
        ref_data.group_on = 1;
        let config = setup_config();
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config),
            vec!["geneA", "seq2"]
        );
    }

    #[test]
    fn test_discard_multi_hits() {
        let ref_data = setup_reference();
        let mut config = setup_config();
        config.discard_multi_hits = 1;
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config),
            Vec::<String>::new()
        );
    }

    #[test]
    fn test_empty_equivalence_class() {
        let ref_data = setup_reference();
        let config = setup_config();
        let equivalence_class = vec![];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config),
            Vec::<String>::new()
        );
    }

    #[test]
    fn test_list_stability_and_order() {
        let mut ref_data = setup_reference();
        ref_data.group_on = 1;
        let config = setup_config();

        let equivalence_class_1 = vec![2, 0, 1];
        let equivalence_class_2 = vec![0, 1, 2];

        let result_1 = process_equivalence_class_to_feature_list(&equivalence_class_1, &ref_data, &config);
        let result_2 = process_equivalence_class_to_feature_list(&equivalence_class_2, &ref_data, &config);

        assert_eq!(result_1, result_2);
        assert_eq!(result_1, vec!["geneA", "geneB"]);
    }

    #[test]
    fn test_alignment_orientation() {
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::First, PairState::First), AlignmentOrientation::FF);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::First, PairState::Second), AlignmentOrientation::FR);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Second, PairState::First), AlignmentOrientation::RF);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Second, PairState::Second), AlignmentOrientation::RR);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::First, PairState::None), AlignmentOrientation::FU);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Second, PairState::None), AlignmentOrientation::RU);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::None, PairState::First), AlignmentOrientation::UF);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::None, PairState::Second), AlignmentOrientation::UR);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Both, PairState::First), AlignmentOrientation::FF);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::First, PairState::Both), AlignmentOrientation::FF);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Both, PairState::Second), AlignmentOrientation::RR);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Second, PairState::Both), AlignmentOrientation::RR);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Both, PairState::Both), AlignmentOrientation::UU);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::Both, PairState::None), AlignmentOrientation::UU);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::None, PairState::Both), AlignmentOrientation::UU);
        assert_eq!(AlignmentOrientation::get_alignment_orientation(PairState::None, PairState::None), AlignmentOrientation::UU);
    }

    #[test]
    fn test_filter_unstranded() {
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::FF), true);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::RR), true);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::UU), true);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::FR), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::FU), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::RF), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::RU), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::UF), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::UR), false);
        assert_eq!(AlignmentOrientation::filter_unstranded(AlignmentOrientation::None), true);
    }

    #[test]
    fn test_filter_fiveprime() {
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::FF), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::RR), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::UU), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::FR), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::FU), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::RF), false);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::RU), false);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::UF), false);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::UR), true);
        assert_eq!(AlignmentOrientation::filter_fiveprime(AlignmentOrientation::None), true);
    }

    #[test]
    fn test_filter_threeprime() {
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::FF), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::RR), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::UU), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::FR), false);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::FU), false);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::RF), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::RU), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::UF), true);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::UR), false);
        assert_eq!(AlignmentOrientation::filter_threeprime(AlignmentOrientation::None), true);
    }

    #[test]
    fn test_filter_orientation_on_library_chemistry() {
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FF, &LibraryChemistry::Unstranded), true);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FR, &LibraryChemistry::Unstranded), false);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FF, &LibraryChemistry::FivePrime), true);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::RF, &LibraryChemistry::FivePrime), false);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FF, &LibraryChemistry::ThreePrime), true);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FR, &LibraryChemistry::ThreePrime), false);
        assert_eq!(AlignmentOrientation::filter_orientation_on_library_chemistry(AlignmentOrientation::FF, &LibraryChemistry::None), false);
    }

    #[test]
    fn test_coerce_sequence_callset_both_present() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta1".to_string()],
            vec!["mate_meta1".to_string()],
        ));
        let result = coerce_sequence_callset(&sequence_call);
        assert_eq!(result, (vec![1, 2, 3, 4, 5, 6], 0.9, vec!["seq_meta1".to_string()], vec!["mate_meta1".to_string()]));
    }

    #[test]
    fn test_coerce_sequence_callset_only_mate() {
        let sequence_call = Some((
            PairState::None,
            None,
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let result = coerce_sequence_callset(&sequence_call);
        assert_eq!(result, (vec![4, 5, 6], 0.8, vec!["seq_meta2".to_string()], vec!["mate_meta2".to_string()]));
    }

    #[test]
    fn test_coerce_sequence_callset_only_sequence() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            None,
            vec!["seq_meta3".to_string()],
            vec!["mate_meta3".to_string()],
        ));
        let result = coerce_sequence_callset(&sequence_call);
        assert_eq!(result, (vec![1, 2, 3], 0.9, vec!["seq_meta3".to_string()], vec!["mate_meta3".to_string()]));
    }

    #[test]
    fn test_coerce_sequence_callset_none_present() {
        let sequence_call = Some((
            PairState::None,
            None,
            None,
            vec!["seq_meta4".to_string()],
            vec!["mate_meta4".to_string()],
        ));
        let result = coerce_sequence_callset(&sequence_call);
        assert_eq!(result, (vec![], 0.0, vec!["seq_meta4".to_string()], vec!["mate_meta4".to_string()]));
    }

    #[test]
    fn test_coerce_sequence_callset_none() {
        let result = coerce_sequence_callset(&None);
        assert_eq!(result, (vec![], 0.0, vec![], vec![]));
    }

    #[test]
    fn test_get_all_calls_both_present() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta1".to_string()],
            vec!["mate_meta1".to_string()],
        ));
        let sequence_call_revcomp_library = Some((
            PairState::None,
            Some((vec![7, 8, 9], 0.7)),
            Some((vec![10, 11, 12], 0.6)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let result = get_all_calls(&sequence_call, &sequence_call_revcomp_library);
        assert_eq!(result, (vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], vec!["seq_meta1".to_string()], vec!["mate_meta1".to_string()]));
    }

    #[test]
    fn test_get_all_calls_empty_metadata() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec![],
            vec![],
        ));
        let sequence_call_revcomp_library = Some((
            PairState::None,
            Some((vec![7, 8, 9], 0.7)),
            Some((vec![10, 11, 12], 0.6)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let result = get_all_calls(&sequence_call, &sequence_call_revcomp_library);
        assert_eq!(result, (vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], vec!["seq_meta2".to_string()], vec!["mate_meta2".to_string()]));
    }

    #[test]
    fn test_get_intersecting_reads_intersect_success() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3, 4], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta1".to_string()],
            vec!["mate_meta1".to_string()],
        ));
        let sequence_call_revcomp_library = Some((
            PairState::None,
            Some((vec![4, 5, 6, 7], 0.7)),
            Some((vec![8, 9, 10], 0.6)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(&sequence_call, &sequence_call_revcomp_library, false, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, (vec![4, 5, 6], vec!["seq_meta1".to_string()], vec!["mate_meta1".to_string()]));
        assert!(filtered_keys.is_empty());
    }

    #[test]
    fn test_get_intersecting_reads_fallback_on_intersect_fail() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta1".to_string()],
            vec!["mate_meta1".to_string()],
        ));
        let sequence_call_revcomp_library = Some((
            PairState::None,
            Some((vec![7, 8, 9], 0.7)),
            Some((vec![10, 11, 12], 0.6)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(&sequence_call, &sequence_call_revcomp_library, true, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, (vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], vec!["seq_meta1".to_string()], vec!["mate_meta1".to_string()]));
        assert!(filtered_keys.is_empty());
    }

    #[test]
    fn test_get_intersecting_reads_intersect_fail_no_fallback() {
        let sequence_call = Some((
            PairState::None,
            Some((vec![1, 2, 3], 0.9)),
            Some((vec![4, 5, 6], 0.8)),
            vec!["seq_meta1".to_string()],
            vec!["mate_meta1".to_string()],
        ));
        let sequence_call_revcomp_library = Some((
            PairState::None,
            Some((vec![7, 8, 9], 0.7)),
            Some((vec![10, 11, 12], 0.6)),
            vec!["seq_meta2".to_string()],
            vec!["mate_meta2".to_string()],
        ));
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(&sequence_call, &sequence_call_revcomp_library, false, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, (vec![], vec!["seq_meta1".to_string()], vec!["mate_meta1".to_string()]));
        assert_eq!(filtered_keys.get("read_key"), Some(&(FilterReason::ForceIntersectFailure, AlignmentOrientation::None)));
    }
}