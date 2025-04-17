use crate::filter;
use crate::reference_library;
use crate::utils::shannon_entropy;

use std::collections::HashMap;
use std::fmt;
use std::io::Error;

use std::collections::HashSet;
use std::fmt::{Display, Formatter};

use array_tool::vec::Intersect;
use array_tool::vec::Uniq;
use debruijn::dna_string::DnaString;
use lexical_sort::{natural_lexical_cmp, StringSort};
use reference_library::Reference;

const MIN_READ_LENGTH: usize = 40;
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
    AboveMismatchThreshold,
    SkippedAlignDueToUnpairedDummy,
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
            FilterReason::AboveMismatchThreshold => write!(f, "Above Mismatch Threshold"),
            FilterReason::SkippedAlignDueToUnpairedDummy => write!(f, "SKipped Align Due To Unpaired Dummy Read"),
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
    pub trim_strictness: f64,
    pub trim_target_length: usize,
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
    fn filter_read_calls_with_orientation(class: Vec<String>) -> Vec<String> {
        let mut seen = HashSet::new();
        let mut to_remove = HashSet::new();
        
        for feature in &class {
            let base_name = if let Some(base) = feature.strip_suffix(&(reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR.to_string() + "rev")) {
                base
            } else {
               feature 
            };

            if seen.contains(base_name) {
                to_remove.insert(base_name.to_string());
            } else {
                seen.insert(base_name.to_string());
            }
        }

        class.into_iter()
            .filter(|call| {
                if let Some(base_name) = call.strip_suffix(&(reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR.to_string() + "rev")) {
                    !to_remove.contains(base_name)
                } else {
                    !to_remove.contains(call)
                }
            })
            .collect()
    }

    /* The value in sequence_call contains a pair_state that reports whether one, both, or neither mates in the given read-pair aligned
     * to the given feature. sequence_call_revcomp_library also contains this information, but for the revcomp version of the feature.
     * As such, we can compute a "combined orientation" for all 4 alignments (both read-mates for both versions of the feature),
     * which we use to filter out calls generated by alignments in a "suspect" combined orientation. The specific set of orientations
     * considered suspect -- i.e., should be filtered, are determined based on the library chemistry flag in the aligner configuration. */
    fn filter_and_coerce_sequence_call_orientations(
        call: (
            PairState,
            Option<(Vec<u32>, f64)>,
            Option<(Vec<u32>, f64)>,
            Vec<String>,
            Vec<String>,
        ),
        results: &mut HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)>,
        reference_metadata: &Reference,
        config: &AlignFilterConfig,
        read_key: String,
        mut filtered_keys: &mut HashMap<String, (FilterReason, AlignmentOrientation)>
    ) {
        let (_, sequence_call_option, mate_sequence_call_option, sequence_metadata, mate_sequence_metadata) = call;

        // Map the equivalence classes to feature calls so we can parse the revcomped features and begin orientation filtration
        let mut sequence_features = Vec::new();
        if let Some((sequence_call, _)) = sequence_call_option {
            sequence_features = process_equivalence_class_to_feature_list(&sequence_call, reference_metadata, &config, true);
        }

        let mut mate_sequence_features = Vec::new();
        if let Some((mate_sequence_call, _)) = mate_sequence_call_option {
            mate_sequence_features = process_equivalence_class_to_feature_list(&mate_sequence_call, reference_metadata, &config, true);
        }

        // Check within each class for feature + rev feature. If it hits both, remove the hit
        sequence_features = AlignmentOrientation::filter_read_calls_with_orientation(sequence_features);
        mate_sequence_features = AlignmentOrientation::filter_read_calls_with_orientation(mate_sequence_features);

        /* For each hit in the class determine an orientation based on if it's also in the other class
            * There are orientations that are valid, and those that are not. If a particular call doesn't comprise a valid orientation, it is
            * removed here. */
        let (sequence_features, mate_sequence_features) = AlignmentOrientation::filter_orientation_on_library_chemistry(sequence_features, mate_sequence_features, &config.strand_filter);

        /* We have up to four different potential calls for a given read-pair, which are coerced into a single callset here.
         * There are a few different options:
         *  - No Intersect: merge all the calls together across normal and revcomped versions of the reference features
         *  - Intersect With Fallback: intersect calls between the orientations, but fall back to No Intersect if the equivalence class is destroyed
         *  - Force Intersect: intersect calls between the orientations, discarding the read upon intersection failure */
        let final_callset = match config.intersect_level {
            IntersectLevel::NoIntersect => get_all_calls(sequence_features, mate_sequence_features),
            IntersectLevel::IntersectWithFallback => {
                get_intersecting_reads(sequence_features, mate_sequence_features, true, read_key.clone(), &mut filtered_keys)
            }
            IntersectLevel::ForceIntersect => {
                get_intersecting_reads(sequence_features, mate_sequence_features, false, read_key.clone(), &mut filtered_keys)
            }
        };

        let final_callset = unmap(&final_callset, reference_metadata);
        let feature_callset = process_equivalence_class_to_feature_list(&final_callset, reference_metadata, &config, false);

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
    fn filter_orientation_on_library_chemistry(sequence_calls: Vec<String>, mate_sequence_calls: Vec<String>, lib_type: &LibraryChemistry) -> (Vec<String>, Vec<String>) {
        let parsed_sequence_calls = AlignmentOrientation::parse_calls(sequence_calls);
        let parsed_mate_sequence_calls = AlignmentOrientation::parse_calls(mate_sequence_calls);        

        match lib_type {
            LibraryChemistry::None => (
                parsed_sequence_calls.into_iter().map(|(feat, _)| feat).collect(),
                parsed_mate_sequence_calls.into_iter().map(|(feat, _)| feat).collect()
            ),
            LibraryChemistry::Unstranded => {
                let (calls, mate_calls) = AlignmentOrientation::filter_unstranded(parsed_sequence_calls, parsed_mate_sequence_calls);
                (
                    calls.into_iter().map(|(feat, _)| feat).collect(),
                    mate_calls.into_iter().map(|(feat, _)| feat).collect()
                )
            },
            LibraryChemistry::FivePrime => AlignmentOrientation::filter_five_prime(parsed_sequence_calls, parsed_mate_sequence_calls),
            LibraryChemistry::ThreePrime => AlignmentOrientation::filter_three_prime(parsed_sequence_calls, parsed_mate_sequence_calls),
        } 
    }

    fn parse_calls(calls: Vec<String>) -> Vec<(String, bool)> {
        calls.into_iter().map(|call| {
            if call.ends_with("rev") {
                let base_feature = call.trim_end_matches("rev").trim_end_matches(reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR).to_string();
                (base_feature, true)
            } else {
                (call, false)
            }
        }).collect()
    }

    fn filter_unstranded(sequence_calls: Vec<(String, bool)>, mate_sequence_calls: Vec<(String, bool)>) -> (Vec<(String, bool)>, Vec<(String, bool)>) {
        let mut sequence_calls_filtered = vec![];
        let mut mate_sequence_calls_filtered = vec![];
    
        let sequence_set: HashSet<_> = sequence_calls.iter().collect();
        let mate_sequence_set: HashSet<_> = mate_sequence_calls.iter().collect();
    
        for call in &sequence_calls {
            if mate_sequence_set.contains(call) {
                continue;
            } else {
                sequence_calls_filtered.push(call.clone());
            }
        }
    
        for call in &mate_sequence_calls {
            if !sequence_set.contains(call) {
                mate_sequence_calls_filtered.push(call.clone());
            }
        }
    
        (sequence_calls_filtered, mate_sequence_calls_filtered)
    }

    fn filter_five_prime(sequence_calls: Vec<(String, bool)>, mate_sequence_calls: Vec<(String, bool)>) -> (Vec<String>, Vec<String>) {
        let (seq_filtered_unstranded, mate_filtered_unstranded) = AlignmentOrientation::filter_unstranded(sequence_calls.clone(), mate_sequence_calls.clone());
    
        let mut seq_filtered = vec![];
        let mut mate_filtered = mate_filtered_unstranded.clone();
    
        for call in seq_filtered_unstranded {
            let (feat, rev) = call.clone();
            if rev {
                // Remove reverse orientation calls from sequence_calls and any corresponding call from mate_sequence_calls
                if let Some(pos) = mate_filtered.iter().position(|(mate_call, _)| mate_call == &feat) {
                    mate_filtered.remove(pos);
                }
            } else {
                seq_filtered.push(call);
            }
        }

        // Remove all cases in mate_sequence_calls where there is no corresponding feature in sequence_calls AND it is forward
        mate_filtered.retain(|(mate_feat, rev)| {
            if !rev {
                seq_filtered.iter().any(|(seq_feat, _)| seq_feat == mate_feat)
            } else {
                true
            }
        });
    
        (
            seq_filtered.into_iter().map(|(feat, _)| feat).collect(),
            mate_filtered.into_iter().map(|(feat, _)| feat).collect()
        )
    }
    
    fn filter_three_prime(sequence_calls: Vec<(String, bool)>, mate_sequence_calls: Vec<(String, bool)>) -> (Vec<String>, Vec<String>) {
        let (seq_filtered_unstranded, mate_filtered_unstranded) = AlignmentOrientation::filter_unstranded(sequence_calls.clone(), mate_sequence_calls.clone());
    
        let mut seq_filtered = vec![];
        let mut mate_filtered = mate_filtered_unstranded.clone();
    
        for call in seq_filtered_unstranded {
            let (feat, rev) = call.clone();
            if !rev {
                // Remove forward orientation calls from sequence_calls and any corresponding call from mate_sequence_calls
                if let Some(pos) = mate_filtered.iter().position(|(mate_call, _)| mate_call == &feat) {
                    mate_filtered.remove(pos);
                }
            } else {
                seq_filtered.push(call);
            }
        }

        // Remove all cases in mate_sequence_calls where there is no corresponding feature in sequence_calls AND it is revcomped
        mate_filtered.retain(|(mate_feat, rev)| {
            if *rev {
                seq_filtered.iter().any(|(seq_feat, _)| seq_feat == mate_feat)
            } else {
                true
            }
        });
    
        (
            seq_filtered.into_iter().map(|(feat, _)| feat).collect(),
            mate_filtered.into_iter().map(|(feat, _)| feat).collect()
        )
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
    index: &PseudoAligner,
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)>
) {
    let mut final_filter_reasons: HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)> = HashMap::new();
    let mut filter_reasons: HashMap<String, ((FilterReason, usize), (FilterReason, usize))> = HashMap::new();
    let mut post_triaged_keys: HashMap<String, (FilterReason, AlignmentOrientation)> = HashMap::new();

    // Unpack the index and sequence pairs
    let (sequences, _) = sequence_iterators;
    let (mate_sequences, _) = match mate_sequence_iterators {
        Some((l, r)) => (Some(l), Some(r)),
        None => (None, None),
    };

    // Generate a set of passing scores for the sequences (mate sequences optional), for both the regular and reverse complemented versions of the reference
    let (sequence_scores, matched_sequences) =
        score_sequences(
            sequences,
            mate_sequences,
            sequence_metadata,
            index,
            reference,
            aligner_config,
            &mut filter_reasons
        );

    // Final results hashmap, which will contain all calls that remain after orientation filtering
    let mut results: HashMap<Vec<String>, (i32, Vec<String>, Vec<String>)> = HashMap::new();

    /* Iterate all of the calls generated from aligning the read-pairs to the reference library in the normal orientation.
     * Send those calls to the alignment orientation filtering pipeline.
     * Once this process is complete, the results hashmap should contain one value per unique feature, and a count of the
     * number of read-pairs that returned a call for that feature. */
    for (read_pair_key, call) in sequence_scores.into_iter() {
        AlignmentOrientation::filter_and_coerce_sequence_call_orientations(
            call,
            &mut results,
            reference,
            aligner_config,
            read_pair_key,
            &mut post_triaged_keys
        );
    }

    /*  Merge all the read-filtration information into a single collection, including normal/revcomp library
     *   filter information, and the alignment orientation filter pipeline's returned info */
    for (key, value) in filter_reasons.into_iter() {
        match post_triaged_keys.get(&key) {
            Some(triage) => final_filter_reasons.insert(key, (value.0, value.1, (FilterReason::None, 0), (FilterReason::None, 0), triage.0, triage.1)),
            None => final_filter_reasons.insert(key, (value.0, value.1, (FilterReason::None, 0), (FilterReason::None, 0), FilterReason::None, AlignmentOrientation::None))
        };
    }

    // Iterate the keys and values in the results hashmap, pushing them to the results vector to return
    let mut ret = Vec::new();
    for (key, value) in results.into_iter() {
        ret.push((key, value));
    }

    (ret, matched_sequences, final_filter_reasons)
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

        // Generate score and equivalence class for this read by aligning the sequence against the reference after trimming it for quality.
        let trimmed_read = if !sequence_metadata.is_empty() {
            trim_sequence(&read, sequence_metadata[1].as_str(), &aligner_config)
        } else {
            read.clone()
        };

        let (sequence_alignment, sequence_filter_reason) = if !sequence_metadata.is_empty() && sequence_metadata[37] == "TRUE" {
            (None, Some((FilterReason::SkippedAlignDueToUnpairedDummy, 0.0, 0)))
        } else {
            pseudoalign(&trimmed_read, index, &aligner_config, MIN_READ_LENGTH)
        };

        // If there's a mate sequence, also perform the alignment for it
        let mut mate_sequence_alignment: Option<AlignmentScore> = None;
        let mut mate_sequence_filter_reason: Filter = None;

        if let Some(itr) = &mut mate_sequences {
            let mate_read = itr
                .next()
                .expect("Error -- read and reverse read files do not have matching lengths: ")
                .expect("Error -- could not parse reverse read. Input R2 data malformed.");

            let trimmed_mate_read = if !mate_sequence_metadata.is_empty() {
                trim_sequence(&mate_read, mate_sequence_metadata[1].as_str(), &aligner_config)
            } else {
                read.clone()
            };

            let (score, filter_reason) = if !mate_sequence_metadata.is_empty() && mate_sequence_metadata[37] == "TRUE" {
                (None, Some((FilterReason::SkippedAlignDueToUnpairedDummy, 0.0, 0)))
            } else {
                pseudoalign(&trimmed_mate_read, index, &aligner_config, MIN_READ_LENGTH)
            };

            mate_sequence_alignment = Some(score);
            read_rev = Some(mate_read);
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
        // Alignments are deterministic, so if this key isn't unique, we would've gotten the same alignment for that read anyways
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
                process_equivalence_class_to_feature_list(&sequence_equivalence_class, reference, &aligner_config, false)
            } else if !mate_sequence_equivalence_class.is_empty() {
                process_equivalence_class_to_feature_list(&mate_sequence_equivalence_class, reference, &aligner_config, false)
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
    sequence_call: Vec<String>,
    mate_sequence_call: Vec<String>,
    fallback_on_intersect_fail: bool,
    read_key: String,
    filtered_keys: &mut HashMap<String, (FilterReason, AlignmentOrientation)>
) -> Vec<String> {
    // Get only the calls shared between normal and revcomped alignments
    let class = sequence_call.intersect(mate_sequence_call.clone());

    // If that destroyed the alignment and fallback is on, use get_call_calls instead
    if class.len() == 0 && fallback_on_intersect_fail {
        get_all_calls(sequence_call, mate_sequence_call)
    
    // Otherwise, if the intersection succeeded, return the class with any existing metadata
    } else if class.len() != 0 {
        class
    // If it failed and fallback is off, add this read-pair to the filtered list
    } else {
        filtered_keys.insert(read_key, (FilterReason::ForceIntersectFailure, AlignmentOrientation::None));
        Vec::new()
    }
}

// Permissively append all the calls for a given set of sequence calls, creating one large combined equivalence class
fn get_all_calls(
    mut sequence_call: Vec<String>,
    mut mate_sequence_call: Vec<String>,
) -> Vec<String> {
    // Take the two coerced equivalence classes and merge them
    sequence_call.append(&mut mate_sequence_call);
    sequence_call.unique();
    sequence_call
}

/* Takes an equivalence class and returns a list of strings. If we're processing allele-level data, the strings will be
 * the sequence name of the relevant feature(s). Otherwise, if we're doing a group_by, the equivalence class will be
 * filtered such that there is only one hit per group_by string (e.g. one hit per lineage) and the corresponding strings
 * (e.g. lineage name) will be returned instead of individual sequence names. */
fn process_equivalence_class_to_feature_list(
    equivalence_class: &Vec<u32>,
    reference: &Reference,
    aligner_config: &AlignFilterConfig,
    ignore_group_rollup: bool,
) -> Vec<String> {
    /* If the group_on column is the default, nt_sequence, no feature roll-up occurs. We map over the equivalence class and translate one-to-one from
     * values in that column. */
    let mut results = if ignore_group_rollup || reference.headers[reference.group_on] == "nt_sequence" {
        equivalence_class
            .into_iter()
            .map(|ref_idx| {
                reference.columns[reference.sequence_name_idx][*ref_idx as usize].clone()
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
    if !ignore_group_rollup && aligner_config.discard_multi_hits > 0 && results.len() > aligner_config.discard_multi_hits {
        Vec::new()
    } else {
        // Sort for deterministic feature lists
        results.string_sort_unstable(natural_lexical_cmp);
        results
    }
}

fn unmap(
    feature_list: &Vec<String>,
    reference: &Reference
) -> Vec<u32> {
    feature_list
        .iter()
        .map(|feature| {
            reference.columns[reference.sequence_name_idx]
                .iter()
                .position(|ref_feature| ref_feature == feature)
                .expect("Feature not found in reference columns") as u32
        })
        .collect()
}

fn trim_sequence(sequence: &DnaString, quality: &str, aligner_config: &AlignFilterConfig) -> DnaString {
    // Calculate the optimal trim position and perform trim
    let trimmed_length = maxinfo(quality, aligner_config.trim_target_length, aligner_config.trim_strictness);
    let trimmed_sequence = &sequence.to_string()[..trimmed_length];
    DnaString::from_dna_string(trimmed_sequence)
}

fn maxinfo(quality: &str, target_length: usize, strictness: f64) -> usize {
    const LONGEST_READ: usize = 1000;
    const MAXQUAL: usize = 60;

    // Precompute length scores
    let mut length_scores = vec![0.0; LONGEST_READ];
    for i in 0..LONGEST_READ {
        let pow1 = (target_length as f64 - i as f64 - 1.0).exp();
        let unique = (1.0 / (1.0 + pow1)).ln();
        let coverage = ((i + 1) as f64).ln() * (1.0 - strictness);
        length_scores[i] = unique + coverage;
    }

    // Precompute quality probabilities
    let mut qual_probs = vec![0.0; MAXQUAL + 1];
    for i in 0..=MAXQUAL {
        let prob_correct = 1.0 - 10f64.powf(-((0.5 + i as f64) / 10.0));
        qual_probs[i] = prob_correct.ln() * strictness;
    }

    let norm_ratio = compute_norm_ratio(&length_scores, LONGEST_READ * 2 as usize)
        .max(compute_norm_ratio(&qual_probs, LONGEST_READ * 2 as usize));

    let length_scores = normalize(&length_scores, norm_ratio);
    let qual_probs = normalize(&qual_probs, norm_ratio);

    // Calculate the optimal trim position
    let mut accum_quality = 0;
    let mut max_score = f64::MIN;
    let mut max_score_position = 0;

    for (i, &q_char) in quality.as_bytes().iter().enumerate() {
        let q = q_char as usize;
        let q = if q > MAXQUAL { MAXQUAL } else { q };
        
        accum_quality += qual_probs[q];
        let ls = length_scores.get(i).copied().unwrap_or(0);
        let score = ls + accum_quality;

        if score as f64 >= max_score {
            max_score = score as f64;
            max_score_position = i + 1;
        }
    }

    if max_score_position < 1 || max_score == 0.0 {
        return 0;
    } else if max_score_position < quality.len() {
        return max_score_position
    } else {
        return quality.len()
    }
}

fn compute_norm_ratio(array: &[f64], margin: usize) -> f64 {
    let mut max_val = array[0].abs();

    for &val in &array[1..] {
        let abs_val = val.abs();
        if abs_val > max_val {
            max_val = abs_val;
        }
    }

    (std::i64::MAX as f64) / (max_val * margin as f64)
}

fn normalize(array: &[f64], ratio: f64) -> Vec<i64> {
    array.iter().map(|&val| (val * ratio) as i64).collect()
}

// Align a sequence against the reference library, with settings specified by the aligner configuration
fn pseudoalign(
    sequence: &DnaString,
    reference_index: &PseudoAligner,
    config: &AlignFilterConfig,
    min_read_length: usize
) -> (
    AlignmentScore,
    Filter
) {
    // Ensure all reads are above a minimum length
    if sequence.len() < min_read_length {
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
                config.num_mismatches,
                mismatches
            )
        }
        None => (None, Some((FilterReason::NoMatch, 0.0, 0))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use debruijn::kmer::Kmer30;
    use reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR;

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

    fn create_test_sequence() -> DnaString {
        DnaString::from_dna_string("ACGTACGTACGTACGTACGT")
    }
    
    fn create_test_quality() -> String {
        "IIIIIIIIIIIIIIIIIIII".to_string()
    }

    fn adjust_quality(quality: &str) -> String {
        quality.chars()
            .map(|c| (c as u8 - 33) as char)
            .collect()
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
            trim_strictness: 0.5,
            trim_target_length: 15
        }
    }

    #[test]
    fn test_short_read() {
        let sequence = DnaString::from_dna_string("ACG");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config, 12);
        assert_eq!(result.1, Some((FilterReason::ShortRead, 0.0, 0)));
    }

    #[test]
    fn test_high_entropy_read() {
        let sequence = DnaString::from_dna_string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config, 12);
        assert_eq!(result.1, Some((FilterReason::HighEntropy, 0.0, 0)));
    }

    #[test]
    fn test_no_alignment_match() {
        let sequence = DnaString::from_dna_string("CCTGAGATTTCGAGCTCGTAACGTGACCTACGGACAC");
        let config = setup_config();
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config, 12);
        assert_eq!(result.1, Some((FilterReason::NoMatch, 0.0, 0)));
    }

    #[test]
    fn test_valid_alignment() {
        let sequence = DnaString::from_dna_string("TGCATGCATGCATGCATGCATGCATGCATGCA");
        let mut config = setup_config();
        let reference_index = setup_pseudoaligner();
        config.score_threshold = 32;
        let result = pseudoalign(&sequence, &reference_index, &config, 12);
        assert_eq!(result.0, Some((vec![1], 1.0, 32)));
        assert_eq!(result.1, None);
    }

    #[test]
    fn test_score_and_match_threshold_filtering() {
        let sequence = DnaString::from_dna_string("TGCATGCATGCATGCATGCATGCATGCATGCA");
        let mut config = setup_config();
        config.score_threshold = 1000;
        let reference_index = setup_pseudoaligner();
        let result = pseudoalign(&sequence, &reference_index, &config, 12);
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
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, false),
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
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, false),
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
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, false),
            vec!["geneA", "seq2"]
        );
    }

    #[test]
    fn test_ignore_groupby() {
        let mut ref_data = setup_reference();
        ref_data.columns[1] = vec!["geneA".to_string(), "".to_string(), "geneA".to_string()];
        ref_data.group_on = 1;
        let config = setup_config();
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, true),
            vec!["seq1", "seq2", "seq3"]
        );
    }

    #[test]
    fn test_discard_multi_hits() {
        let ref_data = setup_reference();
        let mut config = setup_config();
        config.discard_multi_hits = 1;
        let equivalence_class = vec![0, 1, 2];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, false),
            Vec::<String>::new()
        );
    }

    #[test]
    fn test_empty_equivalence_class() {
        let ref_data = setup_reference();
        let config = setup_config();
        let equivalence_class = vec![];
        assert_eq!(
            process_equivalence_class_to_feature_list(&equivalence_class, &ref_data, &config, false),
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

        let result_1 = process_equivalence_class_to_feature_list(&equivalence_class_1, &ref_data, &config, false);
        let result_2 = process_equivalence_class_to_feature_list(&equivalence_class_2, &ref_data, &config, false);

        assert_eq!(result_1, result_2);
        assert_eq!(result_1, vec!["geneA", "geneB"]);
    }

    #[test]
    fn test_parse_calls() {
        let calls = vec![
            "feat1".to_string(),
            format!("feat2{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            "feat3".to_string(),
            format!("feat4{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            format!("feat4{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            "feat4".to_string(),
        ];
        let expected = vec![
            ("feat1".to_string(), false),
            ("feat2".to_string(), true),
            ("feat3".to_string(), false),
            ("feat4".to_string(), true),
            ("feat4".to_string(), true),
            ("feat4".to_string(), false),
        ];
        assert_eq!(AlignmentOrientation::parse_calls(calls), expected);
    }

    #[test]
    fn test_filter_unstranded() {
        let sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat2".to_string(), true),
            ("feat4".to_string(), true),
            ("feat5".to_string(), true),
        ];
        let mate_sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat3".to_string(), false),
            ("feat4".to_string(), false),
            ("feat5".to_string(), true),
        ];
        let expected_sequence = vec![
            ("feat2".to_string(), true),
            ("feat4".to_string(), true),
        ];
        let expected_mate = vec![
            ("feat3".to_string(), false),
            ("feat4".to_string(), false),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_unstranded(sequence_calls, mate_sequence_calls);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_five_prime() {
        let sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat2".to_string(), true),
            ("feat4".to_string(), false),
            ("feat5".to_string(), true),
            ("feat6".to_string(), false),
        ];
        let mate_sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat3".to_string(), true),
            ("feat4".to_string(), true),
            ("feat5".to_string(), false),
            ("feat7".to_string(), false),
        ];
        let expected_sequence = vec![
            "feat4".to_string(),
            "feat6".to_string(),
        ];
        let expected_mate = vec![
            "feat3".to_string(),
            "feat4".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_five_prime(sequence_calls, mate_sequence_calls);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_three_prime() {
        let sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat2".to_string(), true),
            ("feat4".to_string(), false),
            ("feat5".to_string(), true),
            ("feat6".to_string(), false),
        ];
        let mate_sequence_calls = vec![
            ("feat1".to_string(), false),
            ("feat3".to_string(), false),
            ("feat4".to_string(), true),
            ("feat5".to_string(), false),
            ("feat7".to_string(), true),
        ];
        let expected_sequence = vec![
            "feat2".to_string(),
            "feat5".to_string(),
        ];
        let expected_mate = vec![
            "feat3".to_string(),
            "feat5".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_three_prime(sequence_calls, mate_sequence_calls);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_orientation_on_library_chemistry_none() {
        let sequence_calls = vec![
            "feat1".to_string(),
            format!("feat2{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let mate_sequence_calls = vec![
            "feat3".to_string(),
            format!("feat4{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let lib_type = LibraryChemistry::None;
        let expected_sequence = vec![
            "feat1".to_string(),
            "feat2".to_string(),
        ];
        let expected_mate = vec![
            "feat3".to_string(),
            "feat4".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_orientation_on_library_chemistry(sequence_calls, mate_sequence_calls, &lib_type);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_orientation_on_library_chemistry_unstranded() {
        let sequence_calls = vec![
            "feat1".to_string(),
            "feat2".to_string(),
            format!("feat4{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            "feat5".to_string(),
        ];
        let mate_sequence_calls = vec![
            "feat1".to_string(),
            "feat3".to_string(),
            "feat4".to_string(),
            format!("feat5{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let lib_type = LibraryChemistry::Unstranded;
        let expected_sequence = vec![
            "feat2".to_string(),
            "feat4".to_string(),
            "feat5".to_string(),
        ];
        let expected_mate = vec![
            "feat3".to_string(),
            "feat4".to_string(),
            "feat5".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_orientation_on_library_chemistry(sequence_calls, mate_sequence_calls, &lib_type);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_orientation_on_library_chemistry_five_prime() {
        let sequence_calls = vec![
            "feat1".to_string(),
            format!("feat2{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            "feat3".to_string(),
            "feat5".to_string(),
            "feat6".to_string(),
            format!("feat8{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let mate_sequence_calls = vec![
            "feat1".to_string(),
            "feat3".to_string(),
            "feat8".to_string(),
            "feat4".to_string(),
            format!("feat5{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            format!("feat7{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let lib_type = LibraryChemistry::FivePrime;
        let expected_sequence = vec![
            "feat5".to_string(),
            "feat6".to_string(),
        ];
        let expected_mate = vec![
            "feat5".to_string(),
            "feat7".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_orientation_on_library_chemistry(sequence_calls, mate_sequence_calls, &lib_type);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_filter_orientation_on_library_chemistry_three_prime() {
        let sequence_calls = vec![
            "feat1".to_string(),
            format!("feat2{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            "feat3".to_string(),
            format!("feat5{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let mate_sequence_calls = vec![
            "feat7".to_string(),
            "feat1".to_string(),
            "feat5".to_string(),
            format!("feat6{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
            format!("feat4{}rev", SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR),
        ];
        let lib_type = LibraryChemistry::ThreePrime;
        let expected_sequence = vec![
            "feat2".to_string(),
            "feat5".to_string(),
        ];
        let expected_mate = vec![
            "feat7".to_string(),
            "feat5".to_string(),
        ];
        let (filtered_sequence, filtered_mate) = AlignmentOrientation::filter_orientation_on_library_chemistry(sequence_calls, mate_sequence_calls, &lib_type);
        assert_eq!(filtered_sequence, expected_sequence);
        assert_eq!(filtered_mate, expected_mate);
    }

    #[test]
    fn test_no_duplicates() {
        let calls = vec![
            "name1".to_string(),
            "name2".to_string(),
            "name3".to_string(),
            "name4".to_string()
        ];
        let expected = vec![
            "name1".to_string(),
            "name2".to_string(),
            "name3".to_string(),
            "name4".to_string()
        ];
        assert_eq!(AlignmentOrientation::filter_read_calls_with_orientation(calls), expected);
    }

    #[test]
    fn test_with_duplicates() {
        let calls = vec![
            "name1".to_string(),
            "name1".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name2".to_string(),
            "name3".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name3".to_string(),
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
        ];
        let expected = vec![
            "name2".to_string(),
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev"
        ];
        assert_eq!(AlignmentOrientation::filter_read_calls_with_orientation(calls), expected);
    }

    #[test]
    fn test_all_revs() {
        let calls = vec![
            "name1".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name2".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name3".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev"
        ];
        let expected = vec![
            "name1".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name2".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name3".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev"
        ];
        assert_eq!(AlignmentOrientation::filter_read_calls_with_orientation(calls), expected);
    }

    #[test]
    fn test_mixed() {
        let calls = vec![
            "name1".to_string(),
            "name2".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name1".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name3".to_string(),
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name3".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name5".to_string(),
            "name6".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name7".to_string(),
            "name8".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name9".to_string(),
            "name8".to_string()
        ];
        let expected = vec![
            "name2".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name4".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name5".to_string(),
            "name6".to_string() + reference_library::SPECIAL_REVCOMP_FEATURE_NAME_SEPARATOR + "rev",
            "name7".to_string(),
            "name9".to_string()
        ];
        assert_eq!(AlignmentOrientation::filter_read_calls_with_orientation(calls), expected);
    }

    #[test]
    fn test_unmap() {
        let reference = Reference {
            headers: vec!["nt_sequence".to_string()],
            group_on: 0,
            columns: vec![
                vec![
                    "feature1".to_string(),
                    "feature2".to_string(),
                    "feature3".to_string(),
                ],
            ],
            sequence_name_idx: 0,
            sequence_idx: 1,
        };

        let feature_list = vec!["feature1".to_string(), "feature2".to_string(), "feature3".to_string()];

        let result = unmap(&feature_list, &reference);

        assert_eq!(result, vec![0, 1, 2]);
    }

    #[test]
    fn test_unmap_unorder() {
        let reference = Reference {
            headers: vec!["nt_sequence".to_string()],
            group_on: 0,
            columns: vec![
                vec![
                    "feature1".to_string(),
                    "feature2".to_string(),
                    "feature3".to_string(),
                ],
            ],
            sequence_name_idx: 0,
            sequence_idx: 1,
        };

        let feature_list = vec!["feature2".to_string(), "feature1".to_string(), "feature3".to_string()];

        let result = unmap(&feature_list, &reference);

        assert_eq!(result, vec![1, 0, 2]);
    }

    #[test]
    fn test_process_and_unmap() {
        let reference = Reference {
            headers: vec!["nt_sequence".to_string()],
            group_on: 0,
            columns: vec![
                vec![
                    "feature1".to_string(),
                    "feature2".to_string(),
                    "feature3".to_string(),
                ],
            ],
            sequence_name_idx: 0,
            sequence_idx: 1,
        };

        let equivalence_class = vec![0, 1, 2];
        let aligner_config = setup_config();
        let ignore_group_rollup = true;

        let feature_list = process_equivalence_class_to_feature_list(
            &equivalence_class,
            &reference,
            &aligner_config,
            ignore_group_rollup,
        );

        let result = unmap(&feature_list, &reference);

        assert_eq!(result, equivalence_class);
    }

    #[test]
    fn test_get_all_calls_both_present() {
        let sequence_call = vec!["1".to_string(), "2".to_string(), "3".to_string()];
        let mate_sequence_call = vec!["4".to_string(), "5".to_string(), "6".to_string()];
        let result = get_all_calls(sequence_call, mate_sequence_call);
        assert_eq!(result, vec!["1".to_string(), "2".to_string(), "3".to_string(), "4".to_string(), "5".to_string(), "6".to_string()]);
    }

    #[test]
    fn test_get_all_calls_empty_metadata() {
        let sequence_call = vec!["1".to_string(), "2".to_string(), "3".to_string()];
        let mate_sequence_call = vec!["4".to_string(), "5".to_string(), "6".to_string()];
        let result = get_all_calls(sequence_call, mate_sequence_call);
        assert_eq!(result, vec!["1".to_string(), "2".to_string(), "3".to_string(), "4".to_string(), "5".to_string(), "6".to_string()]);
    }

    #[test]
    fn test_get_intersecting_reads_intersect_success() {
        let sequence_call = vec!["1".to_string(), "2".to_string(), "3".to_string(), "4".to_string()];
        let mate_sequence_call = vec!["4".to_string(), "5".to_string(), "6".to_string()];
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(sequence_call, mate_sequence_call, false, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, vec!["4".to_string()]);
        assert!(filtered_keys.is_empty());
    }

    #[test]
    fn test_get_intersecting_reads_fallback_on_intersect_fail() {
        let sequence_call = vec!["1".to_string(), "2".to_string(), "3".to_string()];
        let mate_sequence_call = vec!["4".to_string(), "5".to_string(), "6".to_string()];
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(sequence_call, mate_sequence_call, true, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, vec!["1".to_string(), "2".to_string(), "3".to_string(), "4".to_string(), "5".to_string(), "6".to_string()]);
        assert!(filtered_keys.is_empty());
    }

    #[test]
    fn test_get_intersecting_reads_intersect_fail_no_fallback() {
        let sequence_call = vec!["1".to_string(), "2".to_string(), "3".to_string()];
        let mate_sequence_call = vec!["4".to_string(), "5".to_string(), "6".to_string()];
        let mut filtered_keys = HashMap::new();
        let result = get_intersecting_reads(sequence_call, mate_sequence_call, false, "read_key".to_string(), &mut filtered_keys);
        assert_eq!(result, Vec::<String>::new());
        assert_eq!(filtered_keys.get("read_key"), Some(&(FilterReason::ForceIntersectFailure, AlignmentOrientation::None)));
    }

    #[test]
    fn test_trim_sequence_high_quality() {
        let sequence = create_test_sequence();
        let quality = create_test_quality();
        let config = setup_config();
        
        let adjusted_quality = adjust_quality(&quality);
        let trimmed_sequence = trim_sequence(&sequence, &adjusted_quality, &config);
        
        assert_eq!(trimmed_sequence.to_string(), "ACGTACGTACGTACGTACGT");
    }
    
    #[test]
    fn test_trim_sequence_low_quality() {
        let sequence = create_test_sequence();
        let quality = "!!!!!!!!!!!!!!!!!!!!".to_string(); // Low quality scores
        let mut config = setup_config();
        config.trim_strictness = 0.9;
        
        let adjusted_quality = adjust_quality(&quality);
        let trimmed_sequence = trim_sequence(&sequence, &adjusted_quality, &config);
        
        assert_eq!(trimmed_sequence.to_string(), "A");
    }
    
    #[test]
    fn test_trim_sequence_mixed_quality() {
        let sequence = create_test_sequence();
        let quality = "IIIIII!!!!!!IIIIII".to_string(); // Mixed quality scores
        let mut config = setup_config();
        config.trim_strictness = 0.8;
        
        let adjusted_quality = adjust_quality(&quality);
        let trimmed_sequence = trim_sequence(&sequence, &adjusted_quality, &config);
        
        assert_eq!(trimmed_sequence.to_string(), "ACGTAC");
    }
    
    #[test]
    fn test_maxinfo_all_high_quality() {
        let quality = "IIIIIIIIIIIIIIIIIIII".to_string();
        let target_length = 15;
        let strictness = 0.5;
        
        let adjusted_quality = adjust_quality(&quality);
        let trim_length = maxinfo(&adjusted_quality, target_length, strictness);
        
        assert_eq!(trim_length, 20);
    }
    
    #[test]
    fn test_maxinfo_all_low_quality() {
        let quality = "!!!!!!!!!!!!!!!!!!!!".to_string();
        let target_length = 15;
        let strictness = 0.9;
        
        let adjusted_quality = adjust_quality(&quality);
        let trim_length = maxinfo(&adjusted_quality, target_length, strictness);
        
        assert_eq!(trim_length, 1);
    }
    
    #[test]
    fn test_maxinfo_mixed_quality() {
        let quality = "IIIIII!!!!!!IIIIII".to_string();
        let target_length = 15;
        let strictness = 0.7;
        
        let adjusted_quality = adjust_quality(&quality);
        let trim_length = maxinfo(&adjusted_quality, target_length, strictness);
        
        assert_eq!(trim_length, 6);
    }
    
    #[test]
    fn test_maxinfo_strictness_1() {
        let quality = "IIIIIIIIIIIIIIIIIIII".to_string();
        let target_length = 15;
        let strictness = 1.0;
        
        let adjusted_quality = adjust_quality(&quality);
        let trim_length = maxinfo(&adjusted_quality, target_length, strictness);
        
        assert_eq!(trim_length, 20);
    }
    
    #[test]
    fn test_maxinfo_strictness_0() {
        let quality = "IIIIIIIIIIIIIIIIIIII".to_string();
        let target_length = 15;
        let strictness = 0.0;
        
        let adjusted_quality = adjust_quality(&quality);
        let trim_length = maxinfo(&adjusted_quality, target_length, strictness);
        
        assert_eq!(trim_length, 20);
    }
}