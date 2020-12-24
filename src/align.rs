use crate::filter;
use crate::reference_library;

use std::io::{Error};
use std::collections::HashMap;

use debruijn::dna_string::DnaString;
use array_tool::vec::Intersect;
use reference_library::ReferenceMetadata;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

pub enum IntersectLevel {
  NoIntersect,
  IntersectWithFallback,
  ForceIntersect
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
  pub discard_multi_hits: usize
}

/* Takes a set of sequences and optionally, reverse sequences, a debrujin map index of the reference
 * genome, the reference library metadata object, and the aligner configuration, and performs a
 * debrujin-graph based pseduoalignment, returning a score for each readable reference in the reference
 * genome.
 * This function does some alignment-time filtration based on the provided configuration. */
pub fn score<I>(sequences: I, mut reverse_sequences: Option<I>, index: PseudoAligner, reference_metadata: &ReferenceMetadata,
  config: &AlignFilterConfig) -> Vec<(Vec<String>, i32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>
  {

  // HashMap of the alignment results. The keys are either strong hits or equivalence classes of hits
  let mut score_map: HashMap<Vec<String>, i32> = HashMap::new();

  // Iterate over every read/reverse read pair and align it, incrementing scores for the matching references/equivalence classes
  for read in sequences {

    let read = read.expect("Error -- could not parse read. Input R1 data malformed.");
    /* Generate score and equivalence class for this read by aligning the sequence against
     * the current reference, if there is a match.*/
    let seq_score = pseudoalign(&read, &index, &config);

    // If there's a reversed sequence, do the paired-end alignment
    let mut rev_seq_score = None;
    if let Some(itr) = &mut reverse_sequences {
      let reverse_read = itr.next().expect(
        "Error -- read and reverse read files do not have matching lengths: "
      ).expect("Error -- could not parse reverse read. Input R2 data malformed.");
      rev_seq_score = Some(pseudoalign(&reverse_read, &index, &config));
    }

    // If there are no reverse sequences, ignore the require_valid_pair filter
    if reverse_sequences.is_some() {

      // Otherwise, check if the user enabled the require_valid_pair filter
      if config.require_valid_pair {
        
        // Unpack the data to check if the pairs have the same eq_class. If they both contain data, do the comparison
        if let (Some(Some((mut rev_eq_class, _))), Some((mut eq_class, _))) = ((&rev_seq_score).clone(), (&seq_score).clone()) {

          // Sort the vectors and compare them by counting matching elements. If they don't match, don't modify the score for this read
          rev_eq_class.sort();
          eq_class.sort();
          let matching = eq_class.iter().zip(&rev_eq_class).filter(|&(classl, classr)| classl == classr).count();

          if matching != eq_class.len() || matching != rev_eq_class.len() {
            continue;
          }
        } else { // Otherwise, require_valid_pair is on and rev_score != seq_score, or they're both None. In either case, don't modify score
          continue
        }
      }
    }

    // Take the "best" alignment. The specific behavior is determined by the intersect level set in the aligner config
    let match_eqv_class = match config.intersect_level {
      IntersectLevel::NoIntersect => get_best_reads(&seq_score, &rev_seq_score),
      IntersectLevel::IntersectWithFallback => get_intersecting_reads(&seq_score, &rev_seq_score, true),
      IntersectLevel::ForceIntersect => get_intersecting_reads(&seq_score, &rev_seq_score, false)
    };

    if !match_eqv_class.is_empty() {
      let key = get_score_map_key(match_eqv_class, reference_metadata, &config); // Process the equivalence class into a score key

      // Add the key to the score map and increment the score
      let accessor = score_map.entry(key).or_insert(0); 
      *accessor += 1;
    }
  }

  // Update the results map
  let mut results = Vec::new();
  for (key, value) in score_map.into_iter() {
    results.push((key, value));
  }

  results
}


// Return matches that match in both seq_score and rev_seq_score; if soft intersection is enabled, fall back to best read if one of the reads is empty
fn get_intersecting_reads(seq_score: &Option<(Vec<u32>, usize)>, rev_seq_score: &Option<Option<(Vec<u32>, usize)>>, fallback_on_intersect_fail: bool) -> Vec<u32> {
  if let (Some((eqv_class_seq, _)), Some(Some((eqv_class_rev_seq, _)))) = (&seq_score, &rev_seq_score) {
    eqv_class_seq.intersect(eqv_class_rev_seq.to_vec())
  } else if fallback_on_intersect_fail {
      get_best_reads(seq_score, rev_seq_score)
  } else {
    Vec::new()
  }
}


// Return matches from seq_score -- otherwise, return matches from rev_seq_score
fn get_best_reads(seq_score: &Option<(Vec<u32>, usize)>, rev_seq_score: &Option<Option<(Vec<u32>, usize)>>) -> Vec<u32> {
  if let Some((eqv_class, _)) = &seq_score {
      (*eqv_class).clone()
  } else if let Some(Some((eqv_class, _))) = &rev_seq_score {
    (*eqv_class).clone()
  } else {
    Vec::new()
  }
}


/* Takes a equivalence class and returns a list of strings. If we're processing allele-level data, the strings will be
 * the nt_sequences of the relevant alleles. Otherwise, if we're doing a group_by, the equivalence class will be
 * filtered such that there is only one hit per group_by string (e.g. one hit per lineage) and the corresponding strings
 * (e.g. lineage name) will be returned. */
fn get_score_map_key(equiv_class: Vec<u32>, reference_metadata: &ReferenceMetadata, config: &AlignFilterConfig) -> Vec<String> {
  if reference_metadata.headers[reference_metadata.group_on] == "nt_sequence" {
    equiv_class.into_iter().map(|ref_idx| reference_metadata.columns[reference_metadata.group_on][ref_idx as usize].clone()).collect()
  } else {
    let mut results = Vec::new();

    for ref_idx in equiv_class {
      let group = &reference_metadata.columns[reference_metadata.group_on][ref_idx as usize];
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
fn pseudoalign(sequence: &DnaString, reference_index: &PseudoAligner, config: &AlignFilterConfig) -> Option<(Vec<u32>, usize)> {
  // Perform alignment
  match reference_index.map_read_with_mismatch(sequence, config.num_mismatches) {
    Some((equiv_class, score, mismatches)) => {
      
      // Filter nonzero mismatch
      if config.discard_nonzero_mismatch && mismatches != 0 {
        return None
      }

      // Filter by score and match threshold
      filter::align::filter_alignment_by_metrics(score, equiv_class, config.score_threshold, config.discard_multiple_matches)
    },
    None => None
  }
}
