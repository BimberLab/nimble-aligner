use crate::filter;
use crate::reference_library;

use std::io::Error;
use std::collections::HashMap;

use debruijn::dna_string::DnaString;
use reference_library::ReferenceMetadata;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

#[derive(Debug)]
pub struct AlignFilterConfig {
  pub reference_genome_size: usize,
  pub score_threshold: usize,
  pub num_mismatches: usize,
  pub discard_differing_read_pairs: bool,   // TODO
  pub discard_nonzero_mismatch: bool,
  pub discard_multiple_matches: bool,
  pub score_filter: i32
}

/* Takes a set of sequences and optionally, reverse sequences, a debrujin map index of the reference
 * genome, the size of the reference genome, and the threshold to match against, and performs a
 * debrujin-graph based pseduoalignment, returning a score for each readable reference in the reference
 * genome.
 * This function does some alignment-time filtration based on the provided parameters. */
pub fn score<I>(sequences: I, mut reverse_sequences: Option<I>, index: PseudoAligner, reference_metadata: &ReferenceMetadata,
  config: &AlignFilterConfig) -> Vec<(String, i32)>
  where 
    I: Iterator<Item = Result<DnaString, Error>>
  {

  let mut score_map: HashMap<String, (i32, bool)> = HashMap::new();

  for read in sequences {
    for (_, value) in score_map.iter_mut() {
      *value = (value.0, false);
    }

    /* Generate score and equivalence class for this read by aligning the sequence against
     * the current reference. This alignment returns any scores that are greater than the match threshold. */
    let seq_score = pseudoalign(read, &index, config.num_mismatches, config.score_threshold,
      config.discard_multiple_matches, config.discard_nonzero_mismatch);
    let mut rev_seq_score = None;

    // If there's a reversed sequence, do the paired-end alignment
    if let Some(itr) = &mut reverse_sequences {
      let reverse_read = itr.next().expect("Error -- read and reverse read files do not have matching lengths: ");
      rev_seq_score = Some(pseudoalign(reverse_read, &index, config.num_mismatches, config.score_threshold, 
        config.discard_multiple_matches, config.discard_nonzero_mismatch));
    }

    // Get the score and the associated equivalence class of the forward sequence
    let mut max_score = 0;
    let mut max_eqv_class = Vec::new();
    if let Some((eqv_class, score)) = seq_score {
      max_score = score;
      max_eqv_class = eqv_class;
    } 

    // If there's a reverse sequence and it matches, compare to the forward score and replace it if this score is higher
    if let Some(Some((eqv_class, score))) = rev_seq_score {
      if score > max_score {
        max_eqv_class = eqv_class;
      }
    }

    // If there was a match, update the results accordingly
    if !max_eqv_class.is_empty() {
      for idx in max_eqv_class {
        let key = &reference_metadata.columns[reference_metadata.group_on][idx as usize];
        let accessor = score_map.entry(key.to_string()).or_insert((1, true));

        if accessor.1 == false {
          accessor.0 += 1;
          accessor.1 = true;
        }
      }
    }
  }

  let mut results = Vec::new();
  for (key, value) in score_map.into_iter() {
    results.push((key, value.0));
  }

  results
}


// Align the given sequence against the given reference with a score threshold
fn pseudoalign(sequence: Result<DnaString, Error>, reference_index: &PseudoAligner,
  match_threshold: usize, allowed_mismatches: usize, discard_multiple_matches: bool, discard_nonzero_mismatch: bool) -> Option<(Vec<u32>, usize)> {
  if sequence.is_err() {
    return None;
  }

  // Perform alignment
  match reference_index.map_read_with_mismatch(&sequence.unwrap(), allowed_mismatches) {
    Some((equiv_class, score, mismatches)) => {
      
      // Filter nonzero mismatch
      if discard_nonzero_mismatch && mismatches != 0 {
        return None
      }

      // Filter by score and match threshold
      filter::align::filter_alignment_by_metrics(score, equiv_class, match_threshold, discard_multiple_matches)
    },
    None => None
  }
}