use crate::filter;

use std::io::Error;
use debruijn::dna_string::DnaString;

type PseudoAligner = debruijn_mapping::pseudoaligner::Pseudoaligner<debruijn::kmer::VarIntKmer<u64, debruijn::kmer::K20>>;

pub struct AlignFilterConfig {
  pub reference_genome_size: usize,
  pub score_threshold: usize,
  pub num_mismatches: usize,
  pub discard_differing_read_pairs: bool,   // TODO
  pub discard_nonzero_mismatch: bool,
  pub discard_multiple_matches: bool
}

/* Takes a set of sequences and optionally, reverse sequences, a debrujin map index of the reference
 * genome, the size of the reference genome, and the threshold to match against, and performs a
 * debrujin-graph based pseduoalignment, returning a score for each readable reference in the reference
 * genome.
 * This function does some alignment-time filtration based on the provided parameters. */
pub fn score(sequences: impl Iterator<Item = Result<DnaString, Error>>, 
  mut reverse_sequences: Option<impl Iterator<Item = Result<DnaString, Error>>>,
  index: PseudoAligner, config: AlignFilterConfig) -> Vec<i32> {
  sequences.fold(vec![0; config.reference_genome_size],
    |mut acc, read| {
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
          acc[idx as usize] += 1;
        }
      }

      acc
    }
  )
}


// Align the given sequence against the given reference with a score threshold
fn pseudoalign(sequence: Result<DnaString, Error>, reference_index: &PseudoAligner,
  match_threshold: usize, allowed_mismatches: usize, discard_multiple_matches: bool, discard_nonzero_mismatch: bool) -> Option<(Vec<u32>, usize)> {
  if sequence.is_err() {
    return None;
  }

  match reference_index.map_read_with_mismatch(&sequence.unwrap(), allowed_mismatches) {
    Some((equiv_class, score, mismatches)) => {
      if discard_nonzero_mismatch && mismatches != 0 {
        return None
      }

      filter::align::filter_by_alignment_score(score, equiv_class, match_threshold, discard_multiple_matches)
    },
    None => None
  }
}