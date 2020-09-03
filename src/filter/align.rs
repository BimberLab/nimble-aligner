/* Takes a given alignment and returns it if it's above the given match threshold and
 * has an equivalence class. Can be configured to discard alignments with more than one match. */
pub fn filter_by_alignment_score(score: usize, equiv_class: Vec<u32>, score_threshold: usize, discard_multiple_matches: bool) -> Option<(Vec<u32>, usize)> {
  if score >= score_threshold && !equiv_class.is_empty() { 
    if discard_multiple_matches && equiv_class.len() > 1 {
      None
    } else {
      Some((equiv_class, score)) 
    }
  } else { 
    None 
  }
}

#[cfg(test)]
mod tests {
    // Tests for filter_by_alignment_score
    // Case where the score is higher than the threshold -- should not filter the alignment
    #[test]
    fn do_not_filter() {
      let score = 100;
      let equiv_class = vec![1, 2];

      let results = super::filter_by_alignment_score(score, equiv_class, 50, false);
      let expected_results = Some((vec![1, 2], 100));

      assert_eq!(results, expected_results);
    }

    // Case where the score is lower than the threshold -- should filter the alignment
    #[test]
    fn filter() {
      let score = 25;
      let equiv_class = vec![1, 2];

      let results = super::filter_by_alignment_score(score, equiv_class, 50, false);
      let expected_results = None;

      assert_eq!(results, expected_results);
    }

    /* Case where the score is higher than the threshold, but matches more than one
     * reference and has discard_multiple_matches = true -- should filter the alignment */
    #[test]
    fn filter_multiple_matches() {
      let score = 100;
      let equiv_class = vec![1, 2];

      let results = super::filter_by_alignment_score(score, equiv_class, 50, true);
      let expected_results = None;

      assert_eq!(results, expected_results);
    }
}