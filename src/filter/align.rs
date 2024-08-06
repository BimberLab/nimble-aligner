use crate::align::FilterReason;

// Takes the results of an alignment and returns it if it passes various thresholds, otherwise filtering it.
pub fn filter_alignment_by_metrics(
    equivalence_class: Vec<u32>,
    score: usize,
    normalized_score: f64,
    score_threshold: usize,
    normalized_score_threshold: f64,
    discard_multiple_matches: bool,
    mismatch_threshold: usize,
    mismatches: usize
) -> (
    Option<(Vec<u32>, f64, usize)>,
    Option<(FilterReason, f64, usize)>,
) {
    if score >= score_threshold && normalized_score >= normalized_score_threshold && !equivalence_class.is_empty() {
        if discard_multiple_matches && equivalence_class.len() > 1 {
            (
                None,
                Some((
                    FilterReason::DiscardedMultipleMatch,
                    normalized_score,
                    score,
                )),
            )
        } else if mismatches > mismatch_threshold {
            (
                None,
                Some((
                    FilterReason::AboveMismatchThreshold,
                    normalized_score,
                    score,
                )),
            )
        } else {
            (Some((equivalence_class, normalized_score, score)), None)
        }
    } else {
        (
            None,
            Some((FilterReason::ScoreBelowThreshold, normalized_score, score)),
        )
    }
}

#[cfg(test)]
mod tests {
    // Tests for filter_by_alignment_score
    // Case where the score is higher than the threshold -- should not filter the alignment
    #[test]
    fn do_not_filter() {
        let score = 50;
        let normalized_score = 1.0;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (results, _) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            false,
            0,
            0
        );
        let expected_results = Some((vec![1, 2], 1.0, 50));

        assert_eq!(results, expected_results);
    }

    // Case where the score is lower than the threshold -- should filter the alignment
    #[test]
    fn filter() {
        let score = 10;
        let normalized_score = 0.10;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (_, results) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            false,
            0,
            0
        );
        let expected_results = Some((super::FilterReason::ScoreBelowThreshold, 0.10, 10));

        assert_eq!(results, expected_results);
    }

    /* Case where the score is higher than the threshold, but matches more than one
     * reference and has discard_multiple_matches = true -- should filter the alignment */
    #[test]
    fn filter_multiple_matches() {
        let score = 50;
        let normalized_score = 1.0;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (_, results) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            true,
            0,
            0
        );

        let expected_results = Some((super::FilterReason::DiscardedMultipleMatch, 1.0, 50));

        assert_eq!(results, expected_results);
    }

    // Case where the mismatches are lower than the threshold -- should not filter the alignment
    #[test]
    fn do_not_filter_mismatches() {
        let score = 50;
        let normalized_score = 1.0;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (results, _) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            false,
            1,
            0
        );
        let expected_results = Some((vec![1, 2], 1.0, 50));

        assert_eq!(results, expected_results);
    }

    // Case where the mismatches are equal to the threshold -- should not filter the alignment
    #[test]
    fn do_not_filter_mismatches_equal() {
        let score = 50;
        let normalized_score = 1.0;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (results, _) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            false,
            1,
            1
        );
        let expected_results = Some((vec![1, 2], 1.0, 50));

        assert_eq!(results, expected_results);
    }

    // Case where the mismatches are equal to the threshold -- should not filter the alignment
    #[test]
    fn filter_mismatches() {
        let score = 50;
        let normalized_score = 1.0;
        let score_threshold = 20;
        let score_percent = 0.5;
        let equiv_class = vec![1, 2];

        let (_, results) = super::filter_alignment_by_metrics(
            equiv_class,
            score,
            normalized_score,
            score_threshold,
            score_percent,
            false,
            1,
            2
        );
        let expected_results = Some((super::FilterReason::AboveMismatchThreshold, 1.0, 50));

        assert_eq!(results, expected_results);
    }
}
