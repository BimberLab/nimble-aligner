// Takes a result vector from the filtration pipeline and returns all results above a given score threshold
pub fn threshold_percentage(scores: Vec<(String, f32)>, threshold: f32) -> Vec<(String, f32)> {
    let mut results = Vec::new();

    for (name, score) in scores {
        if score >= threshold {
            results.push((name, score));
        }
    }

    results
}

#[cfg(test)]
mod tests {
    // Tests for threshold_percentage
    // Case where the threshold is 0.0 -- nothing gets thresholded out
    #[test]
    fn threshold_percentage_no_threshold() {
        let scores = vec![
            (String::from("name1"), 50.5),
            (String::from("name2"), 17.2),
            (String::from("name3"), 98.3),
        ];
        let results = super::threshold_percentage(scores, 0.0);

        let mut expected_results: Vec<(String, f32)> = Vec::new();
        expected_results.push((String::from("name1"), 50.5));
        expected_results.push((String::from("name2"), 17.2));
        expected_results.push((String::from("name3"), 98.3));

        assert_eq!(results.len(), 3);
        assert_eq!(results, expected_results);
    }

    // Case where the threshold is 100% -- nothing should get included
    #[test]
    fn threshold_percentage_max_threshold() {
        let scores = vec![
            (String::from("name1"), 50.5),
            (String::from("name2"), 17.2),
            (String::from("name3"), 98.2),
        ];
        let results = super::threshold_percentage(scores, 100.0);

        let expected_results: Vec<(String, f32)> = Vec::new();

        assert_eq!(results.len(), 0);
        assert_eq!(results, expected_results);
    }

    // Case where the threshold is an arbitrary value -- some should be filtered out
    #[test]
    fn threshold_percentage_half_threshold() {
        let scores = vec![
            (String::from("name1"), 50.5),
            (String::from("name2"), 17.2),
            (String::from("name3"), 98.3),
        ];
        let results = super::threshold_percentage(scores, 25.0);

        let mut expected_results: Vec<(String, f32)> = Vec::new();
        expected_results.push((String::from("name1"), 50.5));
        expected_results.push((String::from("name3"), 98.3));

        assert_eq!(results.len(), 2);
        assert_eq!(results, expected_results);
    }
}
