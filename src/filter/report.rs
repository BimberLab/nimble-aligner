use std::io::Read;
use std::slice::Iter;
use csv::StringRecordsIntoIter;

// Collection for managing the state of the grouping algorithm
struct GroupCollapseState {
  results: Vec<(String, i32)>,
  curr_group: String,
  curr_score: i32
}

/* Takes a list of references sorted by grouping string (e.g. by lineage) and a vector of the associated
 * scores. Collapses the reference score list by merging the scores for by group
 * Returns a score vector with a length equal to the number of groups in the references
 */
pub fn collapse_results_by_lineage<R: Read>(mut reference_library: StringRecordsIntoIter<R>, mut scores: Iter<i32>) -> Vec<(String, i32)> {
  const GROUP_COLUMN: usize = 4;
  
  // Initialize group collapse state manager to the first element of the reference and score iterators
  let mut group_collapse_state = GroupCollapseState {
    results: Vec::new(),
    curr_group: reference_library.next()
                                 .expect("Error -- reference library empty: ")
                                 .expect("Error -- cannot parse first element of reference library: ")[GROUP_COLUMN].to_string(),
    curr_score: *scores.next().expect("Error -- empty scores list: ")
  };

  for (i, reference) in reference_library.enumerate() {
    if reference.is_err() {
      println!("Warning -- could not read reference {}", i);
      continue;
    }

    let reference = reference.unwrap();
    let next_score = scores.next().expect("Error -- scores and reference lists are not the same length");

    /* If the group hasn't changed, add the score to the current score. Otherwise, push the current 
     * score to the results vector, change to the new group and the new starting score, and continue
     */
    if reference[GROUP_COLUMN] == group_collapse_state.curr_group {
      group_collapse_state.curr_score += next_score;
    } else {
      group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
      group_collapse_state.curr_group = reference[GROUP_COLUMN].to_string();
      group_collapse_state.curr_score = *next_score;
    }
  }

  // Add the score of the remaining group to the score vector
  group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
  group_collapse_state.results
}


// Takes a result vector from the filtration pipeline and returns all results above a given score threshold
pub fn threshold_percentage(scores: Vec<(String, f32)>, threshold: f32) -> Vec<(String, f32)> {
  let mut results = Vec::new();

  for (name, score) in scores {
    if score > threshold {
      results.push((name, score));
    }
  }

  results
}


#[cfg(test)]
mod tests {
    use crate::utils;

    // Tests for collapse_results_by_lineage
    // Case where each reference has a unique group -- filtering should do nothing
    #[test]
    fn collapse_all_unique_groups() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest14
test20\ttest21\ttest22\ttest23\ttest24
test30\ttest31\ttest32\ttest33\ttest34
test40\ttest41\ttest42\ttest43\ttest44";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![100, 200, 50, 2000];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, i32)> = Vec::new();
      expected_results.push((String::from("test14"), 100));
      expected_results.push((String::from("test24"), 200));
      expected_results.push((String::from("test34"), 50));
      expected_results.push((String::from("test44"), 2000));

      assert_eq!(results.len(), 4);
      assert_eq!(results, expected_results);
    }

    // Case where each reference belongs to the same group -- filtering should collapse to a single group
    #[test]
    fn collapse_to_single_group() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest
test20\ttest21\ttest22\ttest23\ttest
test30\ttest31\ttest32\ttest33\ttest
test40\ttest41\ttest42\ttest43\ttest";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![100, 500, 200, 300];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, i32)> = Vec::new();
      expected_results.push((String::from("test"), 1100));

      assert_eq!(results.len(), 1);
      assert_eq!(results, expected_results);
    }

    // Case where there are several groups of the same size
    #[test]
    fn collapse_multiple_groups_same_length() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest1
test20\ttest21\ttest22\ttest23\ttest1
test30\ttest31\ttest32\ttest33\ttest2
test40\ttest41\ttest42\ttest43\ttest2
test50\ttest51\ttest52\ttest53\ttest3
test60\ttest61\ttest62\ttest63\ttest3";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![100, 200, 50, 20, 30, 6000];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, i32)> = Vec::new();
      expected_results.push((String::from("test1"), 300));
      expected_results.push((String::from("test2"), 70));
      expected_results.push((String::from("test3"), 6030));

      assert_eq!(results.len(), 3);
      assert_eq!(results, expected_results);
    }

    // Case where there are several groups, each of different sizes
    #[test]
    fn collapse_multiple_groups_varying_length() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest1
test20\ttest21\ttest22\ttest23\ttest1
test30\ttest31\ttest32\ttest33\ttest1
test40\ttest41\ttest42\ttest43\ttest2
test50\ttest51\ttest52\ttest53\ttest2
test60\ttest61\ttest62\ttest63\ttest2
test70\ttest71\ttest72\ttest73\ttest2
test80\ttest81\ttest82\ttest83\ttest2
test90\ttest91\ttest92\ttest93\ttest3
test100\ttest101\ttest102\ttest103\ttest3
test110\ttest111\ttest112\ttest113\ttest4";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![50, 60, 20, 1000, 200, 300, 100, 200, 200, 300, 1000];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, i32)> = Vec::new();
      expected_results.push((String::from("test1"), 130));
      expected_results.push((String::from("test2"), 1800));
      expected_results.push((String::from("test3"), 500));
      expected_results.push((String::from("test4"), 1000));

      assert_eq!(results.len(), 4);
      assert_eq!(results, expected_results);
    }


    // Tests for threshold_percentage
    // Case where the threshold is 0.0 -- nothing gets thresholded out
    #[test]
    fn threshold_percentage_no_threshold() {
      let scores = vec![(String::from("name1"), 50.5), (String::from("name2"), 17.2), (String::from("name3"), 98.3)];
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
      let scores = vec![(String::from("name1"), 50.5), (String::from("name2"), 17.2), (String::from("name3"), 98.2)];
      let results = super::threshold_percentage(scores, 100.0);

      let expected_results: Vec<(String, f32)> = Vec::new();

      assert_eq!(results.len(), 0);
      assert_eq!(results, expected_results);
    }

    // Case where the threshold is an arbitrary value -- some should be filtered out
    #[test]
    fn threshold_percentage_half_threshold() {
      let scores = vec![(String::from("name1"), 50.5), (String::from("name2"), 17.2), (String::from("name3"), 98.3)];
      let results = super::threshold_percentage(scores, 25.0);

      let mut expected_results: Vec<(String, f32)> = Vec::new();
      expected_results.push((String::from("name1"), 50.5));
      expected_results.push((String::from("name3"), 98.3));

      assert_eq!(results.len(), 2);
      assert_eq!(results, expected_results);

    }
}
