use std::io::Read;
use std::slice::Iter;
use csv::StringRecordsIntoIter;

// Collection for managing the state of the grouping algorithm
struct GroupCollapseState {
  results: Vec<(String, f32)>,
  curr_group: String,
  curr_score: f32
}

/* Takes a list of references sorted by grouping string (e.g. by lineage) and a vector of the associated
 * scores. Collapses the reference score list by merging the scores for by group
 * Returns a score vector with a length equal to the number of groups in the references
 */
// TODO: Look at using itertools for this
pub fn collapse_results_by_lineage<R: Read>(mut reference_library: StringRecordsIntoIter<R>, mut scores: Iter<f32>) -> Vec<(String, f32)> {
  const GROUP_COLUMN: usize = 4;
  const ROUND_FACTOR: f32 = 100000.0;
  
  // Initialize group collapse state manager to the first element of the reference and score iterators
  let mut group_collapse_state = GroupCollapseState {
    results: Vec::new(),
    curr_group: reference_library.next().unwrap().unwrap()[GROUP_COLUMN].to_string(),
    curr_score: *scores.next().unwrap()
  };

  for reference in reference_library {
    let reference = reference.unwrap();

    /* If the group hasn't changed, add the score to the current score. Otherwise, push the current 
     * score to the results vector, change to the new group and the new starting score, and continue
     */
    if reference[GROUP_COLUMN] == group_collapse_state.curr_group {
      group_collapse_state.curr_score += scores.next().unwrap();
    } else {
      // Cut results to 5 digits to remove any imprecision
      group_collapse_state.curr_score = (group_collapse_state.curr_score * ROUND_FACTOR).round() / ROUND_FACTOR;

      group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
      group_collapse_state.curr_group = reference[GROUP_COLUMN].to_string();
      group_collapse_state.curr_score = *scores.next().unwrap();
    }
  }

  // Cut results to 5 digits to remove any imprecision
  group_collapse_state.curr_score = (group_collapse_state.curr_score * ROUND_FACTOR).round() / ROUND_FACTOR;

  // Add the score of the remaining group to the score vector
  group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
  group_collapse_state.results
}


#[cfg(test)]
mod tests {
    use crate::utils;

    // Case where each reference has a unique group
    #[test]
    fn filter_all_unique_groups() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest14
test20\ttest21\ttest22\ttest23\ttest24
test30\ttest31\ttest32\ttest33\ttest34
test40\ttest41\ttest42\ttest43\ttest44";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![0.3, 0.5, 0.4, 0.1];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, f32)> = Vec::new();
      expected_results.push((String::from("test14"), 0.3));
      expected_results.push((String::from("test24"), 0.5));
      expected_results.push((String::from("test34"), 0.4));
      expected_results.push((String::from("test44"), 0.1));

      assert_eq!(results.len(), 4);
      assert_eq!(results, expected_results);
    }

    // Case where each reference belongs to the same group
    #[test]
    fn filter_single_group() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest
test20\ttest21\ttest22\ttest23\ttest
test30\ttest31\ttest32\ttest33\ttest
test40\ttest41\ttest42\ttest43\ttest";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![0.3, 0.5, 0.4, 0.1];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, f32)> = Vec::new();
      expected_results.push((String::from("test"), 1.3));

      assert_eq!(results.len(), 1);
      assert_eq!(results, expected_results);
    }

    // Case where there are several groups of the same size
    #[test]
    fn filter_multiple_groups_same_length() {
      let reference_library_data = "\
header1\theader2\theader3\theader4\theader5
test10\ttest11\ttest12\ttest13\ttest1
test20\ttest21\ttest22\ttest23\ttest1
test30\ttest31\ttest32\ttest33\ttest2
test40\ttest41\ttest42\ttest43\ttest2
test50\ttest51\ttest52\ttest53\ttest3
test60\ttest61\ttest62\ttest63\ttest3";

      let reference_library = utils::get_tsv_reader(reference_library_data.as_bytes());

      let scores = vec![0.3, 0.5, 0.4, 0.1, 0.2, 0.7];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, f32)> = Vec::new();
      expected_results.push((String::from("test1"), 0.8));
      expected_results.push((String::from("test2"), 0.5));
      expected_results.push((String::from("test3"), 0.9));

      assert_eq!(results.len(), 3);
      assert_eq!(results, expected_results);
    }

    // Case where there are several groups, each of different sizes
    #[test]
    fn filter_multiple_groups_varying_length() {
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

      let scores = vec![0.3, 0.5, 0.4, 0.1, 0.2, 0.7, 0.4, 1.2, 1.5, 1.3, 0.1];

      let results = super::collapse_results_by_lineage(reference_library.into_records(), scores.iter());

      let mut expected_results: Vec<(String, f32)> = Vec::new();
      expected_results.push((String::from("test1"), 1.2));
      expected_results.push((String::from("test2"), 2.6));
      expected_results.push((String::from("test3"), 2.8));
      expected_results.push((String::from("test4"), 0.1));

      assert_eq!(results.len(), 4);
      assert_eq!(results, expected_results);
    }
}
