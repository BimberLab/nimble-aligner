use std::fs::File;
use std::path;
use std::slice::Iter;

// Collection for managing the state of the grouping algorithm
struct GroupCollapseState {
  results: Vec<(String, f64)>,
  curr_group: String,
  curr_score: f64
}

/* Takes a list of references sorted by grouping string (e.g. by lineage) and a vector of the associated
 * scores. Collapses the reference score list by merging the scores for by group
 * Returns a score vector with a length equal to the number of groups in the references
 */
pub fn collapse_results_by_lineage(reference_library_path: &str, mut scores: Iter<f64>) -> Vec<(String, f64)> {
  const GROUP_COLUMN: usize = 4;

  // Get iterator to the reference library
  let mut reference_library = csv::ReaderBuilder::new()
      .delimiter(b'\t')
      .from_reader(File::open(path::Path::new(reference_library_path)).unwrap());
  let mut reference_library = reference_library.records();

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
      group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
      group_collapse_state.curr_group = reference[GROUP_COLUMN].to_string();
      group_collapse_state.curr_score = *scores.next().unwrap();
    }
  }

  // Add the score of the remaining group to the score vector
  group_collapse_state.results.push((group_collapse_state.curr_group, group_collapse_state.curr_score));
  group_collapse_state.results
}

// TODO: Write tests
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
