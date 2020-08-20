
// Takes a given alignment and returns it if it's above the given match threshold and has an equivalence class
pub fn filter_by_alignment_score(score: usize, equiv_class: Vec<u32>, match_threshold: usize) -> Option<(Vec<u32>, usize)>{
  if score >= match_threshold && !equiv_class.is_empty() { Some((equiv_class, score)) } else { None }
}