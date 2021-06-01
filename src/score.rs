use crate::align;
use crate::reference_library;
use crate::utils;
use reference_library::ReferenceMetadata;

/* Takes a list of sequences and optionally reverse sequences, a reference library index, reference library metadata,
 * and an aligner configuration object, and returns a vector of scores and relative match percentages generated from an alignment
 * of the sequences to the reference library. */
pub fn score(
    sequences: (utils::DNAStringIter, utils::DNAStringIter),
    reverse_sequences: Option<(utils::DNAStringIter, utils::DNAStringIter)>,
    reference_index: (align::PseudoAligner, align::PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: align::AlignFilterConfig,
) -> Vec<(Vec<String>, i32)>
{
    // Perform filtered pseudoalignment
    let reference_scores = align::score(
        sequences,
        reverse_sequences,
        reference_index,
        &reference_metadata,
        &align_config,
    );

    // Remove scores below the score threshold
    let reference_scores: Vec<(Vec<String>, i32)> = reference_scores
        .into_iter()
        .filter(|(_, val)| val > &align_config.score_filter)
        .collect();

    utils::sort_score_vector(reference_scores)
}
