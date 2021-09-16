use crate::align;
use crate::reference_library;
use crate::utils;
use reference_library::ReferenceMetadata;

use debruijn::dna_string::DnaString;
use std::io::Error;

/* Takes a list of sequences and optionally reverse sequences, a reference library index, reference library metadata,
 * and an aligner configuration object, and returns a vector of scores and relative match percentages generated from an alignment
 * of the sequences to the reference library. */
pub fn score<'a>(
    sequences: (
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    ),
    reverse_sequences: Option<(
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    )>,
    reference_index: &(align::PseudoAligner, align::PseudoAligner),
    reference_metadata: &ReferenceMetadata,
    align_config: &align::AlignFilterConfig,
    debug_info: Option<&mut align::AlignDebugInfo>
) -> Vec<(Vec<String>, i32)> {
    // Perform filtered pseudoalignment
    let reference_scores = align::score(
        sequences,
        reverse_sequences,
        reference_index,
        reference_metadata,
        align_config,
        debug_info
    );

    // Remove scores below the score threshold
    let reference_scores: Vec<(Vec<String>, i32)> = reference_scores
        .into_iter()
        .filter(|(_, val)| val > &align_config.score_filter)
        .collect();

    utils::sort_score_vector(reference_scores)
}
