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
) -> (Vec<(Vec<String>, i32)>, Vec<(Vec<String>, String, f64, usize)>) {
    // Perform filtered pseudoalignment
    let (reference_scores, alignment_metadata) = align::score(
        sequences,
        reverse_sequences,
        reference_index,
        reference_metadata,
        align_config,
        debug_info
    );

    (utils::sort_score_vector(reference_scores), alignment_metadata)
}
