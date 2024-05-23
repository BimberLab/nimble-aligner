use crate::align;
use crate::align::{FilterReason, AlignmentOrientation};
use crate::reference_library;
use crate::utils;
use reference_library::Reference;

use debruijn::dna_string::DnaString;
use std::collections::HashMap;
use std::io::Error;

/* Takes a list of sequences and optionally mate sequences, associated sequence metadata, a reference library index, reference library metadata,
 * and an aligner configuration object, and returns a vector of scores and relative match percentages generated from an alignment
 * of the sequences to the reference library. */
pub fn call<'a>(
    sequences: (
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    ),
    mate_sequences: Option<(
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
        Box<dyn Iterator<Item = Result<DnaString, Error>> + 'a>,
    )>,
    per_sequence_metadata: &'a Vec<Vec<String>>,
    reference_index: &(align::PseudoAligner, align::PseudoAligner),
    reference: &Reference,
    aligner_config: &align::AlignFilterConfig,
) -> (
    Vec<(Vec<String>, (i32, Vec<String>, Vec<String>))>,
    Vec<(Vec<String>, String, f64, usize, String)>,
    HashMap<String, ((FilterReason, usize), (FilterReason, usize), (FilterReason, usize), (FilterReason, usize), FilterReason, AlignmentOrientation)>
) {
    // Perform filtered pseudoalignment on the given data, and return the scores sorted by name
    let (reference_scores, alignment_metadata, filter_reasons) = align::get_calls(
        sequences,
        mate_sequences,
        per_sequence_metadata,
        reference_index,
        reference,
        aligner_config,
    );
    (
        utils::sort_score_vector(reference_scores),
        alignment_metadata,
        filter_reasons
    )
}
