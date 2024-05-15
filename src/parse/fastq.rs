use bio::io::fastq;
use debruijn::dna_string::DnaString;
use std::io::{Error, ErrorKind};
use std::path;
use unwrap::unwrap;

// Takes the path to a fastq file and returns two error-checked iterators of the sequences in the file
pub fn get_error_checked_fastq_readers(
    file_path: String,
) -> (
    Box<dyn Iterator<Item = Result<DnaString, Error>>>,
    Box<dyn Iterator<Item = Result<DnaString, Error>>>,
) {
    (
        get_error_checked_fastq_reader(file_path.clone()),
        get_error_checked_fastq_reader(file_path),
    )
}

// Takes the path to a fastq file and returns an error-checked iterator of the sequences in the file
fn get_error_checked_fastq_reader(
    file_path: String,
) -> Box<dyn Iterator<Item = Result<DnaString, Error>>> {
    // Decompress the file at the input path
    let (reader, _) = unwrap!(
        niffler::from_path(path::Path::new(&file_path)),
        "Error -- could not determine compression format for {}",
        file_path
    );

    // Instantiate a fastq reader that maps the records to DnaStrings
    Box::new(
        fastq::Reader::new(reader)
            .records()
            .map(|record| match record {
                Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
                _ => Err(Error::new(
                    ErrorKind::InvalidData,
                    "Unable to read sequence",
                )),
            }),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_fastq_file() {
        let file_path_r1 = "tests/test-sequences/reads/fastq_pipeline_test_r1.fastq";
        let file_path_r2 = "tests/test-sequences/reads/fastq_pipeline_test_r2.fastq";

        let (reader1, reader2) = get_error_checked_fastq_readers(file_path_r1.to_string());

        let seqs1: Vec<_> = reader1.map(|r| r.unwrap().to_string()).collect();
        let seqs2: Vec<_> = reader2.map(|r| r.unwrap().to_string()).collect();

        assert_eq!(seqs1, vec!["ATGCGTAC", "CGTAGCTA"]);
        assert_eq!(seqs2, vec!["ATGCGTAC", "CGTAGCTA"]);

        let (reader3, reader4) = get_error_checked_fastq_readers(file_path_r2.to_string());
        let seqs3: Vec<_> = reader3.map(|r| r.unwrap().to_string()).collect();
        let seqs4: Vec<_> = reader4.map(|r| r.unwrap().to_string()).collect();

        assert_eq!(seqs3, vec!["TACGTCAT", "TAGCTACG"]);
        assert_eq!(seqs4, vec!["TACGTCAT", "TAGCTACG"]);
    }

    #[test]
    #[should_panic(expected = "Error -- could not determine compression format")]
    fn test_file_read_error() {
        let non_existing_file_path = "tests/test-sequences/reads/nonexistent.fastq";
        let _ = get_error_checked_fastq_readers(non_existing_file_path.to_string());
    }

    #[test]
    #[should_panic(expected = "Unable to read sequence")]
    fn test_invalid_sequence_data() {
        let file_path_invalid = "tests/test-sequences/reads/fastq_invalid_data.fastq";
        let (reader1, _) = get_error_checked_fastq_readers(file_path_invalid.to_string());
        let _ = reader1.map(|r| r.unwrap()).collect::<Vec<_>>();
    }
}