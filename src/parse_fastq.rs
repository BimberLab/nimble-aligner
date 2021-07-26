use std::path;
use bio::io::fastq;
use std::io::{Error, ErrorKind};
use debruijn::dna_string::DnaString;
use unwrap::unwrap;

pub type DNAStringIter = impl Iterator<Item = Result<DnaString, Error>>;


// Takes the path to a fastq.gz file and returns an error-checked iterator of the DnaStrings of the file
pub fn get_error_checked_fastq_readers(file_path: &str) -> (DNAStringIter, DNAStringIter) {
    (
        get_error_checked_fastq_reader(file_path),
        get_error_checked_fastq_reader(file_path),
    )
}

fn get_error_checked_fastq_reader(file_path: &str) -> DNAStringIter {
    let (reader, _) = unwrap!(
        niffler::from_path(path::Path::new(file_path)),
        "Error -- could not determine compression format for {}",
        file_path
    );

    fastq::Reader::new(reader)
        .records()
        .map(|record| match record {
            Ok(rec) => Ok(DnaString::from_acgt_bytes(rec.seq())),
            _ => Err(Error::new(
                ErrorKind::InvalidData,
                "Unable to read sequence",
            )),
        })
}
