use crate::parse_fastq::DNAStringIter;
use std::io::{Error};
use rust_htslib::{bam, bam::Read, bam::Reader, bam::record::Aux};
use debruijn::dna_string::DnaString;

pub struct UMIReader {
    reader: bam::Reader,
    pub current_umi_group: Vec<DnaString>,
    pub current_umi: String,
    next_umi_group: Vec<DnaString>,
    next_umi: String
}


impl UMIReader {
    pub fn new(file_path: &str) -> UMIReader {
        UMIReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi_group: Vec::new(),
            current_umi: String::new(),
            next_umi_group: Vec::new(),
            next_umi: String::new()
        }
    }

    pub fn get_read_iterators(self) -> (DNAStringIter, DNAStringIter) {
        (
            self.current_umi_group.iter().step_by(2).map(|rec| Ok(rec.to_owned())),
            self.current_umi_group.iter().skip(1).step_by(2).map(|rec| Ok(rec.to_owned()))
        )
    }

    fn get_umi_from_bam(&mut self) -> Option<bool>{
        self.current_umi_group = self.next_umi_group.clone();
        self.current_umi = self.next_umi.clone();
        self.next_umi_group.clear();
        self.next_umi.clear();

        for r in self.reader.records() {
            let record = r.unwrap();

            let read_umi = if let Ok(Aux::String(s)) = record.aux(b"UR") {
                s.to_owned()
            } else {
                panic!("Error -- Could not read UMI, internal error.");
            };

            if self.current_umi == "" {
                self.current_umi = read_umi.clone();
            }

            if self.current_umi == read_umi {
                self.current_umi_group.push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
            } else {
                self.next_umi_group.push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
                self.next_umi = read_umi.clone();
                return Some(true)
            }
        }

        None
    }
}
