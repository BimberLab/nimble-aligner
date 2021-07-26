use rust_htslib::{bam, bam::Read, bam::Reader, bam::record::Aux};
use debruijn::dna_string::DnaString;


struct UMIReader {
    reader: bam::Reader,
    pub current_umi_group: Vec<DnaString>,
    pub current_umi: String,
    next_umi_group: Vec<DnaString>
}


impl UMIReader {
    pub fn new(&self, file_path: &str) -> UMIReader {
        UMIReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi_group: Vec::new(),
            current_umi: String::new(),
            next_umi_group: Vec::new()
        }
    }

    fn get_umi_from_bam(&mut self) {
        let mut current_umi = String::new();
        let mut next_umi_group = Vec::new();

        for r in self.reader.records() {
            let record = r.unwrap();

            let read_umi = if let Ok(Aux::String(s)) = record.aux(b"UR") {
                s.to_owned()
            } else {
                panic!("Error -- Could not read UMI, internal error.");
            };

            if current_umi == "" {
                current_umi = read_umi.clone();
            }

            if current_umi == read_umi {
                self.current_umi_group.push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
            } else {
                next_umi_group.push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
                return
            }
        }
    }
}
