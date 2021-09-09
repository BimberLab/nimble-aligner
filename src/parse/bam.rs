use debruijn::dna_string::DnaString;
use rust_htslib::{bam, bam::record::Aux, bam::Read, bam::Reader};

pub struct UMIReader {
    reader: bam::Reader,
    pub current_umi_group: Vec<DnaString>,
    pub current_umi: String,
    pub current_cell_barcode: String,
    pub next_umi_group: Vec<DnaString>,
    next_umi: String
}

impl UMIReader {
    pub fn new(file_path: &str) -> UMIReader {
        UMIReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi_group: Vec::new(),
            current_umi: String::new(),
            current_cell_barcode: String::new(),
            next_umi_group: Vec::new(),
            next_umi: String::new()
        }
    }

    pub fn next(&mut self) -> bool {
        let mut final_umi = false;

        if self.get_umi_from_bam().is_none() {
            final_umi = true;
        }

        final_umi
    }

    fn get_umi_from_bam(&mut self) -> Option<bool> {
        self.current_umi_group = self.next_umi_group.clone();
        self.current_umi = self.next_umi.clone();
        self.current_cell_barcode.clear();
        self.next_umi_group.clear();
        self.next_umi.clear();

        for r in self.reader.records() {
            let record = r.unwrap();

            let read_umi = if let Ok(Aux::String(s)) = record.aux(b"UR") {
                s.to_owned()
            } else {
                panic!("Error -- Could not read UMI, internal error.");
            };

            let current_cell_barcode = if let Ok(Aux::String(s)) = record.aux(b"CR") {
                s.to_owned()
            } else {
                panic!("Error -- Could not read cell barcode, internal error.");
            };


            if self.current_umi == "" {
                self.current_umi = read_umi.clone();
            }

            if self.current_umi == read_umi {
                self.current_umi_group
                    .push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
                self.current_cell_barcode = current_cell_barcode.clone();
            } else {
                self.next_umi_group
                    .push(DnaString::from_acgt_bytes(&record.seq().as_bytes()[..]));
                self.next_umi = read_umi.clone();
                return Some(true);
            }
        }

        None
    }
}
