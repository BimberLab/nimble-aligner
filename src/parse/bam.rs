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

            let seq = UMIReader::strip_nonbio_regions(&record.seq().as_bytes()[..]);

            if self.current_umi == read_umi {
                self.current_umi_group
                    .push(seq);
                self.current_cell_barcode = current_cell_barcode.clone();
            } else {
                self.next_umi_group
                    .push(seq);
                self.next_umi = read_umi.clone();
                return Some(true);
            }
        }

        None
    }

    // based on
    // https://assets.ctfassets.net/an68im79xiti/6yYLKUTpokvZvs4pVwnd0a/c40790cd90b7d57bc4e457e670ae3561/CG000207_ChromiumNextGEMSingleCellV_D_J_ReagentKits_v1.1_UG_RevF.pdf,
    // page 76, figure 3.1
    fn strip_nonbio_regions(seq: &[u8]) -> DnaString {
        // Convert seq to string for easy search operations
        let seq = String::from_utf8(seq.to_owned()).unwrap();

        // Find TSO if it exists
        let mut tso_idx = seq.find("TTTCTTATATGGG"); // forward case

        if tso_idx.is_none() {
            tso_idx = seq.find("AAAGAATATACCC"); // reverse case
        };

        // If the TSO exists, strip the front of the sequence, removing the forward primer and any
        // UMI/Cell Barcode metadata, as well as the TSO itself, which has length 13.
        // If it doesn't exist, then there can't be any metadata, so we don't process the pre-cDna
        // portion of the sequence.
        let seq = if tso_idx.is_some() {
            String::from_utf8(seq.as_bytes()[tso_idx.unwrap()+13..].to_vec()).unwrap()
        } else {
            seq
        };

        // Remove the poly-T tail if it exists
        let poly_t_tail_idx = seq.find("TTTTTTTTTTTTTTTTTTTTTTTTT");

        let seq = if poly_t_tail_idx.is_some() {
            String::from_utf8(seq.as_bytes()[..poly_t_tail_idx.unwrap()].to_vec()).unwrap()
        } else {
            seq
        };

        // Find the reverse primer if it exists
        let mut reverse_primer_idx = seq.find("GTACTCTGCGTTGATACCACTGCTT"); // forward case

        if reverse_primer_idx.is_none() {
            reverse_primer_idx = seq.find("CATGAGACGCAACTATGGTGACGAA"); // reverse case
        };

        // If the reverse primer exists, remove it
        let seq = if reverse_primer_idx.is_some() {
            String::from_utf8(seq.as_bytes()[..reverse_primer_idx.unwrap()].to_vec()).unwrap()
        } else {
            seq
        };

        DnaString::from_dna_string(&seq)
    }
}
