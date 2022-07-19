use debruijn::dna_string::DnaString;
use rust_htslib::{bam, bam::record::Aux, bam::Read, bam::Reader};

pub struct UMIReader {
    reader: bam::Reader,
    pub current_umi_group: Vec<DnaString>,
    pub current_metadata_group: Vec<(u8, String, String, bool)>,
    pub current_umi: String,
    pub current_cell_barcode: String,
    pub next_umi_group: Vec<DnaString>,
    pub next_metadata_group: Vec<(u8, String, String, bool)>,
    next_umi: String,
    next_cell_barcode: String
}

impl UMIReader {
    pub fn new(file_path: &str) -> UMIReader {
        UMIReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi_group: Vec::new(),
            current_metadata_group: Vec::new(),
            current_umi: String::new(),
            current_cell_barcode: String::new(),
            next_umi_group: Vec::new(),
            next_metadata_group: Vec::new(),
            next_umi: String::new(),
            next_cell_barcode: String::new()
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
        self.current_metadata_group = self.next_metadata_group.clone();
        self.current_umi = self.next_umi.clone();
        self.current_cell_barcode = self.next_cell_barcode.clone();
        self.next_umi_group.clear();
        self.next_metadata_group.clear();
        self.next_umi.clear();
        self.next_cell_barcode.clear();

        for r in self.reader.records() {
            let mut record = r.unwrap();

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
            let mapq = record.mapq();
            let orientation = String::from(record.strand().strand_symbol());
            let pair = match record.is_secondary() {
                true => String::from("T"),
                false => String::from("F")
            };
            let rev_comp = record.is_reverse();

            if self.current_umi == read_umi {
                self.current_umi_group
                    .push(seq);
                self.current_metadata_group
                    .push((mapq, orientation, pair, rev_comp));
                self.current_cell_barcode = current_cell_barcode.clone();
            } else {
                self.next_umi_group
                    .push(seq);
                self.next_metadata_group
                    .push((mapq, orientation, pair, rev_comp));
                self.next_umi = read_umi.clone();
                self.next_cell_barcode = current_cell_barcode;
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

        // Remove the poly-T/poly-A tail if it exists
        let mut poly_tail_idx = seq.find("TTTTTTTTTTTTTTTTTTT");

        if poly_tail_idx.is_none() {
            poly_tail_idx = seq.find("AAAAAAAAAAAAAAAAAAA");
        };

        let seq = if poly_tail_idx.is_some() {
            String::from_utf8(seq.as_bytes()[..poly_tail_idx.unwrap()].to_vec()).unwrap()
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

#[cfg(test)]
mod tests {
    #[test]
    fn strip_tso_forward() {
        let expected_results = String::from("CCCC");
        let input = "TTTCTTATATGGGCCCC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_tso_forward_and_meta() {
        let expected_results = String::from("CCCC");
        let input = "GGGGGGGGGGGGAGAGGAGAGACCACACACATTTCTTATATGGGCCCC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_tso_reverse() {
        let expected_results = String::from("CCCC");
        let input = "AAAGAATATACCCCCCC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_tso_reverse_and_meta() {
        let expected_results = String::from("CCCC");
        let input = "TTTTTTTATATATAAGAGAGAGAGAGAGAGAGAAAGAATATACCCCCCC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }
    
    #[test]
    fn strip_tail_t() {
        let expected_results = String::from("CCCC");
        let input = "CCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }
    
    #[test]
    fn strip_tail_a() {
        let expected_results = String::from("CCCC");
        let input = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_tail_t_and_meta() {
        let expected_results = String::from("CCCC");
        let input = "CCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGAGAGAGAGAG".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_tail_a_and_meta() {
        let expected_results = String::from("CCCC");
        let input = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAGAGAGAG".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_forward_all() {
        let expected_results = String::from("CCCC");
        let input = "GGGGGGGGGGGGAGAGGAGAGACCACACACATTTCTTATATGGGCCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGAGAGAGAGAG".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_reverse_all() {
        let expected_results = String::from("CCCC");
        let input = "AGAGGAGAGAGTATATATATATTAAAAGAATATACCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGAGAGAGGA".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_reverse_primer_forward() {
        let expected_results = String::from("CCCC");
        let input = "CCCCGTACTCTGCGTTGATACCACTGCTT".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_reverse_primer_reverse() {
        let expected_results = String::from("CCCC");
        let input = "CCCCCATGAGACGCAACTATGGTGACGAA".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_nothing_one() {
        let expected_results = String::from("CCCC");
        let input = "CCCC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }

    #[test]
    fn strip_nothing_two() {
        let expected_results = String::from("CAGACTAGCTAGCTAGCTACGCTACGACTAGCGCATCGAGAGGGCATAGCTCTAGCTACTAC");
        let input = "CAGACTAGCTAGCTAGCTACGCTACGACTAGCGCATCGAGAGGGCATAGCTCTAGCTACTAC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }
}
