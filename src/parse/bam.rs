use super::sorted_bam_reader::SortedBamReader;
use debruijn::dna_string::DnaString;
use rust_htslib::bam::record::Aux;
use crate::ALLOCATOR;

const READ_BLOCK_REPORT_SIZE: usize = 1000000;
const MAX_RECORD_ERROR_REPORT_SIZE: usize = 100;
const CLIP_LENGTH: usize = 13;

pub struct UMIReader {
    reader: SortedBamReader,
    read_counter: usize,
    pub current_umi_group: Vec<DnaString>,
    pub current_metadata_group: Vec<(
        u8,
        String,
        String,
        bool,
        String,
        Vec<u8>,
        Vec<u8>,
        String,
        String,
        String,
        String,
    )>,
    pub current_umi: String,
    pub current_cell_barcode: String,
    pub next_umi_group: Vec<DnaString>,
    pub next_metadata_group: Vec<(
        u8,
        String,
        String,
        bool,
        String,
        Vec<u8>,
        Vec<u8>,
        String,
        String,
        String,
        String,
    )>,
    next_umi: String,
    next_cell_barcode: String,
    terminate_on_error: bool,
    pub number_error_reports: usize,
    pub number_cr_skipped: usize,
    current_iteration_key: String,
    next_iteration_key: String,
}

impl UMIReader {
    pub fn new(file_path: &str, terminate_on_error: bool) -> UMIReader {
        UMIReader {
            reader: SortedBamReader::from_path(file_path),
            read_counter: 0,
            current_umi_group: Vec::new(),
            current_metadata_group: Vec::new(),
            current_umi: String::new(),
            current_cell_barcode: String::new(),
            next_umi_group: Vec::new(),
            next_metadata_group: Vec::new(),
            next_umi: String::new(),
            next_cell_barcode: String::new(),
            terminate_on_error,
            number_error_reports: 0,
            number_cr_skipped: 0,
            current_iteration_key: String::new(),
            next_iteration_key: String::new(),
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
        self.current_iteration_key = self.next_iteration_key.clone();
        self.current_cell_barcode = self.next_cell_barcode.clone();
        self.next_umi_group.clear();
        self.next_metadata_group.clear();
        self.next_umi.clear();
        self.next_cell_barcode.clear();
        self.next_iteration_key.clear();

        loop {
            let r = self.reader.next();

            if r.is_err() {
                return None;
            }

            self.read_counter = self.read_counter + 1;

            if self.read_counter % READ_BLOCK_REPORT_SIZE == 0 && self.read_counter != 0 {
                println!(
                    "Aligned reads {}-{}",
                    self.read_counter - READ_BLOCK_REPORT_SIZE,
                    self.read_counter
                );

                println!(
                    "Current allocated bytes: {}",
                    ALLOCATOR.allocated()
                );
            }

            let record = match r {
                Ok(record) => Ok(record),
                Err(rust_htslib::tpool::Error::BamTruncatedRecord) => {
                    if self.number_error_reports < MAX_RECORD_ERROR_REPORT_SIZE {
                        panic!("{}: Found truncated record", self.number_error_reports);
                    }
                    Err(rust_htslib::tpool::Error::BamTruncatedRecord)
                }
                Err(err) => {
                    if self.number_error_reports < MAX_RECORD_ERROR_REPORT_SIZE {
                        panic!(
                            "{}: Error {:?} when reading record",
                            self.number_error_reports, err
                        );
                        //self.number_error_reports += 1;
                    }
                    Err(err)
                }
            };

            if record.is_err() && self.terminate_on_error {
                panic!("Terminating on error when reading record.")
            } else if record.is_err() && !self.terminate_on_error {
                continue;
            }

            let mut record = record.unwrap();

            let read_umi = if let Ok(Aux::String(corrected)) = record.aux(b"UB") {
                corrected.to_owned()
            } else {
                if let Ok(Aux::String(uncorrected)) = record.aux(b"UR") {
                    uncorrected.to_owned()
                } else {
                    panic!("Error -- Could not read UMI.");
                }
            };

            let current_cell_barcode = if let Ok(Aux::String(corrected)) = record.aux(b"CB") {
                (&corrected[0..corrected.len() - 2]).to_owned()
            } else {
                panic!("Error Read without cell barcode, cannot excise read-mate.");
            };

            let current_iteration_key = read_umi.clone() + current_cell_barcode.as_str();

            if self.current_umi == "" {
                self.current_umi = read_umi.clone();
            }

            if self.current_iteration_key == "" {
                self.current_iteration_key = read_umi.clone() + current_cell_barcode.as_str();
            }

            let seq =
                UMIReader::strip_nonbio_regions(&record.seq().as_bytes()[..], record.is_reverse());
            let mapq = record.mapq();
            let orientation = String::from(record.strand().strand_symbol());
            let pair = match record.is_secondary() {
                true => String::from("T"),
                false => String::from("F"),
            };
            let qname = record.qname().to_vec();

            match self.current_metadata_group.last() {
                Some(record) => {
                    if !(self.current_metadata_group.len() % 2 == 0) {
                        if qname != record.5 {
                            panic!("Error -- QNAME mismatch, mangled read-pair due to incorrect assumptions about sorted order.");
                        }
                    }
                }
                None => {}
            }

            let qual = record.qual().to_vec();
            let rev_comp = record.is_reverse();
            let hit = if let Ok(Aux::String(s)) = record.aux(b"GN") {
                s.to_owned()
            } else {
                String::new()
            };

            let tx = if let Ok(Aux::String(s)) = record.aux(b"TX") {
                s.to_owned()
            } else {
                String::new()
            };

            let an = if let Ok(Aux::String(s)) = record.aux(b"TX") {
                s.to_owned()
            } else {
                String::new()
            };

            if self.current_iteration_key == current_iteration_key {
                self.current_umi_group.push(seq);
                self.current_metadata_group.push((
                    mapq,
                    orientation,
                    pair,
                    rev_comp,
                    hit,
                    qname,
                    qual,
                    tx,
                    read_umi.clone(),
                    current_cell_barcode.clone(),
                    an,
                ));
                self.current_cell_barcode = current_cell_barcode.clone();

                self.current_iteration_key = current_iteration_key;
            } else {
                self.next_umi_group.push(seq);
                self.next_metadata_group.push((
                    mapq,
                    orientation,
                    pair,
                    rev_comp,
                    hit,
                    qname,
                    qual,
                    tx,
                    read_umi.clone(),
                    current_cell_barcode.clone(),
                    an,
                ));
                self.next_umi = read_umi.clone();
                self.next_cell_barcode = current_cell_barcode;
                self.next_iteration_key = current_iteration_key;

                //let duration = start.elapsed();
                //println!("time to push read to NEXT readlist: {:?}", duration);
                return Some(true);
            }
        }
    }

    // based on
    // https://assets.ctfassets.net/an68im79xiti/6yYLKUTpokvZvs4pVwnd0a/c40790cd90b7d57bc4e457e670ae3561/CG000207_ChromiumNextGEMSingleCellV_D_J_ReagentKits_v1.1_UG_RevF.pdf,
    // page 76, figure 3.1
    fn strip_nonbio_regions(seq: &[u8], rev_comp: bool) -> DnaString {
        // Convert seq to string for easy search operations
        //let seq = String::from_utf8(seq.to_owned()).unwrap();

        // Find TSO if it exists
        /*let mut tso_idx = seq.find("TTTCTTATATGGG"); // forward case

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
        };*/

        // Remove the poly-T/poly-A tail if it exists
        /*let mut poly_tail_idx = seq.find("TTTTTTTTTTTTTTTTTTT");

        if poly_tail_idx.is_none() {
            poly_tail_idx = seq.find("AAAAAAAAAAAAAAAAAAA");
        };

        let seq = if poly_tail_idx.is_some() {
            println!("before: {}", seq);
            let after =
                String::from_utf8(seq.as_bytes()[..poly_tail_idx.unwrap()].to_vec()).unwrap();
            println!("after: {}", after);
            after
        } else {
            seq
        };*/

        // Find the reverse primer if it exists
        /*let mut reverse_primer_idx = seq.find("GTACTCTGCGTTGATACCACTGCTT"); // forward case

        if reverse_primer_idx.is_none() {
            reverse_primer_idx = seq.find("CATGAGACGCAACTATGGTGACGAA"); // reverse case
        };

        // If the reverse primer exists, remove it
        let seq = if reverse_primer_idx.is_some() {
            String::from_utf8(seq.as_bytes()[..reverse_primer_idx.unwrap()].to_vec()).unwrap()
        } else {
            seq
        };*/

        //return DnaString::from_acgt_bytes(&seq);

        if seq.len() == 124 {
            if rev_comp {
                DnaString::from_acgt_bytes(&seq[0..seq.len() - CLIP_LENGTH])
            } else {
                DnaString::from_acgt_bytes(&seq[CLIP_LENGTH..])
            }
        } else {
            /*if rev_comp {
                DnaString::from_acgt_bytes(&bio::alphabets::dna::revcomp(seq))
            } else {
                DnaString::from_acgt_bytes(seq)
            }*/
            return DnaString::from_acgt_bytes(&seq);
        }
    }
}

#[cfg(test)]
mod tests {
    /*#[test]
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
        let expected_results =
            String::from("CAGACTAGCTAGCTAGCTACGCTACGACTAGCGCATCGAGAGGGCATAGCTCTAGCTACTAC");
        let input = "CAGACTAGCTAGCTAGCTACGCTACGACTAGCGCATCGAGAGGGCATAGCTCTAGCTACTAC".as_bytes();
        let results = super::UMIReader::strip_nonbio_regions(input).to_string();
        assert_eq!(results, expected_results);
    }*/
}
