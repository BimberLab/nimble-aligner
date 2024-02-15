use super::sorted_bam_reader::SortedBamReader;
use debruijn::dna_string::DnaString;
use rust_htslib::bam::record::Aux;

const READ_BLOCK_REPORT_SIZE: usize = 1000000;
const MAX_RECORD_ERROR_REPORT_SIZE: usize = 100;
const CLIP_LENGTH: usize = 13;

pub const BAM_FIELDS_TO_REPORT: [&'static str; 37] = [
    "QNAME",
    "QUAL",
    "REVERSE",
    "MATE_REVERSE",
    "PAIRED",
    "PROPER_PAIRED",
    "PAIR_ORIENTATION",
    "UNMAPPED",
    "MATE_UNMAPPED",
    "FIRST_IN_TEMPLATE",
    "LAST_IN_TEMPLATE",
    "STRAND",
    "MAPQ",
    "POS",
    "MATE_POS",
    //"SEQ",
    "SEQ_LEN",
    "INSERT_SIZE",
    "CIGAR",
    "QUALITY_FAILED",
    "SECONDARY",
    "DUPLICATE",
    "SUPPLEMENTARY",
    "NH",
    "HI",
    "AS",
    "GN",
    "TX",
    "AN",
    "nM",
    "fx",
    "RE",
    "CR",
    "CY",
    "CB",
    "UR",
    "UY",
    "UB"
];

pub struct UMIReader {
    reader: SortedBamReader,
    read_counter: usize,
    pub current_umi_group: Vec<DnaString>,
    pub current_metadata_group: Vec<Vec<String>>,
    pub current_umi: String,
    pub current_cell_barcode: String,
    pub next_umi_group: Vec<DnaString>,
    pub next_metadata_group: Vec<Vec<String>>,
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

                /*println!(
                    "Current allocated bytes: {}",
                    ALLOCATOR.allocated()
                );*/
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

            let mut record_fields = Vec::new();
            for field in BAM_FIELDS_TO_REPORT {
                let field_value = if let Ok(Aux::String(s)) = record.aux(field.as_bytes()) {
                    s.to_owned()
                } else {
                    match field {
                        "TID" => record.tid().to_string(),
                        "MATE_TID" => record.mtid().to_string(),
                        "QNAME" => String::from_utf8(record.qname().to_vec()).unwrap_or_else(|e| {
                            eprintln!("Error: {}", e);
                            String::new()
                        }),
                        "QUAL" => String::from_utf8(record.qual().to_vec()).unwrap_or_else(|e| {
                            eprintln!("Error: {}", e);
                            String::new()
                        }),
                        "REVERSE" => record.is_reverse().to_string(),
                        "MATE_REVERSE" => record.is_mate_reverse().to_string(),
                        "PAIRED" => record.is_paired().to_string(),
                        "PROPER_PAIRED" => record.is_proper_pair().to_string(),
                        "PAIR_ORIENTATION" => record.read_pair_orientation().to_string(),
                        "UNMAPPED" => record.is_unmapped().to_string(),
                        "MATE_UNMAPPED" => record.is_mate_unmapped().to_string(),
                        "FIRST_IN_TEMPLATE" => record.is_first_in_template().to_string(),
                        "LAST_IN_TEMPLATE" => record.is_last_in_template().to_string(),
                        "STRAND" => record.strand().strand_symbol().to_owned(),
                        "MAPQ" => record.mapq().to_string(),
                        "POS" => record.pos().to_string(),
                        "MATE_POS" => record.mpos().to_string(),
                        "SEQ" => String::from_utf8(record.seq().as_bytes()).unwrap_or_else(|e| {
                            eprintln!("Error: {}", e);
                            String::new()
                        }),
                        "SEQ_LEN" => record.seq_len().to_string(),
                        "INSERT_SIZE" => record.insert_size().to_string(),
                        "CIGAR" => record.cigar().to_string(),
                        "QUALITY_FAILED" => record.is_quality_check_failed().to_string(),
                        "SECONDARY" => record.is_secondary().to_string(),
                        "DUPLICATE" => record.is_duplicate().to_string(),
                        "SUPPLEMENTARY" => record.is_supplementary().to_string(),
                        _ => { String::new() }
                    }
                };

                record_fields.push(field_value);
            }

            if self.current_iteration_key == current_iteration_key {
                self.current_umi_group.push(seq);
                self.current_metadata_group.push(record_fields);
                self.current_cell_barcode = current_cell_barcode.clone();
                self.current_iteration_key = current_iteration_key;
            } else {
                self.next_umi_group.push(seq);
                self.next_metadata_group.push(record_fields);
                self.next_umi = read_umi.clone();
                self.next_cell_barcode = current_cell_barcode;
                self.next_iteration_key = current_iteration_key;

                return Some(true);
            }
        }
    }

    // based on
    // https://assets.ctfassets.net/an68im79xiti/6yYLKUTpokvZvs4pVwnd0a/c40790cd90b7d57bc4e457e670ae3561/CG000207_ChromiumNextGEMSingleCellV_D_J_ReagentKits_v1.1_UG_RevF.pdf,
    // page 76, figure 3.1
    fn strip_nonbio_regions(seq: &[u8], rev_comp: bool) -> DnaString {
        if seq.len() == 124 {
            if rev_comp {
                DnaString::from_acgt_bytes(&seq[0..seq.len() - CLIP_LENGTH])
            } else {
                DnaString::from_acgt_bytes(&seq[CLIP_LENGTH..])
            }
        } else {
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
