use rust_htslib::{bam, bam::record::Aux, bam::Read, bam::Reader, bam::Record};
use std::str;

const TAG_WHITELIST: [&'static str; 1] = ["AAAAAAAAAA"];

pub struct SortedBamReader {
    reader: bam::Reader,
    current_umi: String,
    current_cell_barcode: String,
    next_umi: String,
    next_cell_barcode: String,
    dna_sorted_buffer: Vec<bam::record::Record>,
    next_records: Vec<bam::record::Record>,
    force_bam_paired: bool
}

impl SortedBamReader {
    pub fn from_path(file_path: &str, force_bam_paired: bool) -> SortedBamReader {
        SortedBamReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi: String::new(),
            current_cell_barcode: String::new(),
            dna_sorted_buffer: Vec::new(),
            next_umi: String::new(),
            next_cell_barcode: String::new(),
            next_records: Vec::new(),
            force_bam_paired
        }
    }

    fn fill_buffer(&mut self) {
        self.dna_sorted_buffer.clear();
        self.dna_sorted_buffer.append(&mut self.next_records);
        self.current_umi = self.next_umi.clone();
        //self.current_cell_barcode = self.next_cell_barcode.clone();

        for r in self.reader.records() {
            let record = match r {
                Ok(record) => record,
                Err(_) => {
                    continue;
                }
            };

            if !record.is_paired() && self.force_bam_paired {
                continue;
            }

            match record.aux(b"CB") {
                Err(_) => {
                    continue;
                }
                _ => (),
            };

            let read_umi = if let Ok(Aux::String(corrected)) = record.aux(b"UB") {
                corrected.to_owned()
            } else {
                if let Ok(Aux::String(uncorrected)) = record.aux(b"UR") {
                    uncorrected.to_owned()
                } else {
                    panic!("Error -- Could not read UMI.");
                }
            };

            if read_umi == "AAAAAAAAAA" {
                continue;
            }

            if self.current_umi == "" {
                self.current_umi = read_umi.clone();
            }

            /*let read_cell_barcode = if let Ok(Aux::String(corrected)) = record.aux(b"CB") {
                (&corrected[0..corrected.len() - 2]).to_owned()
            } else {
                panic!("Error Read without cell barcode, cannot excise read-mate.");
            };

            if self.current_cell_barcode == "" {
                self.current_cell_barcode = read_cell_barcode.clone();
            }*/

            if self.current_umi != read_umi {//|| self.current_cell_barcode != read_cell_barcode {
                self.dna_sorted_buffer.sort_by(|a, b| {
                    let cb_a = match a.aux(b"CB") {
                        Ok(Aux::String(cb)) => cb,
                        _ => panic!("Could not read CB"),
                    };

                    let cb_b = match b.aux(b"CB") {
                        Ok(Aux::String(cb)) => cb,
                        _ => panic!("Could not read CB"),
                    };

                    cb_a.cmp(cb_b)
                });

                self.next_records.push(record);
                self.next_umi = read_umi;
                //self.next_cell_barcode = read_cell_barcode;
                return;
            } else {
                self.dna_sorted_buffer.push(record);
            }
        }
    }

    fn add_dummy_paired_reads(&mut self) {
        let mut new_buffer = Vec::new();
        
        for read in &mut self.dna_sorted_buffer {
            let mut modified_read = read.clone();
            modified_read.push_aux(b"SKIP_ALIGN", Aux::String("FALSE")).unwrap();
            new_buffer.push(modified_read);
            
            if !read.is_paired() {
                let mut dummy_duplicate = read.clone();
                dummy_duplicate.push_aux(b"SKIP_ALIGN", Aux::String("TRUE")).unwrap();
                new_buffer.push(dummy_duplicate);
            }
        }
        
        self.dna_sorted_buffer = new_buffer;
    }

    fn filter_paired_reads(&mut self) {
        let mut paired_reads_buffer = Vec::new();
        let mut seen_qnames = std::collections::HashSet::new();
        let mut i = 0;
        while i < self.dna_sorted_buffer.len() {
            if i + 1 < self.dna_sorted_buffer.len() {
                let read1_qname = self.dna_sorted_buffer[i].qname();
                let read2_qname = self.dna_sorted_buffer[i + 1].qname();
                if read1_qname == read2_qname {
                    if self.dna_sorted_buffer[i].is_first_in_template() {
                        paired_reads_buffer.push(self.dna_sorted_buffer[i].clone());
                        paired_reads_buffer.push(self.dna_sorted_buffer[i + 1].clone());
                    } else {
                        paired_reads_buffer.push(self.dna_sorted_buffer[i + 1].clone());
                        paired_reads_buffer.push(self.dna_sorted_buffer[i].clone());
                    }

                    seen_qnames.insert(read1_qname);
                    i += 2;
                } else {
                    println!("Warning: Unpaired qname!");
                    if seen_qnames.contains(&read1_qname) {
                        println!(
                            "Warning: Read with qname '{:?}' has been deleted but was seen before.",
                            String::from_utf8_lossy(read1_qname)
                        );
                    }
                    seen_qnames.insert(read1_qname);
                    i += 1;
                }
            } else {
                break;
            }
        }
        self.dna_sorted_buffer = paired_reads_buffer;
    }

    pub fn next(&mut self) -> Result<Record, rust_htslib::tpool::Error> {
        let record = self.dna_sorted_buffer.pop();

        match record {
            Some(r) => return Ok(r),
            None => {
                self.fill_buffer();
                if !self.force_bam_paired {
                    self.add_dummy_paired_reads();
                }
                self.filter_paired_reads();
                self.dna_sorted_buffer.reverse();
            }
        }

        let record = self.dna_sorted_buffer.pop();

        match record {
            Some(r) => return Ok(r),
            None => return Err(rust_htslib::tpool::Error::BamTruncatedRecord),
        }
    }
}
