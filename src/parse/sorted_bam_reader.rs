use rust_htslib::{bam, bam::record::Aux, bam::Read, bam::Reader, bam::Record};

pub struct SortedBamReader {
    reader: bam::Reader,
    current_umi: String,
    next_umi: String,
    dna_sorted_buffer: Vec<bam::record::Record>,
    next_records: Vec<bam::record::Record>,
}

impl SortedBamReader {
    pub fn from_path(file_path: &str) -> SortedBamReader {
        SortedBamReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi: String::new(),
            dna_sorted_buffer: Vec::new(),
            next_umi: String::new(),
            next_records: Vec::new(),
        }
    }

    fn fill_buffer(&mut self) {
        println!("filling buffer and sorting");
        self.dna_sorted_buffer.clear();
        self.dna_sorted_buffer.append(&mut self.next_records);
        self.current_umi = self.next_umi.clone();

        for r in self.reader.records() {
            let record = match r {
                Ok(record) => record,
                Err(_) => {
                    continue;
                }
            };

            match record.aux(b"CB") {
                Err(_) => continue,
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

            if self.current_umi == "" {
                self.current_umi = read_umi.clone();
            }

            if self.current_umi != read_umi {
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
                return;
            } else {
                self.dna_sorted_buffer.push(record);
            }
        }
    }

    pub fn next(&mut self) -> Result<Record, rust_htslib::tpool::Error> {
        let record = self.dna_sorted_buffer.pop();

        match record {
            Some(r) => return Ok(r),
            None => self.fill_buffer(),
        }

        let record = self.dna_sorted_buffer.pop();

        match record {
            Some(r) => return Ok(r),
            None => return Err(rust_htslib::tpool::Error::BamTruncatedRecord),
        }
    }
}
