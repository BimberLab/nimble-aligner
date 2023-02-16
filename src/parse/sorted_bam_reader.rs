use rust_htslib::{bam, bam::record::Aux, bam::Read, bam::Reader};

pub struct SortedBamReader {
    reader: bam::Reader,
    current_umi: String,
    dna_sorted_buffer: Vec<bam::record::Record>,
}

impl SortedBamReader {
    pub fn from_path(file_path: &str) -> SortedBamReader {
        SortedBamReader {
            reader: Reader::from_path(file_path).unwrap(),
            current_umi: String::new(),
            dna_sorted_buffer: Vec::new(),
        }
    }

    pub fn records(&mut self) -> impl Iterator<Item = &bam::record::Record> {
        let mut error_reads = Vec::new();

        if !(self.dna_sorted_buffer.is_empty()) {
            return self.dna_sorted_buffer.iter();
        }

        for r in self.reader.records() {
            match r {
                Ok(record) => {
                    if let Ok(Aux::String(_)) = record.aux(b"CB") {
                    } else {
                        continue;
                    };

                    self.dna_sorted_buffer.push(record)
                }
                Err(_) => {
                    error_reads.push(r);
                    continue;
                }
            };

            let read_umi = if let Ok(Aux::String(corrected)) =
                self.dna_sorted_buffer.last().unwrap().aux(b"UB")
            {
                corrected.to_owned()
            } else {
                if let Ok(Aux::String(uncorrected)) =
                    self.dna_sorted_buffer.last().unwrap().aux(b"UR")
                {
                    uncorrected.to_owned()
                } else {
                    panic!("Error -- Could not read UMI.");
                }
            };

            if self.current_umi != read_umi {
                self.dna_sorted_buffer.sort_by(|a, b| {
                    let cb_a = match a.aux(b"UB") {
                        Ok(Aux::String(cb)) => cb,
                        Ok(_) => panic!("Could not read UB"),
                        Err(_) => panic!("Could not read UB."),
                    };

                    let cb_b = match b.aux(b"UB") {
                        Ok(Aux::String(cb)) => cb,
                        Ok(_) => panic!("Could not read UB"),
                        Err(_) => panic!("Could not read UB."),
                    };

                    cb_a.cmp(cb_b)
                });
                return self.dna_sorted_buffer.iter();
            }
        }

        return self.dna_sorted_buffer.iter();
    }
}
