use std::env;
use std::path;
use bio::io::fasta;

fn main() {
  let args: Vec<String> = env::args().collect();
  let reader = fasta::Reader::from_file(path::Path::new(&args[1])).unwrap();

  for record in reader.records() {
    let record = record.unwrap();
    print!("{}\n", record.id());
  }
}
