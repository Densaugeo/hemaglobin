extern crate hemoglobin;

fn main() {
  let path = std::path::Path::new("./test_sequence");
  
  let a_sequence = match hemoglobin::Sequence::from_file(path) {
    Ok(val) => val,
    Err(val) => panic!("Error reading file: {}", path.to_str().unwrap())
  };
  
  println!("{}", a_sequence);
  
  let a_kmer = hemoglobin::Sequence::new("GCAT");
  
  println!("{} {:?}", a_kmer, a_kmer);
  
  let matches: Vec<u64> = hemoglobin::find_kmers(&a_sequence, &a_kmer, false);
  
  println!("{:?}", matches);
}

