//! Small bioinformatics library for Rust

use std::fs::File;
use std::io::Read;

/// Represents a base
#[repr(u8)]
#[derive(Debug)]
#[derive(Copy, Clone)]
#[derive(Eq, PartialEq)]
pub enum Base { A, T, C, G }

impl std::fmt::Display for Base {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    f.write_str(match *self {
        Base::A => "A",
        Base::T => "T",
        Base::C => "C",
        Base::G => "G"
    })
  }
}

/// Represents a sequence of bases
#[derive(Debug)]
#[derive(Clone)]
#[derive(Eq, PartialEq)]
pub struct Sequence(pub Vec<Base>);

impl std::fmt::Display for Sequence {
  fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    for base in &self.0 {
      try!(base.fmt(f));
    }
    
    return Ok(());
  }
}

impl Sequence {
  /// Constructs a new Sequence
  ///
  /// # Examples
  /// ```
  /// let a_sequence = hemoglobin::Sequence::new("ATCG");
  /// ```
  pub fn new<S: Into<String>>(s: S) -> Self {
    let mut result: Vec<Base> = Vec::new();
    
    for v in s.into().chars() {
      match v {
        'A' => result.push(Base::A),
        'T' => result.push(Base::T),
        'C' => result.push(Base::C),
        'G' => result.push(Base::G),
        _ => {}
      }
    }
    
    return Sequence(result);
  }
  
  /// Constructs a new Sequence from a file
  ///
  /// # Examples
  /// ```
  /// let path = std::path::Path::new("./tests/test_sequence");
  /// let a_sequence = hemoglobin::Sequence::from_file(path);
  /// ```
  pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Self> {
    let mut file = try!(File::open(path));
    let mut contents: String = String::new();
    try!(file.read_to_string(&mut contents));
    
    let result = Sequence::new(contents);
    return Ok(result);
  }
}

/// Find all appearances of a given kmer in a sequence. Must specify if
/// sequence is from a linear circular strand
///
/// # Examples
/// ```
/// let a_sequence = hemoglobin::Sequence::new("ATCGATCG");
/// let a_kmer = hemoglobin::Sequence::new("ATCG");
///
/// let matches = hemoglobin::find_kmers(&a_sequence, &a_kmer, false);
///
/// assert_eq!(matches, vec![0, 4]);
/// ```
pub fn find_kmers(sequence: &Sequence, kmer: &Sequence, circular: bool) -> Vec<u64> {
  let mut result: Vec<u64> = Vec::new();
  
  let end = match circular {
    false => sequence.0.len() - kmer.0.len() + 1,
    true => sequence.0.len()
  };
  
  for i in 0..end {
    let mut matches = true;
    
    for j in 0..kmer.0.len() {
      if sequence.0[(i + j) % sequence.0.len()] != kmer.0[j] {
        matches = false;
      }
    }
    
    if matches {
      result.push(i as u64);
    }
  }
  
  return result;
}

#[cfg(test)]
mod tests {
  use super::*;
  
  #[test]
  fn find_kmers_linear() {
    let a_sequence = Sequence::new("TACGATCTAGTCTAGGATC");
    let a_kmer = Sequence::new("TCTA");
    
    let matches = find_kmers(&a_sequence, &a_kmer, false);
    
    assert_eq!(matches, vec![5, 10]);
  }
  
  #[test]
  fn find_kmers_circular() {
    let a_sequence = Sequence::new("TACGATCTAGTCTAGGATC");
    let a_kmer = Sequence::new("TCTA");
    
    let matches = find_kmers(&a_sequence, &a_kmer, true);
    
    assert_eq!(matches, vec![5, 10, 17]);
  }
}

