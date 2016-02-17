//! Small bioinformatics library for Rust

use std::fs::File;
use std::io::Read;

/// Represents a base
#[repr(u8)]
#[derive(Debug)]
#[derive(Copy, Clone)]
#[derive(Eq, PartialEq)]
pub enum Base { A = 0, C = 1, G = 2, T = 3 }

impl Base {
  /// Complement of a DNA base
  ///
  /// # Examples
  /// ```
  /// let adenine = hemoglobin::Base::A;
  /// let thyamine = adenine.complement();
  ///
  /// assert_eq!(thyamine, hemoglobin::Base::T);
  /// ```
  pub fn complement(&self) -> Self {
    match *self {
      Base::A => Base::T,
      Base::T => Base::A,
      Base::C => Base::G,
      Base::G => Base::C
    }
  }
  
  /// Convert a u64 number to a Base
  ///
  /// # Examples
  /// ```
  /// let a_base = hemoglobin::Base::from_u64(1).unwrap();
  ///
  /// assert_eq!(a_base, hemoglobin::Base::C);
  /// ```
  pub fn from_u64(num: u64) -> Option<Self> {
    match num {
      0 => Some(Base::A),
      1 => Some(Base::C),
      2 => Some(Base::G),
      3 => Some(Base::T),
      _ => None
    }
  }
}

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

#[derive(Debug)]
#[derive(Copy, Clone)]
#[derive(Eq, PartialEq)]
#[allow(non_snake_case)] // A, C, G, T should be capitalized
pub struct BaseCount {
  pub A: u64,
  pub C: u64,
  pub G: u64,
  pub T: u64
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
  
  /// Convert a bitfield in a u64 to a Sequence
  ///
  /// # Examples
  /// ```
  /// let a_sequence = hemoglobin::Sequence::from_u64(0xBA5E500000000000, 10);
  ///
  /// assert_eq!(a_sequence, hemoglobin::Sequence::new("GTGGCCTGCC"));
  /// ```
  pub fn from_u64(num: u64, length: u8) -> Self {
    let mut result: Vec<Base> = Vec::new();
    
    let end = std::cmp::min(length, 32);
    
    for i in 0..end {
      result.push(Base::from_u64(num << 2*i >> 62).unwrap()); // Result guaranteed to exist for numbers <= 3
    }
    
    return Sequence(result);
  }
  
  /// Retrieve a subsequence, encoded as a bitfield in a u64. Limited to 32 bases
  ///
  /// # Examples
  /// ```
  /// let a_sequence = hemoglobin::Sequence::new("TACGATCTAGTCTAGGATC");
  /// let numbers = a_sequence.subsequence_as_u64(4, 7).unwrap();
  ///
  /// assert_eq!(numbers, 0x372C000000000000);
  /// ```
  pub fn subsequence_as_u64(&self, start: usize, length: usize) -> Option<u64> {
    if length > 32 || start + length > self.0.len() {
      return None;
    }
    
    let mut result: u64 = 0;
    
    for i in 0..length {
      result = result | ((self.0[start + i] as u64) << (62 - 2*i));
    }
    
    return Some(result);
  }
  
  /// Count bases in a sequence
  ///
  /// # Examples
  /// ```
  /// let a_sequence = hemoglobin::Sequence::new("TACGATCTAGTCTAGGATC");
  /// let bases = a_sequence.count_bases();
  /// let adenines = bases.A;
  ///
  /// assert_eq!(adenines, 5);
  /// ```
  pub fn count_bases(&self) -> BaseCount {
    let mut result = BaseCount { A: 0, C: 0, G: 0, T: 0 };
    
    for base in &self.0 {
      match *base {
        Base::A => result.A += 1,
        Base::C => result.C += 1,
        Base::G => result.G += 1,
        Base::T => result.T += 1
      }
    }
    
    return result;
  }
  
  /// Reverse sequence and take its complement
  ///
  /// # Examples
  /// ```
  /// let a_sequence = hemoglobin::Sequence::new("TACGATCTAGTCTAGGATC");
  /// let reverse_complement = a_sequence.reverse_complement();
  ///
  /// assert_eq!(reverse_complement, hemoglobin::Sequence::new("GATCCTAGACTAGATCGTA"));
  /// ```
  pub fn reverse_complement(&self) -> Self {
    let mut result: Vec<Base> = Vec::new();
    
    for i in (0..self.0.len()).rev() {
      result.push(self.0[i].complement());
    }
    
    return Sequence(result);
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
  fn base_complements() {
    assert_eq!(Base::A.complement(), Base::T);
    assert_eq!(Base::T.complement(), Base::A);
    assert_eq!(Base::C.complement(), Base::G);
    assert_eq!(Base::G.complement(), Base::C);
  }
  
  #[test]
  fn base_from_u64() {
    assert_eq!(Base::from_u64(0), Some(Base::A));
    assert_eq!(Base::from_u64(1), Some(Base::C));
    assert_eq!(Base::from_u64(2), Some(Base::G));
    assert_eq!(Base::from_u64(3), Some(Base::T));
    assert_eq!(Base::from_u64(4), None);
  }
  
  #[test]
  fn sequence_reverse_complement() {
    let a_sequence = Sequence::new("TACGATCTAGTCTAGGATC");
    
    assert_eq!(a_sequence.reverse_complement(), Sequence::new("GATCCTAGACTAGATCGTA"));
  }
  
  #[test]
  fn sequence_from_u64() {
    assert_eq!(Sequence::from_u64(0xBA5E500000000000, 10), Sequence::new("GTGGCCTGCC"));
    assert_eq!(Sequence::from_u64(0xBA5E500000000000, 7 ), Sequence::new("GTGGCCT"));
    assert_eq!(Sequence::from_u64(0xBA5E500000000000, 12), Sequence::new("GTGGCCTGCCAA"));
    assert_eq!(Sequence::from_u64(0xBA5E500000000000, 99), Sequence::new("GTGGCCTGCCAAAAAAAAAAAAAAAAAAAAAA"));
  }
  
  #[test]
  fn subsequence_as_u64() {
    let a_sequence = Sequence::new("TACGATCTAGT");
    assert_eq!(a_sequence.subsequence_as_u64(4, 7), Some(0x372C000000000000));
    assert_eq!(a_sequence.subsequence_as_u64(4, 0), Some(0));
    assert_eq!(a_sequence.subsequence_as_u64(4, 8), None);
  }
  
  #[test]
  fn count_bases() {
    assert_eq!(Sequence::new("").count_bases(), BaseCount { A: 0, C: 0, G: 0, T: 0 });
    assert_eq!(Sequence::new("TACGATCTAGTCTAGGATC").count_bases(), BaseCount { A: 5, C: 4, G: 4, T: 6 });
  }
  
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

