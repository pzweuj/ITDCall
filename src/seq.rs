use std::fmt;

/// Complement a single ASCII DNA base.
#[inline]
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Return uppercase reverse complement of a DNA byte slice.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(complement_base(b.to_ascii_uppercase()));
    }
    out
}

/// Normalize ASCII DNA bytes to uppercase (non-ACGT becomes 'N').
pub fn normalize_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'A',
            b'C' => b'C',
            b'G' => b'G',
            b'T' => b'T',
            _ => b'N',
        })
        .collect()
}

/// Standard Illumina adapters matching FLT3_ITD_ext (TruSeq R1/R2 + Nextera).
const ILLUMINA_ADAPTERS: &[&[u8]] = &[
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    b"CTGTCTCTTATACACATCT",
];

/// Trim 3' Illumina adapter sequence from a read (and its quality vector if present).
/// Mirrors bbduk `ktrim=r k=23 mink=11 hdist=1`.
pub fn trim_illumina_adapter(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = seq.len();
    if n < 11 {
        return (seq.to_vec(), qual.to_vec());
    }

    let mut cut_pos = n;
    for adapter in ILLUMINA_ADAPTERS {
        let max_k = adapter.len().min(23);
        let min_k = 11usize;

        // 1) Internal full-length (k=23) match with <= 1 mismatch.
        if n >= max_k {
            for pos in 0..=(n - max_k) {
                if pos >= cut_pos {
                    break;
                }
                let mut mm = 0usize;
                for i in 0..max_k {
                    if seq[pos + i] != adapter[i] {
                        mm += 1;
                        if mm > 1 {
                            break;
                        }
                    }
                }
                if mm <= 1 {
                    cut_pos = cut_pos.min(pos);
                    break;
                }
            }
        }

        // 2) 3' suffix partial match down to mink=11 with <= 1 mismatch (0 mismatch if < 14 bp).
        let start_min = n.saturating_sub(max_k.saturating_sub(1));
        let start_max = n.saturating_sub(min_k);
        for pos in start_min..=start_max {
            if pos >= cut_pos {
                break;
            }
            let overlap = n - pos;
            let allowed_mm = if overlap >= 14 { 1 } else { 0 };
            let mut mm = 0usize;
            for i in 0..overlap {
                if seq[pos + i] != adapter[i] {
                    mm += 1;
                    if mm > allowed_mm {
                        break;
                    }
                }
            }
            if mm <= allowed_mm {
                cut_pos = cut_pos.min(pos);
                break;
            }
        }
    }

    let trimmed_seq = seq[..cut_pos].to_vec();
    let trimmed_qual = if qual.len() >= cut_pos {
        qual[..cut_pos].to_vec()
    } else {
        qual.to_vec()
    };
    (trimmed_seq, trimmed_qual)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match(usize),     // M
    Insertion(usize), // I
    Deletion(usize),  // D
    SoftClip(usize),  // S
}

impl CigarOp {
    pub fn len(self) -> usize {
        match self {
            CigarOp::Match(n)
            | CigarOp::Insertion(n)
            | CigarOp::Deletion(n)
            | CigarOp::SoftClip(n) => n,
        }
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Cigar(pub Vec<CigarOp>);

impl Cigar {
    pub fn push(&mut self, op: CigarOp) {
        if op.is_empty() {
            return;
        }
        if let Some(last) = self.0.last_mut() {
            match (last, op) {
                (CigarOp::Match(a), CigarOp::Match(b)) => {
                    *a += b;
                    return;
                }
                (CigarOp::Insertion(a), CigarOp::Insertion(b)) => {
                    *a += b;
                    return;
                }
                (CigarOp::Deletion(a), CigarOp::Deletion(b)) => {
                    *a += b;
                    return;
                }
                (CigarOp::SoftClip(a), CigarOp::SoftClip(b)) => {
                    *a += b;
                    return;
                }
                _ => {}
            }
        }
        self.0.push(op);
    }

    pub fn leading_clip(&self) -> usize {
        match self.0.first() {
            Some(CigarOp::SoftClip(n)) => *n,
            _ => 0,
        }
    }

    pub fn trailing_clip(&self) -> usize {
        match self.0.last() {
            Some(CigarOp::SoftClip(n)) => *n,
            _ => 0,
        }
    }

    pub fn max_clip(&self) -> usize {
        self.leading_clip().max(self.trailing_clip())
    }

    pub fn total_ins(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| match op {
                CigarOp::Insertion(n) => Some(*n),
                _ => None,
            })
            .sum()
    }

    pub fn total_del(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| match op {
                CigarOp::Deletion(n) => Some(*n),
                _ => None,
            })
            .sum()
    }

    pub fn total_match_ops(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(n) => Some(*n),
                _ => None,
            })
            .sum()
    }

    pub fn ref_consumed(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(n) | CigarOp::Deletion(n) => Some(*n),
                _ => None,
            })
            .sum()
    }

    pub fn query_consumed(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::SoftClip(n) => Some(*n),
                _ => None,
            })
            .sum()
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in &self.0 {
            match op {
                CigarOp::Match(n) => write!(f, "{n}M")?,
                CigarOp::Insertion(n) => write!(f, "{n}I")?,
                CigarOp::Deletion(n) => write!(f, "{n}D")?,
                CigarOp::SoftClip(n) => write!(f, "{n}S")?,
            }
        }
        Ok(())
    }
}

/// Parsimonious left-normalization of a variant against a local reference sequence.
/// Input `pos` is 0-based index into `ref_seq` where `ref_allele` begins.
/// Returns `(normalized_pos_0based, ref_allele, alt_allele)` guaranteed to be non-empty
/// and satisfy `ref_seq[..pos] + alt + ref_seq[pos+ref.len()..] == mut_seq`.
pub fn left_normalize_variant(
    ref_seq: &[u8],
    mut pos: usize,
    mut ref_allele: Vec<u8>,
    mut alt_allele: Vec<u8>,
) -> (usize, Vec<u8>, Vec<u8>) {
    // Ensure we start from non-empty alleles by extending right if needed.
    while (ref_allele.is_empty() || alt_allele.is_empty()) && pos + ref_allele.len() < ref_seq.len()
    {
        let next_b = ref_seq[pos + ref_allele.len()];
        ref_allele.push(next_b);
        alt_allele.push(next_b);
    }

    // Tan et al. (2015) standard left-alignment & parsimony:
    // 1. Trim common trailing bases; if either becomes empty, prepend base to the left.
    loop {
        if ref_allele.last() == alt_allele.last()
            && !ref_allele.is_empty()
            && !alt_allele.is_empty()
        {
            if ref_allele.len() > 1 && alt_allele.len() > 1 {
                ref_allele.pop();
                alt_allele.pop();
                continue;
            } else if pos > 0 {
                ref_allele.pop();
                alt_allele.pop();
                pos -= 1;
                let prev_b = ref_seq[pos];
                ref_allele.insert(0, prev_b);
                alt_allele.insert(0, prev_b);
                continue;
            }
        }
        break;
    }

    // 2. Trim common leading bases while keeping at least 1 base in both REF and ALT.
    while ref_allele.len() > 1 && alt_allele.len() > 1 && ref_allele.first() == alt_allele.first() {
        ref_allele.remove(0);
        alt_allele.remove(0);
        pos += 1;
    }

    (pos, ref_allele, alt_allele)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(b"ACGTacgtN"), b"NACGTACGT");
    }

    #[test]
    fn test_trim_adapter() {
        let seq = b"GCTACAGATGGTACAGGTGACCGGCTCCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let qual = vec![b'I'; seq.len()];
        let (tseq, tqual) = trim_illumina_adapter(seq, &qual);
        assert_eq!(tseq, b"GCTACAGATGGTACAGGTGACCGGCTCCTC");
        assert_eq!(tqual.len(), tseq.len());

        // Partial 12 bp 3' adapter
        let seq2 = b"GCTACAGATGGTACAGGTGACCGGCTCCTCAGATCGGAAGAG";
        let qual2 = vec![b'I'; seq2.len()];
        let (tseq2, _) = trim_illumina_adapter(seq2, &qual2);
        assert_eq!(tseq2, b"GCTACAGATGGTACAGGTGACCGGCTCCTC");
    }

    #[test]
    fn test_left_normalize() {
        // ref: GATATATAG, insert "TA" after pos 5 ('A') -> GATATATATAG
        // left-normalized should anchor at pos 0: G -> GAT
        let reference = b"GATATATAG";
        let (pos, r, a) = left_normalize_variant(reference, 5, b"A".to_vec(), b"ATA".to_vec());
        assert_eq!(pos, 0);
        assert_eq!(r, b"G");
        assert_eq!(a, b"GAT");
    }
}
