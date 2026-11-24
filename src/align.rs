use std::collections::HashMap;

use crate::seq::{Cigar, CigarOp, revcomp};

const SCORE_MATCH: i32 = 1;
const SCORE_MISMATCH: i32 = -4;
const GAP_OPEN: i32 = -6;
const GAP_EXT: i32 = -1;
const NEG_INF: i32 = -1_000_000;

#[inline]
fn encode_base(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Compact k-mer seed index for a local target or candidate allele sequence (typically 0.5–5 kb).
#[derive(Debug, Clone)]
pub struct LocalIndex {
    pub seq: Vec<u8>,
    k_long: usize,
    long_seeds: HashMap<u32, Vec<usize>>,
    k_short: usize,
    short_seeds: Vec<Vec<usize>>, // flat 4^k table for k=6 (4096 buckets)
}

impl LocalIndex {
    pub fn new(seq: Vec<u8>) -> Self {
        let k_long = 11usize;
        let k_short = 6usize;

        let mut long_seeds: HashMap<u32, Vec<usize>> = HashMap::with_capacity(seq.len());
        if seq.len() >= k_long {
            let mask = (1u32 << (2 * k_long)) - 1;
            let mut key = 0u32;
            let mut valid = 0usize;
            for (idx, &b) in seq.iter().enumerate() {
                if let Some(code) = encode_base(b) {
                    key = ((key << 2) | code) & mask;
                    valid += 1;
                    if valid >= k_long {
                        long_seeds.entry(key).or_default().push(idx + 1 - k_long);
                    }
                } else {
                    key = 0;
                    valid = 0;
                }
            }
        }

        let mut short_seeds = vec![Vec::new(); 1usize << (2 * k_short)];
        if seq.len() >= k_short {
            let mask = (1u32 << (2 * k_short)) - 1;
            let mut key = 0u32;
            let mut valid = 0usize;
            for (idx, &b) in seq.iter().enumerate() {
                if let Some(code) = encode_base(b) {
                    key = ((key << 2) | code) & mask;
                    valid += 1;
                    if valid >= k_short {
                        short_seeds[key as usize].push(idx + 1 - k_short);
                    }
                } else {
                    key = 0;
                    valid = 0;
                }
            }
        }

        Self {
            seq,
            k_long,
            long_seeds,
            k_short,
            short_seeds,
        }
    }

    /// Collect candidate diagonals `d = r_pos - q_pos` sorted by seed vote count descending.
    fn top_diagonals(&self, query: &[u8], use_short: bool, max_cands: usize) -> Vec<isize> {
        let k = if use_short { self.k_short } else { self.k_long };
        if query.len() < k {
            return Vec::new();
        }
        let mask = (1u32 << (2 * k)) - 1;
        let mut diag_votes: HashMap<isize, usize> = HashMap::new();
        let mut key = 0u32;
        let mut valid = 0usize;

        for (idx, &b) in query.iter().enumerate() {
            if let Some(code) = encode_base(b) {
                key = ((key << 2) | code) & mask;
                valid += 1;
                if valid >= k {
                    let q_pos = idx + 1 - k;
                    let hits: Option<&[usize]> = if use_short {
                        let v = &self.short_seeds[key as usize];
                        if v.is_empty() {
                            None
                        } else {
                            Some(v.as_slice())
                        }
                    } else {
                        self.long_seeds.get(&key).map(|v| v.as_slice())
                    };
                    if let Some(positions) = hits {
                        // Ignore hyper-repetitive homopolymers on short seeds
                        if use_short && positions.len() > 24 {
                            continue;
                        }
                        for &r_pos in positions {
                            let d = r_pos as isize - q_pos as isize;
                            // Bin diagonals by 4 bp to tolerate small indels
                            let bin = d.div_euclid(4) * 4;
                            *diag_votes.entry(bin).or_insert(0) += 1;
                        }
                    }
                }
            } else {
                key = 0;
                valid = 0;
            }
        }

        if diag_votes.is_empty() {
            return Vec::new();
        }

        let mut ranked: Vec<(isize, usize)> = diag_votes.into_iter().collect();
        ranked.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let mut chosen = Vec::with_capacity(max_cands);
        for (d, _) in ranked {
            if chosen.iter().all(|&prev: &isize| (prev - d).abs() >= 16) {
                chosen.push(d);
                if chosen.len() >= max_cands {
                    break;
                }
            }
        }
        chosen
    }
}

/// Result of a local alignment of a query slice (already oriented on target `+` strand) to a reference.
#[derive(Debug, Clone)]
pub struct LocalAlignment {
    pub ref_start: usize,   // 0-based inclusive
    pub ref_end: usize,     // 0-based exclusive
    pub query_start: usize, // 0-based inclusive
    pub query_end: usize,   // 0-based exclusive
    pub score: i32,
    pub nm: usize,            // edit distance (mismatches + ins + del)
    pub exact_matches: usize, // exact base matches in aligned region
    pub cigar: Cigar,
}

impl LocalAlignment {
    /// Returns true if the alignment satisfies stringent WT / Mut assignment criteria.
    pub fn is_stringent(
        &self,
        max_ins: usize,
        max_del: usize,
        max_clip: usize,
        max_edit_dist: usize,
    ) -> bool {
        self.cigar.total_ins() <= max_ins
            && self.cigar.total_del() <= max_del
            && self.cigar.max_clip() <= max_clip
            && self.nm <= max_edit_dist
    }
}

/// Result of aligning a read on both `+` and `-` strands.
#[derive(Debug, Clone)]
pub struct OrientedAlignment {
    pub is_reverse: bool,
    pub oriented_seq: Vec<u8>,
    pub aln: LocalAlignment,
}

/// Fast-path ungapped alignment check for clean full-length reads.
fn try_ungapped_fastpath(
    reference: &[u8],
    query: &[u8],
    diag: isize,
    max_mm: usize,
) -> Option<LocalAlignment> {
    if diag < 0 {
        return None;
    }
    let r_start = diag as usize;
    let q_len = query.len();
    if r_start + q_len > reference.len() || q_len < 20 {
        return None;
    }
    let r_slice = &reference[r_start..r_start + q_len];
    let mut mm = 0usize;
    let edge_guard = 5usize.min(q_len / 4);
    for i in 0..q_len {
        if query[i] != r_slice[i] {
            // If mismatch is right near read ends, fall back to full SW so soft-clips aren't missed
            if i < edge_guard || i + edge_guard >= q_len {
                return None;
            }
            mm += 1;
            if mm > max_mm {
                return None;
            }
        }
    }
    let exact_matches = q_len - mm;
    let score = (exact_matches as i32) * SCORE_MATCH + (mm as i32) * SCORE_MISMATCH;
    let mut cigar = Cigar::default();
    cigar.push(CigarOp::Match(q_len));
    Some(LocalAlignment {
        ref_start: r_start,
        ref_end: r_start + q_len,
        query_start: 0,
        query_end: q_len,
        score,
        nm: mm,
        exact_matches,
        cigar,
    })
}

/// Banded local Gotoh (affine-gap) Smith-Waterman alignment around diagonal `d_center = r - q`.
pub fn banded_local_sw(
    reference: &[u8],
    query: &[u8],
    d_center: isize,
    band_half: usize,
) -> Option<LocalAlignment> {
    let m = query.len();
    let n = reference.len();
    if m == 0 || n == 0 {
        return None;
    }

    let w = band_half as isize;
    let band_width = 2 * band_half + 1;
    let rows = m + 1;

    // Flat DP buffers indexed as [i * band_width + col] where col = (j - (i + d_center - w))
    let total_cells = rows * band_width;
    let mut h_prev = vec![0i32; band_width];
    let mut h_curr = vec![0i32; band_width];
    let mut f_prev = vec![NEG_INF; band_width];
    let mut f_curr = vec![NEG_INF; band_width];
    let mut tb = vec![0u8; total_cells];

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=m {
        let j_min_raw = (i as isize) + d_center - w;
        let j_max_raw = (i as isize) + d_center + w;
        let j_start = j_min_raw.max(1) as usize;
        let j_end = j_max_raw.min(n as isize);

        h_curr.fill(0);
        f_curr.fill(NEG_INF);

        if j_end < j_start as isize {
            std::mem::swap(&mut h_prev, &mut h_curr);
            std::mem::swap(&mut f_prev, &mut f_curr);
            continue;
        }
        let j_end = j_end as usize;
        let qb = query[i - 1];
        let row_offset = i * band_width;
        let mut e_left = NEG_INF;

        for j in j_start..=j_end {
            let col = (j as isize - j_min_raw) as usize;
            let rb = reference[j - 1];

            // Previous row (i - 1) has j_min_prev = j_min_raw - 1.
            // So:
            //   (i-1, j-1) has col_prev = (j-1) - (j_min_raw - 1) = col
            //   (i-1, j)   has col_prev = j - (j_min_raw - 1)     = col + 1
            let h_diag_prev = h_prev[col];
            let sub = if qb == rb && qb != b'N' {
                SCORE_MATCH
            } else {
                SCORE_MISMATCH
            };
            let diag_score = h_diag_prev + sub;

            // F[i][j]: insertion in query (consumes query, comes from i-1, j -> col + 1)
            let (f_score, f_ext_bit) = if col + 1 < band_width {
                let open = h_prev[col + 1] + GAP_OPEN + GAP_EXT;
                let ext = f_prev[col + 1] + GAP_EXT;
                if ext > open { (ext, 8u8) } else { (open, 0u8) }
            } else {
                (NEG_INF, 0u8)
            };
            f_curr[col] = f_score;

            // E[i][j]: deletion in query (consumes ref, comes from i, j-1 -> col - 1)
            let (e_score, e_ext_bit) = if col > 0 {
                let open = h_curr[col - 1] + GAP_OPEN + GAP_EXT;
                let ext = e_left + GAP_EXT;
                if ext > open { (ext, 4u8) } else { (open, 0u8) }
            } else {
                (NEG_INF, 0u8)
            };
            e_left = e_score;

            let mut best_cell = 0i32;
            let mut h_dir = 0u8;
            if diag_score > best_cell {
                best_cell = diag_score;
                h_dir = 1;
            }
            if e_score > best_cell {
                best_cell = e_score;
                h_dir = 2;
            }
            if f_score > best_cell {
                best_cell = f_score;
                h_dir = 3;
            }

            h_curr[col] = best_cell;
            tb[row_offset + col] = h_dir | e_ext_bit | f_ext_bit;

            if best_cell > best_score {
                best_score = best_cell;
                best_i = i;
                best_j = j;
            }
        }

        std::mem::swap(&mut h_prev, &mut h_curr);
        std::mem::swap(&mut f_prev, &mut f_curr);
    }

    if best_score <= 0 {
        return None;
    }

    // Traceback
    let mut ops_rev: Vec<CigarOp> = Vec::with_capacity(m + 16);
    let mut i = best_i;
    let mut j = best_j;
    let mut state = 0u8; // 0 = H, 1 = E (Del), 2 = F (Ins)
    let mut exact_matches = 0usize;
    let mut nm = 0usize;

    while i > 0 && j > 0 {
        let j_min_raw = (i as isize) + d_center - w;
        let col_isize = j as isize - j_min_raw;
        if col_isize < 0 || (col_isize as usize) >= band_width {
            break;
        }
        let cell = tb[i * band_width + (col_isize as usize)];

        match state {
            0 => match cell & 3 {
                0 => break,
                1 => {
                    ops_rev.push(CigarOp::Match(1));
                    if query[i - 1] == reference[j - 1] && query[i - 1] != b'N' {
                        exact_matches += 1;
                    } else {
                        nm += 1;
                    }
                    i -= 1;
                    j -= 1;
                }
                2 => {
                    state = 1;
                }
                3 => {
                    state = 2;
                }
                _ => break,
            },
            1 => {
                // Deletion (consumes ref)
                ops_rev.push(CigarOp::Deletion(1));
                nm += 1;
                j -= 1;
                if (cell & 4) == 0 {
                    state = 0;
                }
            }
            2 => {
                // Insertion (consumes query)
                ops_rev.push(CigarOp::Insertion(1));
                nm += 1;
                i -= 1;
                if (cell & 8) == 0 {
                    state = 0;
                }
            }
            _ => break,
        }
    }

    let query_start = i;
    let ref_start = j;
    let query_end = best_i;
    let ref_end = best_j;

    let mut cigar = Cigar::default();
    if query_start > 0 {
        cigar.push(CigarOp::SoftClip(query_start));
    }
    for op in ops_rev.into_iter().rev() {
        cigar.push(op);
    }
    if query_end < m {
        cigar.push(CigarOp::SoftClip(m - query_end));
    }

    Some(LocalAlignment {
        ref_start,
        ref_end,
        query_start,
        query_end,
        score: best_score,
        nm,
        exact_matches,
        cigar,
    })
}

/// Align a query (single strand, already on target `+` strand) to an indexed sequence.
pub fn align_single_strand(
    index: &LocalIndex,
    query: &[u8],
    use_short_seeds: bool,
    band_half: usize,
    min_score: i32,
) -> Option<LocalAlignment> {
    if query.is_empty() {
        return None;
    }

    let mut diags = index.top_diagonals(query, use_short_seeds, 6);
    if diags.is_empty() && !use_short_seeds {
        diags = index.top_diagonals(query, true, 6);
    }
    if diags.is_empty() {
        return None;
    }

    // Fast path for clean ungapped full-length reads
    if !use_short_seeds && query.len() >= 25 {
        for &d in diags.iter().take(2) {
            // Check exact diagonal neighborhood (-2..=2 around binned diagonal)
            for offset in 0..4isize {
                if let Some(aln) = try_ungapped_fastpath(&index.seq, query, d + offset, 3)
                    && aln.score >= min_score
                {
                    return Some(aln);
                }
            }
        }
    }

    let mut best: Option<LocalAlignment> = None;
    for d in diags {
        if let Some(aln) = banded_local_sw(&index.seq, query, d, band_half)
            && aln.score >= min_score
        {
            let replace = match &best {
                None => true,
                Some(prev) => {
                    aln.score > prev.score || (aln.score == prev.score && aln.nm < prev.nm)
                }
            };
            if replace {
                best = Some(aln);
            }
        }
    }
    best
}

/// Align a read on both `+` and `-` strands against an indexed target sequence.
pub fn align_read_both_strands(
    index: &LocalIndex,
    read_seq: &[u8],
    band_half: usize,
    min_score: i32,
) -> Option<OrientedAlignment> {
    let fwd = align_single_strand(index, read_seq, false, band_half, min_score);
    let rc_seq = revcomp(read_seq);
    let rev = align_single_strand(index, &rc_seq, false, band_half, min_score);

    match (fwd, rev) {
        (Some(f), Some(r)) => {
            if f.score >= r.score {
                Some(OrientedAlignment {
                    is_reverse: false,
                    oriented_seq: read_seq.to_vec(),
                    aln: f,
                })
            } else {
                Some(OrientedAlignment {
                    is_reverse: true,
                    oriented_seq: rc_seq,
                    aln: r,
                })
            }
        }
        (Some(f), None) => Some(OrientedAlignment {
            is_reverse: false,
            oriented_seq: read_seq.to_vec(),
            aln: f,
        }),
        (None, Some(r)) => Some(OrientedAlignment {
            is_reverse: true,
            oriented_seq: rc_seq,
            aln: r,
        }),
        (None, None) => None,
    }
}

/// Anchor a soft-clipped tail or head of a primary read alignment onto the same strand of WT ref.
/// Returns `(lloc, rloc, q_start, q_end)` in 0-based half-open coordinates for candidate extension:
/// `ext_seq = ref[0..lloc] + oriented_seq[q_start..q_end] + ref[rloc..]`.
pub fn anchor_softclip(
    wt_index: &LocalIndex,
    oriented_seq: &[u8],
    primary: &LocalAlignment,
    min_clip_len: usize,
    min_itd_len: usize,
    max_itd_len: usize,
) -> Option<(usize, usize, usize, usize)> {
    let lead = primary.cigar.leading_clip();
    let trail = primary.cigar.trailing_clip();

    // Require one dominant clip (>= min_clip_len) and at most a small opposite-end clip (<= 5 bp)
    let (is_leading, clip_slice) = if lead >= min_clip_len && trail <= 5 {
        (true, &oriented_seq[..lead])
    } else if trail >= min_clip_len && lead <= 5 {
        (false, &oriented_seq[oriented_seq.len() - trail..])
    } else {
        return None;
    };

    let diags = wt_index.top_diagonals(clip_slice, true, 12);
    if diags.is_empty() {
        return None;
    }

    type AnchorHit = (i32, usize, (usize, usize, usize, usize));
    let min_clip_score = 9i32;
    let mut best_candidate: Option<AnchorHit> = None;

    for d in diags {
        let Some(caln) = banded_local_sw(&wt_index.seq, clip_slice, d, 16) else {
            continue;
        };
        if caln.score < min_clip_score {
            continue;
        }
        // Clip match ratio threshold >= 0.5 (matching FLT3_ITD_ext)
        if caln.exact_matches * 2 < clip_slice.len() {
            continue;
        }

        let (lloc, rloc, q_start, q_end) = if is_leading {
            // Clip must anchor near the 5' start of the read (allowing <= 2 bp low-qual overhang)
            if caln.query_start > 2 {
                continue;
            }
            (
                caln.ref_start,
                primary.ref_end,
                caln.query_start,
                primary.query_end,
            )
        } else {
            // Clip must anchor near the 3' end of the read (allowing <= 2 bp low-qual overhang)
            if caln.query_end + 2 < clip_slice.len() {
                continue;
            }
            (
                primary.ref_start,
                caln.ref_end,
                primary.query_start,
                primary.query_end + caln.query_end,
            )
        };

        if q_end <= q_start {
            continue;
        }
        let used_q_len = q_end - q_start;
        let ref_span = rloc as isize - lloc as isize;
        let itd_size_isize = used_q_len as isize - ref_span;
        if itd_size_isize < min_itd_len as isize || itd_size_isize > max_itd_len as isize {
            continue;
        }
        let itd_size = itd_size_isize as usize;

        let replace = match best_candidate {
            None => true,
            Some((prev_score, prev_size, _)) => {
                caln.score > prev_score || (caln.score == prev_score && itd_size < prev_size)
            }
        };
        if replace {
            best_candidate = Some((caln.score, itd_size, (lloc, rloc, q_start, q_end)));
        }
    }

    best_candidate.map(|(_, _, coords)| coords)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_ref() -> Vec<u8> {
        // 300 bp non-repetitive DNA sequence
        let s = concat!(
            "TGAATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCT",
            "ACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAG",
            "AGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAAT",
            "CTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTT",
            "TTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGG"
        );
        s.as_bytes().to_vec()
    }

    #[test]
    fn test_clean_alignment_both_strands() {
        let reference = sample_ref();
        let idx = LocalIndex::new(reference.clone());
        let read = &reference[40..160];
        let res = align_read_both_strands(&idx, read, 48, 15).unwrap();
        assert!(!res.is_reverse);
        assert_eq!(res.aln.ref_start, 40);
        assert_eq!(res.aln.ref_end, 160);
        assert_eq!(res.aln.nm, 0);

        let rc_read = revcomp(read);
        let res_rc = align_read_both_strands(&idx, &rc_read, 48, 15).unwrap();
        assert!(res_rc.is_reverse);
        assert_eq!(res_rc.aln.ref_start, 40);
        assert_eq!(res_rc.aln.ref_end, 160);
        assert_eq!(res_rc.aln.nm, 0);
    }

    #[test]
    fn test_short_itd_inline_insertion() {
        let reference = sample_ref();
        let idx = LocalIndex::new(reference.clone());
        // Duplicate 18 bp: ref[90..108]
        let mut mut_allele = Vec::new();
        mut_allele.extend_from_slice(&reference[..108]);
        mut_allele.extend_from_slice(&reference[90..108]);
        mut_allele.extend_from_slice(&reference[108..]);

        // Read spanning across the 18 bp ITD
        let read = &mut_allele[40..178];
        let res = align_read_both_strands(&idx, read, 48, 15).unwrap();
        assert_eq!(res.aln.cigar.total_ins(), 18);
        assert_eq!(res.aln.ref_start, 40);
        assert_eq!(res.aln.ref_end, 160);
    }

    #[test]
    fn test_long_itd_split_clip_anchoring() {
        let reference = sample_ref();
        let idx = LocalIndex::new(reference.clone());
        // Duplicate 120 bp (longer than a 100 bp read!): ref[60..180]
        let mut mut_allele = Vec::new();
        mut_allele.extend_from_slice(&reference[..180]);
        mut_allele.extend_from_slice(&reference[60..180]);
        mut_allele.extend_from_slice(&reference[180..]);

        // 100 bp read crossing junction at pos 180 (55 bp left of junction, 45 bp right)
        let read = &mut_allele[125..225];
        let res = align_read_both_strands(&idx, read, 48, 15).unwrap();
        assert!(!res.is_reverse);
        assert_eq!(res.aln.cigar.trailing_clip(), 45);

        let (lloc, rloc, qs, qe) =
            anchor_softclip(&idx, &res.oriented_seq, &res.aln, 9, 3, 500).unwrap();
        let mut ext = Vec::new();
        ext.extend_from_slice(&reference[..lloc]);
        ext.extend_from_slice(&res.oriented_seq[qs..qe]);
        ext.extend_from_slice(&reference[rloc..]);
        assert_eq!(ext, mut_allele);
    }
}
