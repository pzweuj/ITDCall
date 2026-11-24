use std::collections::HashMap;

use crate::align::{LocalIndex, align_single_strand};
use crate::seq::left_normalize_variant;

/// Compute Levenshtein edit distance between `a` and `b` if it is `<= max_dist`, else `None`.
/// Strips identical prefix and suffix first for O((L_diff) * max_dist) speed.
pub fn banded_edit_distance(a: &[u8], b: &[u8], max_dist: usize) -> Option<usize> {
    if a.len().abs_diff(b.len()) > max_dist {
        return None;
    }

    // Strip common prefix
    let mut prefix = 0usize;
    let min_len = a.len().min(b.len());
    while prefix < min_len && a[prefix] == b[prefix] {
        prefix += 1;
    }
    let a = &a[prefix..];
    let b = &b[prefix..];

    // Strip common suffix
    let mut suffix = 0usize;
    let min_rem = a.len().min(b.len());
    while suffix < min_rem && a[a.len() - 1 - suffix] == b[b.len() - 1 - suffix] {
        suffix += 1;
    }
    let a = &a[..a.len() - suffix];
    let b = &b[..b.len() - suffix];

    if a.is_empty() {
        return (b.len() <= max_dist).then_some(b.len());
    }
    if b.is_empty() {
        return (a.len() <= max_dist).then_some(a.len());
    }

    let m = a.len();
    let n = b.len();
    if m.abs_diff(n) > max_dist {
        return None;
    }

    let inf = max_dist + 1;
    let band = 2 * max_dist + 1;
    let mut prev = vec![inf; band];
    let mut curr = vec![inf; band];

    // Initialize row 0: j in 0..=max_dist -> col = j + max_dist
    for j in 0..=max_dist.min(n) {
        prev[j + max_dist] = j;
    }

    for i in 1..=m {
        curr.fill(inf);
        let j_min = i.saturating_sub(max_dist);
        let j_max = (i + max_dist).min(n);
        let mut row_min = inf;

        for j in j_min..=j_max {
            let col = (j as isize - i as isize + max_dist as isize) as usize;
            let mut val = inf;

            if j == 0 {
                val = i;
            } else {
                // Substitution / match from (i-1, j-1) -> same col in prev
                let cost = usize::from(a[i - 1] != b[j - 1]);
                val = val.min(prev[col] + cost);

                // Deletion from a (i-1, j) -> col + 1 in prev
                if col + 1 < band {
                    val = val.min(prev[col + 1] + 1);
                }
                // Insertion into a (i, j-1) -> col - 1 in curr
                if col > 0 {
                    val = val.min(curr[col - 1] + 1);
                }
            }

            curr[col] = val;
            if val < row_min {
                row_min = val;
            }
        }

        if row_min > max_dist {
            return None;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let final_col = (n as isize - m as isize + max_dist as isize) as usize;
    let dist = prev[final_col];
    (dist <= max_dist).then_some(dist)
}

#[derive(Debug, Clone)]
pub struct CandidateCluster {
    pub consensus_seq: Vec<u8>,
    pub total_support: usize,
    pub itd_size: usize,
}

/// Greedy abundance-sorted clustering of extended candidate ITD sequences (replacing sumaclust).
pub fn cluster_extended_reads(
    ext_counts: &HashMap<Vec<u8>, usize>,
    ref_len: usize,
    cluster_edit_dist: usize,
) -> Vec<CandidateCluster> {
    if ext_counts.is_empty() {
        return Vec::new();
    }

    let mut items: Vec<(&Vec<u8>, usize)> = ext_counts.iter().map(|(s, &c)| (s, c)).collect();
    items.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then_with(|| a.0.len().cmp(&b.0.len()))
            .then_with(|| a.0.cmp(b.0))
    });

    struct RawCluster<'a> {
        centroid: &'a [u8],
        total_support: usize,
        members: Vec<(&'a [u8], usize)>,
    }

    let mut clusters: Vec<RawCluster<'_>> = Vec::new();

    for (seq, count) in items {
        let mut best_match: Option<(usize, usize)> = None; // (cluster_idx, dist)
        for (idx, cl) in clusters.iter().enumerate() {
            if let Some(d) = banded_edit_distance(seq, cl.centroid, cluster_edit_dist) {
                let replace = match best_match {
                    None => true,
                    Some((_, prev_d)) => d < prev_d,
                };
                if replace {
                    best_match = Some((idx, d));
                    if d == 0 {
                        break;
                    }
                }
            }
        }

        if let Some((idx, _)) = best_match {
            clusters[idx].total_support += count;
            clusters[idx].members.push((seq.as_slice(), count));
        } else {
            clusters.push(RawCluster {
                centroid: seq.as_slice(),
                total_support: count,
                members: vec![(seq.as_slice(), count)],
            });
        }
    }

    // Build consensus among same-length members of each cluster (using centroid length)
    let mut out = Vec::with_capacity(clusters.len());
    for cl in clusters {
        let target_len = cl.centroid.len();
        if target_len <= ref_len {
            continue;
        }
        let mut consensus = cl.centroid.to_vec();
        let same_len_count: usize = cl
            .members
            .iter()
            .filter(|(s, _)| s.len() == target_len)
            .count();

        if same_len_count > 1 {
            let mut counts = vec![[0usize; 4]; target_len];
            for &(s, w) in &cl.members {
                if s.len() == target_len {
                    for (pos, &b) in s.iter().enumerate() {
                        match b {
                            b'A' => counts[pos][0] += w,
                            b'C' => counts[pos][1] += w,
                            b'G' => counts[pos][2] += w,
                            b'T' => counts[pos][3] += w,
                            _ => {}
                        }
                    }
                }
            }
            const BASES: [u8; 4] = *b"ACGT";
            for (pos, c) in counts.iter().enumerate() {
                let (best_i, &best_w) = c.iter().enumerate().max_by_key(|&(_, w)| w).unwrap();
                if best_w > 0 {
                    consensus[pos] = BASES[best_i];
                }
            }
        }

        out.push(CandidateCluster {
            itd_size: target_len - ref_len,
            consensus_seq: consensus,
            total_support: cl.total_support,
        });
    }

    out.sort_by(|a, b| {
        b.total_support
            .cmp(&a.total_support)
            .then_with(|| a.itd_size.cmp(&b.itd_size))
    });
    out
}

/// Structural characterization of a candidate ITD allele relative to WT target sequence (0-based).
#[derive(Debug, Clone)]
pub struct CharacterizedItd {
    pub id: usize,
    pub candidate_seq: Vec<u8>,
    pub initial_support: usize,
    /// Net insertion size: `candidate_seq.len() - wt_ref.len()`.
    pub itd_size: usize,
    /// 0-based inclusive start of duplicated reference region on WT ref.
    pub dup_r0: usize,
    /// 0-based exclusive end of duplicated reference region on WT ref.
    pub dup_r1: usize,
    /// Length of duplicated reference segment (`dup_r1 - dup_r0`).
    pub dup_len: usize,
    /// Non-templated inserted sequence at junction (if any).
    pub ins_seq: Vec<u8>,
    /// Left internal mutant breakpoint (0-based boundary in `candidate_seq`: end of left copy).
    pub mutbk_s0: usize,
    /// Right internal mutant breakpoint (0-based boundary in `candidate_seq`: start of right copy).
    pub mutbk_s1: usize,
    /// Number of point mismatches inside the candidate arms.
    pub mismatches: usize,
    /// Normalized 0-based local VCF position on WT ref.
    pub vcf_local_pos: usize,
    /// Normalized VCF REF allele (genomic `+` strand).
    pub vcf_ref: Vec<u8>,
    /// Normalized VCF ALT allele (genomic `+` strand).
    pub vcf_alt: Vec<u8>,
}

/// Forward X-drop anchored ungapped extension from `(0, 0)` along `(wt_ref, cand_seq)`.
/// Returns `(best_len, mismatch_positions)`.
fn extend_left_arm(wt_ref: &[u8], cand_seq: &[u8]) -> (usize, Vec<usize>) {
    let max_len = wt_ref.len().min(cand_seq.len());
    let mut score = 0i32;
    let mut best_score = 0i32;
    let mut best_len = 0usize;
    let mut mm_all = Vec::new();
    let x_drop = 24i32;

    for i in 0..max_len {
        if wt_ref[i] == cand_seq[i] && wt_ref[i] != b'N' {
            score += 1;
            if score >= best_score {
                best_score = score;
                best_len = i + 1;
            }
        } else {
            score -= 4;
            mm_all.push(i);
            if best_score - score > x_drop {
                break;
            }
        }
    }

    mm_all.retain(|&p| p < best_len);
    (best_len, mm_all)
}

/// Backward X-drop anchored ungapped extension from `(wt_ref.len(), cand_seq.len())` going left.
/// Returns `(best_suffix_len, mismatch_offsets_from_end)`.
fn extend_right_arm(wt_ref: &[u8], cand_seq: &[u8]) -> (usize, Vec<usize>) {
    let wt_len = wt_ref.len();
    let cand_len = cand_seq.len();
    let max_len = wt_len.min(cand_len);
    let mut score = 0i32;
    let mut best_score = 0i32;
    let mut best_len = 0usize;
    let mut mm_all = Vec::new();
    let x_drop = 24i32;

    for offset in 0..max_len {
        let rb = wt_ref[wt_len - 1 - offset];
        let cb = cand_seq[cand_len - 1 - offset];
        if rb == cb && rb != b'N' {
            score += 1;
            if score >= best_score {
                best_score = score;
                best_len = offset + 1;
            }
        } else {
            score -= 4;
            mm_all.push(offset);
            if best_score - score > x_drop {
                break;
            }
        }
    }

    mm_all.retain(|&o| o < best_len);
    (best_len, mm_all)
}

/// Characterize a candidate ITD sequence against `wt_ref` and derive its normalized VCF allele.
pub fn characterize_candidate(
    id: usize,
    wt_ref: &[u8],
    wt_index: &LocalIndex,
    cluster: &CandidateCluster,
) -> Option<CharacterizedItd> {
    let cand_seq = &cluster.consensus_seq;
    let wt_len = wt_ref.len();
    let cand_len = cand_seq.len();
    if cand_len <= wt_len {
        return None;
    }
    let itd_size = cand_len - wt_len;

    let (mut left_len, left_mm) = extend_left_arm(wt_ref, cand_seq);
    let (right_suffix, right_mm_off) = extend_right_arm(wt_ref, cand_seq);

    if left_len < 15 || right_suffix < 15 {
        return None;
    }

    let mut left_s1 = left_len;
    let mut left_r1 = left_len;
    let right_s0 = cand_len - right_suffix;
    let right_r0 = wt_len - right_suffix;

    // Resolve microhomology overlap between left arm and right arm (left-shift consistently)
    if right_s0 < left_s1 {
        let overlap = left_s1 - right_s0;
        left_s1 -= overlap;
        left_r1 -= overlap;
        left_len = left_s1;
    }

    let mut ins_seq: Vec<u8> = Vec::new();
    let dup_r0: usize;
    let dup_r1: usize;
    let mutbk_s0: usize = left_s1;
    let mut mutbk_s1: usize = right_s0;

    if right_r0 <= left_r1 {
        // Classic tandem duplication (pure or with non-templated spacer `cand_seq[left_s1..right_s0]`)
        dup_r0 = right_r0;
        dup_r1 = left_r1;
        if right_s0 > left_s1 {
            ins_seq = cand_seq[left_s1..right_s0].to_vec();
        }
    } else {
        // Require a templated segment in a complex delins to characterize it as an ITD.
        let mid = &cand_seq[left_s1..right_s0];
        let mid_aln = align_single_strand(wt_index, mid, true, 16, 10)?;
        if mid_aln.exact_matches * 2 < mid.len()
            || (mid_aln.ref_start >= left_r1 && mid_aln.ref_end <= right_r0)
        {
            return None;
        }
        dup_r0 = mid_aln.ref_start;
        dup_r1 = mid_aln.ref_end;
        ins_seq.extend_from_slice(&mid[..mid_aln.query_start]);
        ins_seq.extend_from_slice(&mid[mid_aln.query_end..]);
    }

    if dup_r1 <= dup_r0 {
        return None;
    }
    let dup_len = dup_r1 - dup_r0;

    // Count mismatches inside arms
    let left_mm_count = left_mm.iter().filter(|&&p| p < left_len).count();
    let right_mm_count = right_mm_off.len();
    let mismatches = left_mm_count + right_mm_count;

    // Build clean structural mutant sequence for VCF allele derivation:
    // Revert isolated distant flanking SNPs (> 5 bp outside the ITD event bounds) so they don't bloat REF/ALT.
    let event_r_lo = dup_r0.min(left_r1).min(right_r0).saturating_sub(5);
    let event_r_hi = (dup_r1.max(left_r1).max(right_r0) + 5).min(wt_len);

    let mut clean_cand = cand_seq.clone();
    for &p in &left_mm {
        if p < left_s1 && p < event_r_lo {
            clean_cand[p] = wt_ref[p];
        }
    }
    for &off in &right_mm_off {
        let r_pos = wt_len - 1 - off;
        let c_pos = cand_len - 1 - off;
        if c_pos >= right_s0 && r_pos > event_r_hi {
            clean_cand[c_pos] = wt_ref[r_pos];
        }
    }

    // Derive minimal parsimonious VCF REF/ALT by comparing wt_ref vs clean_cand
    let mut prefix = 0usize;
    let min_l = wt_len.min(cand_len);
    while prefix < min_l && wt_ref[prefix] == clean_cand[prefix] {
        prefix += 1;
    }

    let mut suffix = 0usize;
    let max_suffix = (wt_len - prefix).min(cand_len - prefix);
    while suffix < max_suffix && wt_ref[wt_len - 1 - suffix] == clean_cand[cand_len - 1 - suffix] {
        suffix += 1;
    }

    let raw_pos = prefix;
    let raw_ref = wt_ref[prefix..wt_len - suffix].to_vec();
    let raw_alt = clean_cand[prefix..cand_len - suffix].to_vec();

    // Anchor on preceding base if pure insertion/deletion
    let (anchor_pos, anchor_ref, anchor_alt) =
        if (raw_ref.is_empty() || raw_alt.is_empty()) && raw_pos > 0 {
            let p = raw_pos - 1;
            let mut r = vec![wt_ref[p]];
            r.extend_from_slice(&raw_ref);
            let mut a = vec![wt_ref[p]];
            a.extend_from_slice(&raw_alt);
            (p, r, a)
        } else {
            (raw_pos, raw_ref, raw_alt)
        };

    let (vcf_local_pos, vcf_ref, vcf_alt) =
        left_normalize_variant(wt_ref, anchor_pos, anchor_ref, anchor_alt);

    // Ensure mutbk_s0 <= mutbk_s1
    if mutbk_s1 < mutbk_s0 {
        mutbk_s1 = mutbk_s0;
    }

    Some(CharacterizedItd {
        id,
        candidate_seq: cand_seq.clone(),
        initial_support: cluster.total_support,
        itd_size,
        dup_r0,
        dup_r1,
        dup_len,
        ins_seq,
        mutbk_s0,
        mutbk_s1,
        mismatches,
        vcf_local_pos,
        vcf_ref,
        vcf_alt,
    })
}

/// Second-stage structural deduplication of characterized candidates (removes shadow candidates
/// that collapse to the same normalized VCF variant or nearly identical allele within edit dist).
pub fn deduplicate_candidates(
    mut cands: Vec<CharacterizedItd>,
    cluster_edit_dist: usize,
) -> Vec<CharacterizedItd> {
    cands.sort_by(|a, b| {
        b.initial_support
            .cmp(&a.initial_support)
            .then_with(|| a.mismatches.cmp(&b.mismatches))
            .then_with(|| a.itd_size.cmp(&b.itd_size))
    });

    let mut kept: Vec<CharacterizedItd> = Vec::new();
    for cand in cands {
        let is_dup = kept.iter().any(|existing| {
            // Exact same normalized VCF allele
            if existing.vcf_local_pos == cand.vcf_local_pos
                && existing.vcf_ref == cand.vcf_ref
                && existing.vcf_alt == cand.vcf_alt
            {
                return true;
            }
            // Or nearly identical full allele sequence
            banded_edit_distance(
                &existing.candidate_seq,
                &cand.candidate_seq,
                cluster_edit_dist,
            )
            .is_some()
        });
        if !is_dup {
            kept.push(cand);
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_banded_edit_distance() {
        assert_eq!(banded_edit_distance(b"ACGTACGT", b"ACGTACGT", 5), Some(0));
        assert_eq!(banded_edit_distance(b"ACGTACGT", b"ACGTTCGT", 5), Some(1));
        assert_eq!(banded_edit_distance(b"ACGTACGT", b"ACGTAACGT", 5), Some(1));
        assert_eq!(banded_edit_distance(b"ACGTACGT", b"TTTTTTTT", 3), None);
    }

    #[test]
    fn test_characterize_pure_and_complex_itd_with_flanking_snp() {
        let wt = concat!(
            "TGAATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCT",
            "ACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAG",
            "AGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAAT",
            "CTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTT",
            "TTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGG"
        )
        .as_bytes()
        .to_vec();
        let idx = LocalIndex::new(wt.clone());

        // 30 bp dup of wt[120..150] + 4 bp non-templated insertion "CCCT" + distant flanking SNP at pos 30
        let mut cand = Vec::new();
        cand.extend_from_slice(&wt[..150]);
        cand.extend_from_slice(b"CCCT");
        cand.extend_from_slice(&wt[120..]);
        cand[30] = if cand[30] == b'A' { b'G' } else { b'A' }; // germline SNP far away

        let cl = CandidateCluster {
            consensus_seq: cand,
            total_support: 10,
            itd_size: 34,
        };

        let ch = characterize_candidate(0, &wt, &idx, &cl).unwrap();
        assert_eq!(ch.itd_size, 34);
        assert_eq!(ch.dup_r0, 120);
        assert_eq!(ch.dup_r1, 150);
        assert_eq!(ch.ins_seq, b"CCCT");
        // Verify VCF allele reconstructs clean ITD (without bloating to pos 30!)
        assert!(ch.vcf_local_pos >= 115);
        assert_eq!(ch.vcf_alt.len() - ch.vcf_ref.len(), 34);
    }
}
