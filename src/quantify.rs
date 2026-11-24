use rayon::prelude::*;

use crate::align::{LocalAlignment, LocalIndex, align_read_both_strands};
use crate::bam_io::ReadFragment;
use crate::cli::AssayMethod;
use crate::cluster::CharacterizedItd;

#[derive(Debug, Clone)]
pub struct WtFragmentSpan {
    /// Aligned 0-based half-open intervals `(ref_start, ref_end)` on WT target sequence (1 or 2 mates).
    pub spans: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
pub struct MutFragmentAln {
    pub cand_idx: usize,
    /// Aligned 0-based half-open intervals `(ref_start, ref_end)` on `candidate_seq` (1 or 2 mates).
    pub spans: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
pub struct QuantParams {
    pub method: AssayMethod,
    pub buffer: usize,
    pub max_ins: usize,
    pub max_del: usize,
    pub max_clip: usize,
    pub max_edit_dist: usize,
    pub min_reads: f64,
    pub min_vaf: f64,
    pub in_frame_only: bool,
}

#[derive(Debug, Clone)]
pub struct QuantifiedItd {
    pub cand: CharacterizedItd,
    pub vd: f64,
    pub dp: f64,
    pub af: f64,
    pub ar: f64,
    pub mut_fragments: usize,
    pub coverage_radius: [usize; 6], // 0, 1-5, 6-10, 11-25, 26-50, >50
    pub filters: Vec<String>,
}

#[inline]
fn spans_boundary(ref_start: usize, ref_end: usize, bk: usize, buffer: usize) -> bool {
    ref_start + buffer <= bk && bk + buffer <= ref_end
}

fn any_span_crosses(spans: &[(usize, usize)], bk: usize, buffer: usize) -> bool {
    spans.iter().any(|&(s, e)| spans_boundary(s, e, bk, buffer))
}

/// Realign non-WT fragments against each candidate ITD allele (ordered by initial cluster support)
/// and greedily assign each fragment to at most one ITD allele.
pub fn assign_mutant_fragments(
    non_wt_fragments: &[&ReadFragment],
    candidates: &[CharacterizedItd],
    params: &QuantParams,
) -> Vec<MutFragmentAln> {
    if non_wt_fragments.is_empty() || candidates.is_empty() {
        return Vec::new();
    }

    let mut assigned_flags = vec![false; non_wt_fragments.len()];
    let mut all_mut: Vec<MutFragmentAln> = Vec::new();
    let min_score = 20i32;

    for (cand_idx, cand) in candidates.iter().enumerate() {
        let cand_index = LocalIndex::new(cand.candidate_seq.clone());

        // Evaluate unassigned fragments in parallel against this candidate ITD
        let hits: Vec<(usize, Vec<(usize, usize)>)> = non_wt_fragments
            .par_iter()
            .enumerate()
            .filter_map(|(frag_idx, frag)| {
                if assigned_flags[frag_idx] {
                    return None;
                }
                match frag {
                    ReadFragment::Paired { r1, r2, .. } => {
                        let a1 = align_read_both_strands(&cand_index, &r1.seq, 16, min_score)?;
                        if !a1.aln.is_stringent(
                            params.max_ins,
                            params.max_del,
                            params.max_clip,
                            params.max_edit_dist,
                        ) {
                            return None;
                        }
                        let a2 = align_read_both_strands(&cand_index, &r2.seq, 16, min_score)?;
                        if !a2.aln.is_stringent(
                            params.max_ins,
                            params.max_del,
                            params.max_clip,
                            params.max_edit_dist,
                        ) {
                            return None;
                        }
                        // Require proper opposite-strand pair orientation
                        if a1.is_reverse == a2.is_reverse {
                            return None;
                        }
                        let spans = vec![
                            (a1.aln.ref_start, a1.aln.ref_end),
                            (a2.aln.ref_start, a2.aln.ref_end),
                        ];
                        let crosses_mj = any_span_crosses(&spans, cand.mutbk_s0, params.buffer)
                            || any_span_crosses(&spans, cand.mutbk_s1, params.buffer);
                        crosses_mj.then_some((frag_idx, spans))
                    }
                    ReadFragment::Single { r1, .. } => {
                        let a1 = align_read_both_strands(&cand_index, &r1.seq, 16, min_score)?;
                        if !a1.aln.is_stringent(
                            params.max_ins,
                            params.max_del,
                            params.max_clip,
                            params.max_edit_dist,
                        ) {
                            return None;
                        }
                        let spans = vec![(a1.aln.ref_start, a1.aln.ref_end)];
                        let crosses_mj = any_span_crosses(&spans, cand.mutbk_s0, params.buffer)
                            || any_span_crosses(&spans, cand.mutbk_s1, params.buffer);
                        crosses_mj.then_some((frag_idx, spans))
                    }
                }
            })
            .collect();

        for (frag_idx, spans) in hits {
            assigned_flags[frag_idx] = true;
            all_mut.push(MutFragmentAln { cand_idx, spans });
        }
    }

    all_mut
}

/// Project a span `(s, e)` on candidate `cand` back onto approximate WT target coordinates
/// for cross-clone depth accounting in multi-ITD samples.
fn project_mut_span_to_wt(s: usize, e: usize, cand: &CharacterizedItd) -> (usize, usize) {
    if e <= cand.mutbk_s0 {
        (s, e)
    } else if s >= cand.mutbk_s1 {
        (
            s.saturating_sub(cand.itd_size),
            e.saturating_sub(cand.itd_size),
        )
    } else {
        let wt_s = s.min(cand.dup_r0);
        let wt_e = e.saturating_sub(cand.itd_size).max(cand.dup_r1);
        (wt_s, wt_e)
    }
}

/// Compute VD, DP, AF, AR, coverage radius, and FILTER flags for every candidate ITD.
pub fn quantify_candidates(
    candidates: Vec<CharacterizedItd>,
    wt_fragments: &[WtFragmentSpan],
    mut_fragments: &[MutFragmentAln],
    params: &QuantParams,
) -> Vec<QuantifiedItd> {
    let n_cands = candidates.len();
    if n_cands == 0 {
        return Vec::new();
    }

    // Per-candidate coverage counters
    let mut mut_left = vec![0usize; n_cands];
    let mut mut_right = vec![0usize; n_cands];
    let mut mut_any = vec![0usize; n_cands];
    let mut mut_wt_left = vec![0usize; n_cands];
    let mut mut_wt_right = vec![0usize; n_cands];
    let mut other_mut_left = vec![0usize; n_cands];
    let mut other_mut_right = vec![0usize; n_cands];
    let mut other_mut_any = vec![0usize; n_cands];
    let mut wt_left = vec![0usize; n_cands];
    let mut wt_right = vec![0usize; n_cands];
    let mut wt_any = vec![0usize; n_cands];

    // Coverage radius arrays across duplicated + inserted region on each candidate allele
    let mut cov_radii: Vec<Vec<usize>> = candidates
        .iter()
        .map(|c| {
            let r_start = c.dup_r0;
            let r_end = (c.dup_r1 + c.itd_size).min(c.candidate_seq.len());
            vec![0usize; r_end.saturating_sub(r_start)]
        })
        .collect();

    // 1) Accumulate mutant fragment statistics
    for mf in mut_fragments {
        let cidx = mf.cand_idx;
        let cand = &candidates[cidx];

        let ml = any_span_crosses(&mf.spans, cand.mutbk_s0, params.buffer);
        let mr = any_span_crosses(&mf.spans, cand.mutbk_s1, params.buffer);
        if ml {
            mut_left[cidx] += 1;
        }
        if mr {
            mut_right[cidx] += 1;
        }
        if ml || mr {
            mut_any[cidx] += 1;
        }

        // Outer WT-like boundaries on the mutant allele
        let mwl = any_span_crosses(&mf.spans, cand.dup_r0, params.buffer);
        let mwr = any_span_crosses(&mf.spans, cand.dup_r1 + cand.itd_size, params.buffer);
        if mwl {
            mut_wt_left[cidx] += 1;
        }
        if mwr {
            mut_wt_right[cidx] += 1;
        }

        // Update coverage radius for mutant reads supporting this candidate
        if ml || mr {
            let r_start = cand.dup_r0;
            let r_end = r_start + cov_radii[cidx].len();
            for &(s, e) in &mf.spans {
                let lo = s.max(r_start);
                let hi = e.min(r_end);
                for pos in lo..hi {
                    let rad = (pos - s + 1).min(e - pos);
                    let slot = &mut cov_radii[cidx][pos - r_start];
                    if rad > *slot {
                        *slot = rad;
                    }
                }
            }
        }

        // Cross-clone depth accounting: check if this mutant fragment spans other ITDs' WT boundaries
        if n_cands > 1 && (ml || mr) {
            let proj_spans: Vec<(usize, usize)> = mf
                .spans
                .iter()
                .map(|&(s, e)| project_mut_span_to_wt(s, e, cand))
                .collect();
            for (other_idx, other) in candidates.iter().enumerate() {
                if other_idx == cidx {
                    continue;
                }
                let ol = any_span_crosses(&proj_spans, other.dup_r0, params.buffer);
                let or = any_span_crosses(&proj_spans, other.dup_r1, params.buffer);
                if ol {
                    other_mut_left[other_idx] += 1;
                }
                if or {
                    other_mut_right[other_idx] += 1;
                }
                if ol || or {
                    other_mut_any[other_idx] += 1;
                }
            }
        }
    }

    // 2) Accumulate WT fragment statistics
    for wf in wt_fragments {
        for (cidx, cand) in candidates.iter().enumerate() {
            let wl = any_span_crosses(&wf.spans, cand.dup_r0, params.buffer);
            let wr = any_span_crosses(&wf.spans, cand.dup_r1, params.buffer);
            if wl {
                wt_left[cidx] += 1;
            }
            if wr {
                wt_right[cidx] += 1;
            }
            if wl || wr {
                wt_any[cidx] += 1;
            }
        }
    }

    // 3) Compute VD, DP, AF per methodology (HC vs Amplicon)
    let mut results: Vec<QuantifiedItd> = Vec::with_capacity(n_cands);
    for (cidx, cand) in candidates.into_iter().enumerate() {
        if mut_any[cidx] == 0 {
            continue;
        }

        let (vd, raw_dp, wt_depth) = match params.method {
            AssayMethod::Hc => {
                let vd = (mut_left[cidx] as f64 + mut_right[cidx] as f64) / 2.0;
                let dp = (wt_left[cidx] as f64
                    + wt_right[cidx] as f64
                    + mut_wt_left[cidx] as f64
                    + mut_wt_right[cidx] as f64
                    + other_mut_left[cidx] as f64
                    + other_mut_right[cidx] as f64)
                    / 2.0;
                let wt_depth = (wt_left[cidx] as f64 + wt_right[cidx] as f64) / 2.0;
                (vd, dp, wt_depth)
            }
            AssayMethod::Amplicon => {
                let vd = mut_any[cidx] as f64;
                let dp = (wt_any[cidx] + mut_any[cidx] + other_mut_any[cidx]) as f64;
                (vd, dp, wt_any[cidx] as f64)
            }
        };

        if vd <= 0.0 {
            continue;
        }
        let dp = raw_dp.max(vd);
        let af = if dp > 0.0 { vd / dp } else { 0.0 };

        // Bin coverage radius histogram
        let mut hist = [0usize; 6];
        for &r in &cov_radii[cidx] {
            if r == 0 {
                hist[0] += 1;
            } else if r <= 5 {
                hist[1] += 1;
            } else if r <= 10 {
                hist[2] += 1;
            } else if r <= 25 {
                hist[3] += 1;
            } else if r <= 50 {
                hist[4] += 1;
            } else {
                hist[5] += 1;
            }
        }

        let mut filters = Vec::new();
        if vd < params.min_reads {
            filters.push("LOW_SUPPORT".to_string());
        }
        if af < params.min_vaf {
            filters.push("LOW_VAF".to_string());
        }
        if params.in_frame_only && cand.itd_size % 3 != 0 {
            filters.push("FRAMESHIFT".to_string());
        }
        let cov_len = cov_radii[cidx].len();
        if cand.itd_size >= 200 && cov_len > 0 && (hist[0] as f64) / (cov_len as f64) >= 0.5 {
            filters.push("NOISY_COVERAGE".to_string());
        }

        let ar = if wt_depth > 0.0 {
            vd / wt_depth
        } else {
            f64::INFINITY
        };

        results.push(QuantifiedItd {
            cand,
            vd,
            dp,
            af,
            ar,
            mut_fragments: mut_any[cidx],
            coverage_radius: hist,
            filters,
        });
    }

    results
}

/// Helper to build a stringent WT fragment span if a read pair / single read aligns cleanly to WT ref.
pub fn check_stringent_wt(
    a1: &LocalAlignment,
    is_rev1: bool,
    a2: Option<(&LocalAlignment, bool)>,
    params: &QuantParams,
) -> Option<WtFragmentSpan> {
    if !a1.is_stringent(
        params.max_ins,
        params.max_del,
        params.max_clip,
        params.max_edit_dist,
    ) {
        return None;
    }
    match a2 {
        Some((aln2, is_rev2)) => {
            if is_rev1 == is_rev2 {
                return None;
            }
            if !aln2.is_stringent(
                params.max_ins,
                params.max_del,
                params.max_clip,
                params.max_edit_dist,
            ) {
                return None;
            }
            Some(WtFragmentSpan {
                spans: vec![(a1.ref_start, a1.ref_end), (aln2.ref_start, aln2.ref_end)],
            })
        }
        None => Some(WtFragmentSpan {
            spans: vec![(a1.ref_start, a1.ref_end)],
        }),
    }
}
