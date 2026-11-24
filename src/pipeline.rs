use std::collections::HashMap;

use anyhow::Result;
use rayon::prelude::*;

use crate::align::{LocalIndex, OrientedAlignment, align_read_both_strands, anchor_softclip};
use crate::bam_io::{ReadFragment, extract_reads_for_locus, read_bam_metadata};
use crate::cli::{AssayMethod, Cli, ReadMode};
use crate::cluster::{characterize_candidate, cluster_extended_reads, deduplicate_candidates};
use crate::quantify::{
    QuantParams, WtFragmentSpan, assign_mutant_fragments, check_stringent_wt, quantify_candidates,
};
use crate::target::{ReferenceGenome, TargetLocus, resolve_target_loci};
use crate::vcf::{VcfRecord, write_vcf};

#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub method: AssayMethod,
    pub read_mode: ReadMode,
    pub exclude_duplicates: bool,
    pub trim_adapters: bool,
    pub min_reads: f64,
    pub min_vaf: f64,
    pub min_itd_len: usize,
    pub max_itd_len: usize,
    pub min_clip_len: usize,
    pub buffer: usize,
    pub max_edit_dist: usize,
    pub max_clip: usize,
    pub max_ins: usize,
    pub max_del: usize,
    pub cluster_edit_dist: usize,
    pub in_frame_only: bool,
    pub include_filtered: bool,
}

impl From<&Cli> for PipelineConfig {
    fn from(cli: &Cli) -> Self {
        Self {
            method: cli.method,
            read_mode: cli.read_mode,
            exclude_duplicates: cli.exclude_duplicates,
            trim_adapters: !cli.no_trim_adapters,
            min_reads: cli.min_reads,
            min_vaf: cli.min_vaf,
            min_itd_len: cli.min_itd_len,
            max_itd_len: cli.max_itd_len,
            min_clip_len: cli.min_clip_len,
            buffer: cli.buffer,
            max_edit_dist: cli.max_edit_dist,
            max_clip: cli.max_clip,
            max_ins: cli.max_ins,
            max_del: cli.max_del,
            cluster_edit_dist: cli.cluster_edit_dist,
            in_frame_only: cli.in_frame_only,
            include_filtered: cli.include_filtered,
        }
    }
}

/// Try to extract a full-length extended candidate allele from a single read alignment to WT ref.
fn maybe_extend_read(
    wt_index: &LocalIndex,
    oa: &OrientedAlignment,
    cfg: &PipelineConfig,
) -> Option<Vec<u8>> {
    let aln = &oa.aln;
    let ins_len = aln.cigar.total_ins();
    let del_len = aln.cigar.total_del();
    let max_clip = aln.cigar.max_clip();
    let wt_len = wt_index.seq.len();

    // Path 1: Direct inline net insertion >= min_itd_len (short/medium ITDs bridged by local SW)
    let coords = if ins_len > del_len
        && (ins_len - del_len) >= cfg.min_itd_len
        && max_clip <= cfg.max_clip
    {
        Some((aln.ref_start, aln.ref_end, aln.query_start, aln.query_end))
    } else if max_clip >= cfg.min_clip_len {
        // Path 2: Split-read soft-clip anchoring (medium/long ITDs, including ITDs > read length)
        anchor_softclip(
            wt_index,
            &oa.oriented_seq,
            aln,
            cfg.min_clip_len,
            cfg.min_itd_len,
            cfg.max_itd_len,
        )
    } else {
        None
    }?;

    let (lloc, rloc, qs, qe) = coords;
    if lloc > wt_len || rloc > wt_len || qe <= qs || qe > oa.oriented_seq.len() {
        return None;
    }

    let mut ext_seq = Vec::with_capacity(lloc + (qe - qs) + (wt_len - rloc));
    ext_seq.extend_from_slice(&wt_index.seq[..lloc]);
    ext_seq.extend_from_slice(&oa.oriented_seq[qs..qe]);
    ext_seq.extend_from_slice(&wt_index.seq[rloc..]);

    let itd_size = ext_seq.len().saturating_sub(wt_len);
    if itd_size < cfg.min_itd_len || itd_size > cfg.max_itd_len {
        return None;
    }
    Some(ext_seq)
}

enum FragmentOutcome {
    Wt(WtFragmentSpan),
    NonWt {
        frag_idx: usize,
        extensions: Vec<Vec<u8>>,
    },
    UnmappedToLocus,
}

/// Analyze a single target locus given extracted read fragments.
pub fn analyze_locus_fragments(
    locus: &TargetLocus,
    fragments: &[ReadFragment],
    cfg: &PipelineConfig,
) -> Vec<VcfRecord> {
    if fragments.is_empty() {
        return Vec::new();
    }

    let wt_index = LocalIndex::new(locus.seq.clone());
    let quant_params = QuantParams {
        method: cfg.method,
        buffer: cfg.buffer,
        max_ins: cfg.max_ins,
        max_del: cfg.max_del,
        max_clip: cfg.max_clip,
        max_edit_dist: cfg.max_edit_dist,
        min_reads: cfg.min_reads,
        min_vaf: cfg.min_vaf,
        in_frame_only: cfg.in_frame_only,
    };

    let min_wt_score = 15i32;
    let band_half = 48usize;

    // Step 1: Parallel local alignment of all reads to WT target sequence
    let outcomes: Vec<FragmentOutcome> = fragments
        .par_iter()
        .enumerate()
        .map(|(frag_idx, frag)| match frag {
            ReadFragment::Paired { r1, r2, .. } => {
                let a1 = align_read_both_strands(&wt_index, &r1.seq, band_half, min_wt_score);
                let a2 = align_read_both_strands(&wt_index, &r2.seq, band_half, min_wt_score);

                if a1.is_none() && a2.is_none() {
                    return FragmentOutcome::UnmappedToLocus;
                }

                if let (Some(o1), Some(o2)) = (&a1, &a2)
                    && let Some(wt_span) = check_stringent_wt(
                        &o1.aln,
                        o1.is_reverse,
                        Some((&o2.aln, o2.is_reverse)),
                        &quant_params,
                    )
                {
                    return FragmentOutcome::Wt(wt_span);
                }

                let mut extensions = Vec::new();
                if let Some(o1) = &a1
                    && let Some(ext) = maybe_extend_read(&wt_index, o1, cfg)
                {
                    extensions.push(ext);
                }
                if let Some(o2) = &a2
                    && let Some(ext) = maybe_extend_read(&wt_index, o2, cfg)
                {
                    extensions.push(ext);
                }
                FragmentOutcome::NonWt {
                    frag_idx,
                    extensions,
                }
            }
            ReadFragment::Single { r1, .. } => {
                let Some(o1) = align_read_both_strands(&wt_index, &r1.seq, band_half, min_wt_score)
                else {
                    return FragmentOutcome::UnmappedToLocus;
                };
                if let Some(wt_span) =
                    check_stringent_wt(&o1.aln, o1.is_reverse, None, &quant_params)
                {
                    return FragmentOutcome::Wt(wt_span);
                }
                let mut extensions = Vec::new();
                if let Some(ext) = maybe_extend_read(&wt_index, &o1, cfg) {
                    extensions.push(ext);
                }
                FragmentOutcome::NonWt {
                    frag_idx,
                    extensions,
                }
            }
        })
        .collect();

    let mut wt_spans: Vec<WtFragmentSpan> = Vec::new();
    let mut non_wt_refs: Vec<&ReadFragment> = Vec::new();
    let mut ext_counts: HashMap<Vec<u8>, usize> = HashMap::new();

    for out in outcomes {
        match out {
            FragmentOutcome::Wt(span) => wt_spans.push(span),
            FragmentOutcome::NonWt {
                frag_idx,
                extensions,
            } => {
                non_wt_refs.push(&fragments[frag_idx]);
                for ext in extensions {
                    *ext_counts.entry(ext).or_insert(0) += 1;
                }
            }
            FragmentOutcome::UnmappedToLocus => {}
        }
    }

    if ext_counts.is_empty() {
        return Vec::new();
    }

    // Step 2: Cluster extended candidate sequences (pure-Rust greedy edit-distance clustering)
    let clusters = cluster_extended_reads(&ext_counts, locus.seq.len(), cfg.cluster_edit_dist);
    if clusters.is_empty() {
        return Vec::new();
    }

    // Step 3: Characterize ITD structure & deduplicate equivalent alleles
    let mut characterized = Vec::new();
    // Short matches at an insertion boundary can come from microhomology alone.
    let min_dup_len = cfg.min_itd_len.clamp(1, 3);
    for (idx, cl) in clusters.iter().enumerate() {
        if let Some(ch) = characterize_candidate(idx, &locus.seq, &wt_index, cl)
            && ch.dup_len >= min_dup_len
        {
            characterized.push(ch);
        }
    }
    let candidates = deduplicate_candidates(characterized, cfg.cluster_edit_dist);
    if candidates.is_empty() {
        return Vec::new();
    }

    // Step 4: Realign non-WT fragments to candidate ITD alleles
    let mut_alns = assign_mutant_fragments(&non_wt_refs, &candidates, &quant_params);

    // Step 5: Quantify VD / DP / AF / AR according to assay methodology (HC vs Amplicon)
    let quantified = quantify_candidates(candidates, &wt_spans, &mut_alns, &quant_params);

    // Step 6: Convert to VCF records
    quantified
        .into_iter()
        .filter(|q| cfg.include_filtered || q.filters.is_empty())
        .map(|q| VcfRecord::from_quantified(locus, &q, cfg.method))
        .collect()
}

/// Run the full CLI pipeline from BAM + reference FASTA to output VCF.
pub fn run_cli(cli: &Cli) -> Result<()> {
    if cli.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global();
    }

    let reference = ReferenceGenome::open(&cli.reference)?;
    let bam_meta = read_bam_metadata(&cli.bam, cli.sample.as_deref())?;
    let loci = resolve_target_loci(
        &reference,
        cli.genome,
        cli.target_fasta.as_deref(),
        &cli.regions,
        cli.bed.as_deref(),
    )?;

    let cfg = PipelineConfig::from(cli);
    let mut all_records = Vec::new();

    for locus in &loci {
        eprintln!(
            "[itdcall] Processing locus '{}' ({}:{}-{}, {} bp, method={})",
            locus.name,
            locus.chrom,
            locus.start,
            locus.end,
            locus.seq.len(),
            cfg.method.as_str()
        );
        let fragments = extract_reads_for_locus(
            &cli.bam,
            locus,
            cfg.read_mode,
            cfg.exclude_duplicates,
            cfg.trim_adapters,
        )?;
        eprintln!(
            "[itdcall]   Extracted {} read fragments from BAM",
            fragments.len()
        );
        let mut recs = analyze_locus_fragments(locus, &fragments, &cfg);
        eprintln!("[itdcall]   Called {} ITD variant(s)", recs.len());
        all_records.append(&mut recs);
    }

    // Prefer reference FASTA contigs for VCF header (falling back to BAM contigs)
    let header_contigs: Vec<(String, usize)> = if !reference.contigs.is_empty() {
        reference
            .contigs
            .iter()
            .map(|c| (c.name.clone(), c.length))
            .collect()
    } else {
        bam_meta.contigs
    };

    write_vcf(
        &cli.output,
        &cli.reference,
        &bam_meta.sample_name,
        &header_contigs,
        all_records,
    )?;

    Ok(())
}
