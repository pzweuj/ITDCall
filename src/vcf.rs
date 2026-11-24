use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::cli::AssayMethod;
use crate::quantify::QuantifiedItd;
use crate::target::TargetLocus;

#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: usize, // 1-based genomic POS
    pub ref_allele: String,
    pub alt_allele: String,
    pub filter: String,
    pub svlen: usize,
    pub dup_start: usize, // 1-based inclusive
    pub dup_end: usize,   // 1-based inclusive
    pub dup_len: usize,
    pub ins_len: usize,
    pub method: AssayMethod,
    pub locus_name: String,
    pub af: f64,
    pub ar: f64,
    pub dp: f64,
    pub vd: f64,
    pub coverage_radius: [usize; 6],
}

impl VcfRecord {
    pub fn from_quantified(locus: &TargetLocus, q: &QuantifiedItd, method: AssayMethod) -> Self {
        let pos = locus.start + q.cand.vcf_local_pos;
        let dup_start = locus.start + q.cand.dup_r0;
        let dup_end = if q.cand.dup_r1 > q.cand.dup_r0 {
            locus.start + q.cand.dup_r1 - 1
        } else {
            dup_start
        };
        let filter = if q.filters.is_empty() {
            "PASS".to_string()
        } else {
            q.filters.join(";")
        };

        Self {
            chrom: locus.chrom.clone(),
            pos,
            ref_allele: String::from_utf8_lossy(&q.cand.vcf_ref).into_owned(),
            alt_allele: String::from_utf8_lossy(&q.cand.vcf_alt).into_owned(),
            filter,
            svlen: q.cand.itd_size,
            dup_start,
            dup_end,
            dup_len: q.cand.dup_len,
            ins_len: q.cand.ins_seq.len(),
            method,
            locus_name: locus.name.clone(),
            af: q.af,
            ar: q.ar,
            dp: q.dp,
            vd: q.vd,
            coverage_radius: q.coverage_radius,
        }
    }
}

fn fmt_float(v: f64) -> String {
    if v.is_infinite() {
        "Inf".to_string()
    } else if v.is_nan() {
        ".".to_string()
    } else if (v - v.round()).abs() < 1e-9 && v.abs() >= 1.0 {
        format!("{:.1}", v)
    } else {
        format!("{:.4}", v)
    }
}

fn fmt_depth(v: f64) -> String {
    if (v - v.round()).abs() < 1e-9 {
        format!("{:.0}", v)
    } else {
        format!("{:.1}", v)
    }
}

pub fn write_vcf(
    output_path: &Path,
    reference_path: &Path,
    sample_name: &str,
    contigs: &[(String, usize)],
    mut records: Vec<VcfRecord>,
) -> Result<()> {
    let contig_order: HashMap<&str, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.as_str(), i))
        .collect();

    records.sort_by(|a, b| {
        let ca = contig_order
            .get(a.chrom.as_str())
            .copied()
            .unwrap_or(usize::MAX);
        let cb = contig_order
            .get(b.chrom.as_str())
            .copied()
            .unwrap_or(usize::MAX);
        ca.cmp(&cb)
            .then_with(|| a.chrom.cmp(&b.chrom))
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            // Prefer PASS, then higher support and depth, for overlapping target calls.
            .then_with(|| (b.filter == "PASS").cmp(&(a.filter == "PASS")))
            .then_with(|| b.vd.total_cmp(&a.vd))
            .then_with(|| b.dp.total_cmp(&a.dp))
            .then_with(|| a.locus_name.cmp(&b.locus_name))
    });
    records.dedup_by(|a, b| {
        a.chrom == b.chrom
            && a.pos == b.pos
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });

    if output_path == Path::new("-") {
        let stdout = io::stdout();
        let mut writer = BufWriter::new(stdout.lock());
        write_vcf_to(&mut writer, reference_path, sample_name, contigs, &records)?;
        writer.flush()?;
    } else {
        let file = File::create(output_path)
            .with_context(|| format!("Failed to create output VCF {}", output_path.display()))?;
        let mut writer = BufWriter::new(file);
        write_vcf_to(&mut writer, reference_path, sample_name, contigs, &records)?;
        writer.flush()?;
    }

    Ok(())
}

fn write_vcf_to<W: Write>(
    w: &mut W,
    reference_path: &Path,
    sample_name: &str,
    contigs: &[(String, usize)],
    records: &[VcfRecord],
) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=itdcall v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference={}", reference_path.display())?;

    for (name, len) in contigs {
        writeln!(w, "##contig=<ID={name},length={len}>")?;
    }

    writeln!(w, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    writeln!(
        w,
        "##FILTER=<ID=LOW_SUPPORT,Description=\"Effective supporting mutant depth below --min-reads\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=LOW_VAF,Description=\"Variant allele frequency below --min-vaf\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=FRAMESHIFT,Description=\"ITD net insertion length is not a multiple of 3 (--in-frame-only)\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=NOISY_COVERAGE,Description=\"Large ITD (>=200 bp) with >=50% uncovered duplicated sequence\">"
    )?;

    writeln!(
        w,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Net insertion length of the ITD (bp)\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DUPSTART,Number=1,Type=Integer,Description=\"1-based genomic start of duplicated reference segment\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DUPEND,Number=1,Type=Integer,Description=\"1-based genomic end of duplicated reference segment\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DUPLEN,Number=1,Type=Integer,Description=\"Length of duplicated reference segment (bp)\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of non-templated inserted sequence (bp)\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=METHOD,Number=1,Type=String,Description=\"Quantification methodology (hc or amplicon)\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=LOCUS,Number=1,Type=String,Description=\"Target locus identifier\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Effective total depth at ITD locus\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=VD,Number=1,Type=Float,Description=\"Effective supporting variant depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=AR,Number=1,Type=Float,Description=\"Allelic ratio (Mut/WT)\">"
    )?;

    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Effective total depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=VD,Number=1,Type=Float,Description=\"Effective variant depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=AR,Number=1,Type=Float,Description=\"Allelic ratio (Mut/WT)\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=CR,Number=6,Type=Integer,Description=\"Coverage radius histogram (0,1-5,6-10,11-25,26-50,>50)\">"
    )?;

    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )?;

    for r in records {
        let gt = if r.af >= 0.85 { "1/1" } else { "0/1" };
        let dp_s = fmt_depth(r.dp);
        let vd_s = fmt_depth(r.vd);
        let af_s = fmt_float(r.af);
        let ar_s = fmt_float(r.ar);
        let cr_s = format!(
            "{},{},{},{},{},{}",
            r.coverage_radius[0],
            r.coverage_radius[1],
            r.coverage_radius[2],
            r.coverage_radius[3],
            r.coverage_radius[4],
            r.coverage_radius[5]
        );

        let info = format!(
            "SVTYPE=ITD;SVLEN={};DUPSTART={};DUPEND={};DUPLEN={};INSLEN={};METHOD={};LOCUS={};DP={};VD={};AF={};AR={}",
            r.svlen,
            r.dup_start,
            r.dup_end,
            r.dup_len,
            r.ins_len,
            r.method.as_str(),
            r.locus_name,
            dp_s,
            vd_s,
            af_s,
            ar_s
        );

        writeln!(
            w,
            "{}\t{}\t.\t{}\t{}\t.\t{}\t{}\tGT:DP:VD:AF:AR:CR\t{}:{}:{}:{}:{}:{}",
            r.chrom,
            r.pos,
            r.ref_allele,
            r.alt_allele,
            r.filter,
            info,
            gt,
            dp_s,
            vd_s,
            af_s,
            ar_s,
            cr_s
        )?;
    }

    Ok(())
}
