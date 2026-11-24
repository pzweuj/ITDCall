use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow};
use noodles::bam;
use noodles::core::{Position, Region};
use noodles::sam::alignment::Record as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

use crate::cli::ReadMode;
use crate::seq::{normalize_dna, revcomp, trim_illumina_adapter};
use crate::target::TargetLocus;

#[derive(Debug, Clone)]
pub struct ReadData {
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone)]
pub enum ReadFragment {
    Paired {
        name: String,
        r1: ReadData,
        r2: ReadData,
    },
    Single {
        name: String,
        r1: ReadData,
    },
}

#[derive(Debug, Clone)]
pub struct BamMetadata {
    pub sample_name: String,
    pub contigs: Vec<(String, usize)>,
}

/// Inspect BAM header to extract contig list and default sample name (`@RG SM` or filename stem).
pub fn read_bam_metadata(bam_path: &Path, sample_override: Option<&str>) -> Result<BamMetadata> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .with_context(|| format!("Failed to open BAM {}", bam_path.display()))?;
    let header = reader.read_header()?;

    let contigs: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            (
                String::from_utf8_lossy(name.as_ref()).into_owned(),
                map.length().get(),
            )
        })
        .collect();

    let sample_name = if let Some(s) = sample_override {
        s.to_string()
    } else {
        extract_rg_sample(&header).unwrap_or_else(|| {
            bam_path
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "SAMPLE".to_string())
        })
    };

    Ok(BamMetadata {
        sample_name,
        contigs,
    })
}

fn extract_rg_sample(header: &noodles::sam::Header) -> Option<String> {
    for (_, rg) in header.read_groups() {
        if let Some(sm) = rg.other_fields().get(&rg_tag::SAMPLE) {
            let s = String::from_utf8_lossy(sm.as_ref()).trim().to_string();
            if !s.is_empty() {
                return Some(s);
            }
        }
    }
    None
}

/// Resolve BAM contig name matching `target_chrom` (handles `chr13` <-> `13`).
fn resolve_bam_contig_index(
    header: &noodles::sam::Header,
    target_chrom: &str,
) -> Option<(usize, String)> {
    let refs = header.reference_sequences();
    if let Some((idx, name, _)) = refs.get_full(target_chrom.as_bytes()) {
        return Some((idx, String::from_utf8_lossy(name.as_ref()).into_owned()));
    }
    let alt = if let Some(stripped) = target_chrom.strip_prefix("chr") {
        stripped.to_string()
    } else {
        format!("chr{target_chrom}")
    };
    refs.get_full(alt.as_bytes())
        .map(|(idx, name, _)| (idx, String::from_utf8_lossy(name.as_ref()).into_owned()))
}

fn require_bam_contig(
    header: &noodles::sam::Header,
    locus: &TargetLocus,
) -> Result<(usize, String)> {
    resolve_bam_contig_index(header, &locus.chrom).ok_or_else(|| {
        anyhow!(
            "Target contig '{}' (or chr/no-chr alias) not found in BAM header",
            locus.chrom
        )
    })
}

fn find_bai_path(bam_path: &Path) -> Option<PathBuf> {
    let p1 = bam_path.with_added_extension("bai");
    if p1.exists() {
        return Some(p1);
    }
    let p2 = bam_path.with_extension("bai");
    if p2.exists() {
        return Some(p2);
    }
    None
}

#[derive(Default)]
struct PairBucket {
    r1: Option<ReadData>,
    r2: Option<ReadData>,
    single: Option<ReadData>,
}

fn ingest_record(
    record: &bam::Record,
    exclude_duplicates: bool,
    trim_adapters: bool,
    buckets: &mut HashMap<String, PairBucket>,
) {
    let flags = record.flags();
    if flags.contains(Flags::SECONDARY)
        || flags.contains(Flags::SUPPLEMENTARY)
        || flags.contains(Flags::QC_FAIL)
    {
        return;
    }
    if exclude_duplicates && flags.contains(Flags::DUPLICATE) {
        return;
    }

    let raw_seq: Vec<u8> = record.sequence().iter().collect();
    if raw_seq.len() < 15 {
        return;
    }
    let raw_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

    // Restore original read strand before adapter trimming
    let (seq_orig, qual_orig) = if flags.is_reverse_complemented() {
        let mut q = raw_qual;
        q.reverse();
        (revcomp(&raw_seq), q)
    } else {
        (normalize_dna(&raw_seq), raw_qual)
    };

    let final_seq = if trim_adapters {
        let (tseq, _) = trim_illumina_adapter(&seq_orig, &qual_orig);
        if tseq.len() < 15 {
            return;
        }
        tseq
    } else {
        seq_orig
    };

    let qname = match record.name() {
        Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
        None => return,
    };

    let bucket = buckets.entry(qname).or_default();
    let rd = ReadData { seq: final_seq };

    if flags.is_segmented() {
        if flags.is_first_segment() {
            bucket.r1 = Some(rd);
        } else if flags.is_last_segment() {
            bucket.r2 = Some(rd);
        } else {
            bucket.single = Some(rd);
        }
    } else {
        bucket.single = Some(rd);
    }
}

/// Extract read fragments (paired-end and single-end/orphan) for a target locus from a BAM file.
pub fn extract_reads_for_locus(
    bam_path: &Path,
    locus: &TargetLocus,
    read_mode: ReadMode,
    exclude_duplicates: bool,
    trim_adapters: bool,
) -> Result<Vec<ReadFragment>> {
    let mut buckets: HashMap<String, PairBucket> = HashMap::new();
    let bai_path = match read_mode {
        ReadMode::Targeted => find_bai_path(bam_path),
        ReadMode::Loose | ReadMode::All => None,
    };
    if let Some(bai_path) = bai_path {
        let index = bam::bai::fs::read(&bai_path)
            .with_context(|| format!("Failed to read BAM index {}", bai_path.display()))?;
        let mut reader = bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(bam_path)
            .with_context(|| format!("Failed to open indexed BAM {}", bam_path.display()))?;
        let header = reader.read_header()?;

        let (_, bam_chrom) = require_bam_contig(&header, locus)?;

        let start_pos = Position::try_from(locus.start)
            .with_context(|| format!("Invalid 1-based start {}", locus.start))?;
        let end_pos = Position::try_from(locus.end)
            .with_context(|| format!("Invalid 1-based end {}", locus.end))?;
        let region = Region::new(bam_chrom.as_str(), start_pos..=end_pos);

        let query = reader.query(&header, &region)?;
        for res in query.records() {
            let record = res?;
            if !record.flags().is_unmapped() {
                ingest_record(&record, exclude_duplicates, trim_adapters, &mut buckets);
            }
        }
    } else {
        // Loose scans for all unmapped reads, including those placed beside mapped mates.
        let mut reader = bam::io::reader::Builder
            .build_from_path(bam_path)
            .with_context(|| format!("Failed to open BAM {}", bam_path.display()))?;
        let header = reader.read_header()?;
        let target_tid = match read_mode {
            ReadMode::All => None,
            ReadMode::Targeted | ReadMode::Loose => Some(require_bam_contig(&header, locus)?.0),
        };

        let mut record = bam::Record::default();
        while reader.read_record(&mut record)? != 0 {
            let keep = match read_mode {
                ReadMode::All => true,
                ReadMode::Targeted | ReadMode::Loose => {
                    let flags = record.flags();
                    if flags.is_unmapped() {
                        read_mode == ReadMode::Loose
                    } else if record.reference_sequence_id().transpose()? == target_tid {
                        let start = record.alignment_start().transpose()?;
                        let end = record.alignment_end().transpose()?;
                        matches!((start, end), (Some(s), Some(e))
                            if s.get() <= locus.end && e.get() >= locus.start)
                    } else {
                        false
                    }
                }
            };
            if keep {
                ingest_record(&record, exclude_duplicates, trim_adapters, &mut buckets);
            }
        }
    }

    let mut fragments = Vec::with_capacity(buckets.len());
    for (name, b) in buckets {
        match (b.r1, b.r2, b.single) {
            (Some(r1), Some(r2), _) => fragments.push(ReadFragment::Paired { name, r1, r2 }),
            (Some(r1), None, _) => fragments.push(ReadFragment::Single { name, r1 }),
            (None, Some(r2), _) => fragments.push(ReadFragment::Single { name, r1: r2 }),
            (None, None, Some(s)) => fragments.push(ReadFragment::Single { name, r1: s }),
            (None, None, None) => {}
        }
    }

    // Deterministic ordering by fragment name
    fragments.sort_by(|a, b| {
        let na = match a {
            ReadFragment::Paired { name, .. } | ReadFragment::Single { name, .. } => name,
        };
        let nb = match b {
            ReadFragment::Paired { name, .. } | ReadFragment::Single { name, .. } => name,
        };
        na.cmp(nb)
    });

    Ok(fragments)
}
