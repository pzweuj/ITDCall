use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use memchr::memmem;

use crate::cli::GenomeBuild;
use crate::seq::{normalize_dna, revcomp};

/// Classic FLT3 exons 14-15 locus coordinates (1-based inclusive, 1500 bp window, matching FLT3_ITD_ext).
pub const FLT3_HG19_START: usize = 28_607_438;
pub const FLT3_HG19_END: usize = 28_608_937;
pub const FLT3_HG38_START: usize = 28_033_301;
pub const FLT3_HG38_END: usize = 28_034_800;

/// 40 bp signature motif from FLT3 exon 14 (minus strand in FLT3_dna_e14e15.fa).
const FLT3_EXON14_MOTIF_MINUS: &[u8] = b"GTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCC";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

/// A resolved genomic target locus ready for BAM extraction, local alignment, and VCF emission.
#[derive(Debug, Clone)]
pub struct TargetLocus {
    pub name: String,
    /// Contig name matching the reference genome FASTA.
    pub chrom: String,
    /// 1-based inclusive genomic start coordinate on `+` strand.
    pub start: usize,
    /// 1-based inclusive genomic end coordinate on `+` strand.
    pub end: usize,
    /// Original orientation of user-provided target FASTA (always Plus when extracted from reference).
    pub source_strand: Strand,
    /// Target local reference sequence normalized to genomic `+` strand (length = end - start + 1).
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct ContigInfo {
    pub name: String,
    pub length: usize,
    offset: u64,
    line_bases: usize,
    line_width: usize,
}

/// Reference genome FASTA accessor supporting `.fai` random access and fallback streaming scan.
pub struct ReferenceGenome {
    path: PathBuf,
    pub contigs: Vec<ContigInfo>,
    contig_map: HashMap<String, usize>,
    has_fai: bool,
    /// Fallback in-memory sequences only populated if `.fai` is absent.
    in_memory: Option<HashMap<String, Vec<u8>>>,
}

impl ReferenceGenome {
    pub fn open(path: &Path) -> Result<Self> {
        let fai_path = PathBuf::from(format!("{}.fai", path.display()));
        if fai_path.exists() {
            let (contigs, contig_map) = parse_fai(&fai_path)?;
            Ok(Self {
                path: path.to_path_buf(),
                contigs,
                contig_map,
                has_fai: true,
                in_memory: None,
            })
        } else {
            // Fallback: parse FASTA directly into memory (fine for targeted/small refs or missing .fai)
            let records = read_fasta_records(path)?;
            let mut contigs = Vec::with_capacity(records.len());
            let mut contig_map = HashMap::with_capacity(records.len());
            let mut mem = HashMap::with_capacity(records.len());
            for (idx, (header, seq)) in records.into_iter().enumerate() {
                let name = header
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| anyhow!("Missing contig identifier in {}", path.display()))?
                    .to_string();
                if contig_map.insert(name.clone(), idx).is_some() {
                    bail!("Duplicate contig '{name}' in {}", path.display());
                }
                contigs.push(ContigInfo {
                    name: name.clone(),
                    length: seq.len(),
                    offset: 0,
                    line_bases: 0,
                    line_width: 0,
                });
                mem.insert(name, seq);
            }
            Ok(Self {
                path: path.to_path_buf(),
                contigs,
                contig_map,
                has_fai: false,
                in_memory: Some(mem),
            })
        }
    }

    /// Resolve contig name handling `chr13` <-> `13` aliases.
    pub fn resolve_contig<'a>(&'a self, name: &str) -> Option<&'a ContigInfo> {
        if let Some(&idx) = self.contig_map.get(name) {
            return Some(&self.contigs[idx]);
        }
        let alt = if let Some(stripped) = name.strip_prefix("chr") {
            stripped.to_string()
        } else {
            format!("chr{name}")
        };
        self.contig_map.get(&alt).map(|&idx| &self.contigs[idx])
    }

    /// Fetch 1-based inclusive genomic interval `[start, end]` from `chrom` (uppercase ACGTN).
    pub fn fetch(&self, chrom: &str, start_1based: usize, end_1based: usize) -> Result<Vec<u8>> {
        let info = self
            .resolve_contig(chrom)
            .ok_or_else(|| anyhow!("Contig '{chrom}' not found in reference FASTA"))?;
        if start_1based == 0 || end_1based < start_1based || end_1based > info.length {
            bail!(
                "Region {}:{}-{} out of bounds (contig length {})",
                info.name,
                start_1based,
                end_1based,
                info.length
            );
        }

        if let Some(mem) = &self.in_memory {
            let seq = mem
                .get(&info.name)
                .ok_or_else(|| anyhow!("Contig '{}' missing from memory store", info.name))?;
            return Ok(seq[start_1based - 1..end_1based].to_vec());
        }

        assert!(self.has_fai);
        let len = end_1based - start_1based + 1;
        let zero_start = start_1based - 1;
        let line_idx = zero_start / info.line_bases;
        let line_off = zero_start % info.line_bases;
        let byte_start =
            info.offset + (line_idx as u64) * (info.line_width as u64) + (line_off as u64);

        let lines_needed = (line_off + len).div_ceil(info.line_bases);
        let bytes_to_read = lines_needed * info.line_width;

        let mut file = File::open(&self.path)
            .with_context(|| format!("Failed to open reference FASTA {}", self.path.display()))?;
        file.seek(SeekFrom::Start(byte_start))?;
        let mut raw = vec![0u8; bytes_to_read];
        let n_read = read_up_to(&mut file, &mut raw)?;
        raw.truncate(n_read);

        let mut out = Vec::with_capacity(len);
        for &b in &raw {
            if !b.is_ascii_whitespace() {
                out.push(match b.to_ascii_uppercase() {
                    b'A' => b'A',
                    b'C' => b'C',
                    b'G' => b'G',
                    b'T' => b'T',
                    _ => b'N',
                });
                if out.len() == len {
                    break;
                }
            }
        }
        if out.len() != len {
            bail!(
                "Truncated FASTA read for {}:{}-{} (expected {} bp, got {} bp)",
                info.name,
                start_1based,
                end_1based,
                len,
                out.len()
            );
        }
        Ok(out)
    }

    /// Stream one contig sequence at a time to locate a user-provided target FASTA sequence
    /// (such as `FLT3_dna_e14e15.fa`) without loading the entire genome into RAM.
    pub fn locate_target_sequence(
        &self,
        name: &str,
        raw_seq: &[u8],
        hint_region: Option<&(String, usize, usize)>,
    ) -> Result<TargetLocus> {
        let query_plus = normalize_dna(raw_seq);
        if query_plus.len() < 50 {
            bail!(
                "Target FASTA sequence '{}' is too short ({} bp; minimum 50 bp)",
                name,
                query_plus.len()
            );
        }
        let query_minus = revcomp(&query_plus);
        let locate_in_window = |info: &ContigInfo, start: usize, window: &[u8]| {
            let (offset, strand) = match_target_in_slice(window, &query_plus, &query_minus)?;
            let exact_start = start + offset;
            Some(TargetLocus {
                name: name.to_string(),
                chrom: info.name.clone(),
                start: exact_start,
                end: exact_start + query_plus.len() - 1,
                source_strand: strand,
                seq: window[offset..offset + query_plus.len()].to_vec(),
            })
        };

        // 1) Fast path: explicit coordinate hint (from --region or FASTA header)
        if let Some((chrom, start, end)) = hint_region
            && let Some(info) = self.resolve_contig(chrom)
        {
            let ref_slice = self.fetch(&info.name, *start, *end)?;
            if let Some(locus) = locate_in_window(info, *start, &ref_slice) {
                return Ok(locus);
            }
        }

        // 2) Fast path for classic FLT3 target FASTA: check known hg38 / hg19 FLT3 windows first!
        for &(s, e) in &[
            (FLT3_HG38_START, FLT3_HG38_END),
            (FLT3_HG19_START, FLT3_HG19_END),
        ] {
            if let Some(info) = self.resolve_contig("chr13") {
                let win_s = s.saturating_sub(2000).max(1);
                let win_e = (e + 2000).min(info.length);
                if win_e >= win_s
                    && win_e - win_s + 1 >= query_plus.len()
                    && let Ok(ref_slice) = self.fetch(&info.name, win_s, win_e)
                    && let Some(locus) = locate_in_window(info, win_s, &ref_slice)
                {
                    return Ok(locus);
                }
            }
        }

        // 3) General genome-wide seed scan across all contigs in reference FASTA
        for info in &self.contigs {
            if info.length < query_plus.len() {
                continue;
            }
            let contig_seq = self.fetch(&info.name, 1, info.length)?;
            if let Some(locus) = locate_in_window(info, 1, &contig_seq) {
                return Ok(locus);
            }
        }

        bail!(
            "Could not locate target FASTA sequence '{}' ({} bp) on reference genome {}. \
             Ensure the sequence matches the reference build (or pass --region chr:start-end).",
            name,
            query_plus.len(),
            self.path.display()
        )
    }

    /// Resolve the built-in classic FLT3 locus (auto-detecting hg19 vs hg38 and chr13 vs 13).
    pub fn resolve_builtin_flt3(&self, build: GenomeBuild) -> Result<TargetLocus> {
        let motif_plus = revcomp(FLT3_EXON14_MOTIF_MINUS);
        let finder = memmem::Finder::new(&motif_plus);

        let try_coords = |start: usize, end: usize, label: &str| -> Result<Option<TargetLocus>> {
            let Some(info) = self.resolve_contig("chr13") else {
                return Ok(None);
            };
            if end > info.length {
                return Ok(None);
            }
            let seq = self.fetch(&info.name, start, end)?;
            if finder.find(&seq).is_some() {
                Ok(Some(TargetLocus {
                    name: format!("FLT3_{label}"),
                    chrom: info.name.clone(),
                    start,
                    end,
                    source_strand: Strand::Plus,
                    seq,
                }))
            } else {
                Ok(None)
            }
        };

        match build {
            GenomeBuild::Hg19 => {
                let info = self
                    .resolve_contig("chr13")
                    .ok_or_else(|| anyhow!("Contig chr13/13 not found in reference FASTA"))?;
                let seq = self.fetch(&info.name, FLT3_HG19_START, FLT3_HG19_END)?;
                Ok(TargetLocus {
                    name: "FLT3_hg19".to_string(),
                    chrom: info.name.clone(),
                    start: FLT3_HG19_START,
                    end: FLT3_HG19_END,
                    source_strand: Strand::Plus,
                    seq,
                })
            }
            GenomeBuild::Hg38 => {
                let info = self
                    .resolve_contig("chr13")
                    .ok_or_else(|| anyhow!("Contig chr13/13 not found in reference FASTA"))?;
                let seq = self.fetch(&info.name, FLT3_HG38_START, FLT3_HG38_END)?;
                Ok(TargetLocus {
                    name: "FLT3_hg38".to_string(),
                    chrom: info.name.clone(),
                    start: FLT3_HG38_START,
                    end: FLT3_HG38_END,
                    source_strand: Strand::Plus,
                    seq,
                })
            }
            GenomeBuild::Auto => {
                if let Some(locus) = try_coords(FLT3_HG38_START, FLT3_HG38_END, "hg38")? {
                    return Ok(locus);
                }
                if let Some(locus) = try_coords(FLT3_HG19_START, FLT3_HG19_END, "hg19")? {
                    return Ok(locus);
                }
                // Fallback: search reference for FLT3 exon 14 motif (handles custom/cropped reference FASTAs)
                let rc_finder = memmem::Finder::new(FLT3_EXON14_MOTIF_MINUS);
                for info in &self.contigs {
                    let seq = self.fetch(&info.name, 1, info.length)?;
                    if let Some(pos) = finder.find(&seq).or_else(|| rc_finder.find(&seq)) {
                        // Extract up to 1500 bp centered on the motif
                        let half = 750usize;
                        let start0 = pos.saturating_sub(half);
                        let end0 = (start0 + 1500).min(seq.len());
                        let start0 = end0.saturating_sub(1500);
                        return Ok(TargetLocus {
                            name: "FLT3_auto".to_string(),
                            chrom: info.name.clone(),
                            start: start0 + 1,
                            end: end0,
                            source_strand: Strand::Plus,
                            seq: seq[start0..end0].to_vec(),
                        });
                    }
                }
                bail!(
                    "Could not auto-detect built-in FLT3 locus on {}. Specify --genome hg19|hg38, --region, or --target-fasta.",
                    self.path.display()
                )
            }
        }
    }
}

/// Search `haystack` for `query_plus` or `query_minus` allowing up to 5% substitutions (user calibration SNPs).
/// Returns `Some((0based_offset, strand))` for the best match.
fn match_target_in_slice(
    haystack: &[u8],
    query_plus: &[u8],
    query_minus: &[u8],
) -> Option<(usize, Strand)> {
    let qlen = query_plus.len();
    if haystack.len() < qlen {
        return None;
    }

    let seed_len = 25usize.min(qlen);
    let seed_offsets: Vec<usize> = if qlen >= 5 * seed_len {
        vec![
            0,
            qlen / 4,
            qlen / 2 - seed_len / 2,
            (3 * qlen) / 4 - seed_len,
            qlen - seed_len,
        ]
    } else {
        vec![0, (qlen - seed_len) / 2, qlen - seed_len]
    };

    let max_mm = (qlen / 20).max(5); // allow up to 5% SNPs or 5 bp
    let mut best: Option<(usize, usize, Strand)> = None; // (mismatches, offset, strand)

    for (strand, query) in [(Strand::Plus, query_plus), (Strand::Minus, query_minus)] {
        for &s_off in &seed_offsets {
            let seed = &query[s_off..s_off + seed_len];
            if seed.contains(&b'N') {
                continue;
            }
            for hit in memmem::find_iter(haystack, seed) {
                if hit < s_off {
                    continue;
                }
                let start = hit - s_off;
                if start + qlen > haystack.len() {
                    continue;
                }
                let window = &haystack[start..start + qlen];
                let mut mm = 0usize;
                for i in 0..qlen {
                    if window[i] != query[i] {
                        mm += 1;
                        if mm > max_mm {
                            break;
                        }
                    }
                }
                if mm <= max_mm {
                    let replace = match best {
                        None => true,
                        Some((prev_mm, _, _)) => mm < prev_mm,
                    };
                    if replace {
                        best = Some((mm, start, strand));
                        if mm == 0 {
                            return Some((start, strand));
                        }
                    }
                }
            }
        }
    }

    best.map(|(_, offset, strand)| (offset, strand))
}

/// Parse a `chr:start-end` region string (1-based inclusive).
pub fn parse_region_spec(spec: &str) -> Result<(String, usize, usize)> {
    let (chrom, coords) = spec
        .rsplit_once(':')
        .ok_or_else(|| anyhow!("Invalid region '{spec}': expected format chr:start-end"))?;
    let (s_str, e_str) = coords
        .split_once('-')
        .ok_or_else(|| anyhow!("Invalid region '{spec}': expected format chr:start-end"))?;
    let start: usize = s_str
        .replace(',', "")
        .parse()
        .with_context(|| format!("Invalid start coordinate in region '{spec}'"))?;
    let end: usize = e_str
        .replace(',', "")
        .parse()
        .with_context(|| format!("Invalid end coordinate in region '{spec}'"))?;
    if start == 0 || end < start {
        bail!("Invalid region '{spec}': start must be >= 1 and <= end");
    }
    Ok((chrom.to_string(), start, end))
}

/// Try to extract an optional `chr:start-end` hint from a FASTA header line.
pub fn extract_header_region_hint(header: &str) -> Option<(String, usize, usize)> {
    for token in header.split_whitespace() {
        let cleaned = token.trim_end_matches("(-)").trim_end_matches("(+)");
        if let Ok(reg) = parse_region_spec(cleaned) {
            return Some(reg);
        }
    }
    None
}

/// Read all records `(name, normalized_seq)` from a FASTA file.
pub fn read_fasta_records(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let file =
        File::open(path).with_context(|| format!("Failed to open FASTA {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                records.push((name, std::mem::take(&mut current_seq)));
            }
            let rec_name = header.trim().to_string();
            current_name = Some(rec_name);
        } else {
            for &b in trimmed.as_bytes() {
                if !b.is_ascii_whitespace() {
                    current_seq.push(match b.to_ascii_uppercase() {
                        b'A' => b'A',
                        b'C' => b'C',
                        b'G' => b'G',
                        b'T' => b'T',
                        _ => b'N',
                    });
                }
            }
        }
    }
    if let Some(name) = current_name {
        records.push((name, current_seq));
    }
    if records.is_empty() {
        bail!("No FASTA sequences found in {}", path.display());
    }
    Ok(records)
}

/// Parsed BED interval `(chrom, start_1based, end_1based, optional_name)`.
pub type BedInterval = (String, usize, usize, Option<String>);

/// Parse a 3+ column BED file (0-based half-open `[start, end)` -> 1-based inclusive `[start+1, end]`).
pub fn read_bed_regions(path: &Path) -> Result<Vec<BedInterval>> {
    let file =
        File::open(path).with_context(|| format!("Failed to open BED file {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    for (line_no, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        let cols = if cols.len() < 3 {
            trimmed.split_whitespace().collect::<Vec<_>>()
        } else {
            cols
        };
        if cols.len() < 3 {
            bail!(
                "Invalid BED line {} in {}: expected >= 3 columns",
                line_no + 1,
                path.display()
            );
        }
        let chrom = cols[0].to_string();
        let start0: usize = cols[1]
            .parse()
            .with_context(|| format!("Invalid BED start at line {}", line_no + 1))?;
        let end: usize = cols[2]
            .parse()
            .with_context(|| format!("Invalid BED end at line {}", line_no + 1))?;
        if end <= start0 {
            bail!(
                "Invalid BED interval {}:{}-{} at line {}",
                chrom,
                start0,
                end,
                line_no + 1
            );
        }
        let label = cols.get(3).map(|s| (*s).to_string());
        regions.push((chrom, start0 + 1, end, label));
    }
    if regions.is_empty() {
        bail!("No target intervals found in BED file {}", path.display());
    }
    Ok(regions)
}

/// Resolve all target loci from CLI inputs (`--target-fasta`, `--region`, `--bed`, or built-in FLT3).
pub fn resolve_target_loci(
    reference: &ReferenceGenome,
    genome_build: GenomeBuild,
    target_fasta: Option<&Path>,
    regions: &[String],
    bed: Option<&Path>,
) -> Result<Vec<TargetLocus>> {
    let mut loci = Vec::new();

    // Case 1: User-provided target FASTA (e.g. FLT3_dna_e14e15.fa)
    if let Some(tf_path) = target_fasta {
        let records = read_fasta_records(tf_path)?;
        let parsed_hints: Vec<(String, usize, usize)> = regions
            .iter()
            .map(|r| parse_region_spec(r))
            .collect::<Result<_>>()?;

        for (idx, (header, seq)) in records.into_iter().enumerate() {
            let short_name = header
                .split_whitespace()
                .next()
                .unwrap_or("target")
                .to_string();
            let hint = parsed_hints
                .get(idx)
                .cloned()
                .or_else(|| extract_header_region_hint(&header));
            let locus = reference.locate_target_sequence(&short_name, &seq, hint.as_ref())?;
            loci.push(locus);
        }
        return Ok(loci);
    }

    // Case 2: Explicit --region arguments
    for reg_str in regions {
        let (chrom, start, end) = parse_region_spec(reg_str)?;
        let info = reference
            .resolve_contig(&chrom)
            .ok_or_else(|| anyhow!("Contig '{chrom}' not found in reference FASTA"))?;
        let seq = reference.fetch(&info.name, start, end)?;
        loci.push(TargetLocus {
            name: format!("{}:{}-{}", info.name, start, end),
            chrom: info.name.clone(),
            start,
            end,
            source_strand: Strand::Plus,
            seq,
        });
    }

    // Case 3: Explicit --bed file
    if let Some(bed_path) = bed {
        for (chrom, start, end, label) in read_bed_regions(bed_path)? {
            let info = reference
                .resolve_contig(&chrom)
                .ok_or_else(|| anyhow!("Contig '{chrom}' not found in reference FASTA"))?;
            let seq = reference.fetch(&info.name, start, end)?;
            let name = label.unwrap_or_else(|| format!("{}:{}-{}", info.name, start, end));
            loci.push(TargetLocus {
                name,
                chrom: info.name.clone(),
                start,
                end,
                source_strand: Strand::Plus,
                seq,
            });
        }
    }

    // Case 4: Default built-in classic FLT3 locus
    if loci.is_empty() {
        loci.push(reference.resolve_builtin_flt3(genome_build)?);
    }

    Ok(loci)
}

fn parse_fai(path: &Path) -> Result<(Vec<ContigInfo>, HashMap<String, usize>)> {
    let file =
        File::open(path).with_context(|| format!("Failed to open FAI {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut contigs = Vec::new();
    let mut map = HashMap::new();
    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            bail!("Malformed FAI line {} in {}", idx + 1, path.display());
        }
        let name = cols[0].to_string();
        let length: usize = cols[1].parse()?;
        let offset: u64 = cols[2].parse()?;
        let line_bases: usize = cols[3].parse()?;
        let line_width: usize = cols[4].parse()?;
        map.insert(name.clone(), contigs.len());
        contigs.push(ContigInfo {
            name,
            length,
            offset,
            line_bases,
            line_width,
        });
    }
    Ok((contigs, map))
}

fn read_up_to<R: Read>(reader: &mut R, buf: &mut [u8]) -> Result<usize> {
    let mut total = 0usize;
    while total < buf.len() {
        let n = reader.read(&mut buf[total..])?;
        if n == 0 {
            break;
        }
        total += n;
    }
    Ok(total)
}
