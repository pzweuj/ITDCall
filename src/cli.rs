use std::path::PathBuf;

use clap::{Parser, ValueEnum};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AssayMethod {
    /// Hybrid capture / probe-based enrichment (random fragment ends; boundary-averaged unbiased VAF)
    #[value(alias = "capture")]
    Hc,
    /// Amplicon sequencing (fixed primer ends; fragment-level Mut / (WT + Mut) VAF)
    #[value(alias = "amp")]
    Amplicon,
}

impl AssayMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            AssayMethod::Hc => "hc",
            AssayMethod::Amplicon => "amplicon",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum GenomeBuild {
    /// Auto-detect hg19 vs hg38 from reference genome / BAM header
    Auto,
    /// GRCh37 / hg19 built-in FLT3 locus (chr13:28607438-28608937)
    #[value(alias = "grch37")]
    Hg19,
    /// GRCh38 / hg38 built-in FLT3 locus (chr13:28033301-28034800)
    #[value(alias = "grch38")]
    Hg38,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ReadMode {
    /// Extract only reads aligned to the target locus (fastest; uses BAI index if present)
    Targeted,
    /// Extract target-locus reads plus unmapped reads (recovers heavily diverged/long ITD reads)
    Loose,
    /// Scan all reads in the BAM file
    All,
}

#[derive(Debug, Parser)]
#[command(
    name = "itdcall",
    version,
    about = "Pure-Rust NGS Internal Tandem Duplication (ITD) caller for Hybrid Capture and Amplicon BAMs"
)]
pub struct Cli {
    /// Input BAM file (.bai index recommended for targeted mode)
    #[arg(short = 'b', long = "bam")]
    pub bam: PathBuf,

    /// Reference genome FASTA file (.fai index recommended; auto-indexed if missing)
    #[arg(short = 'r', long = "reference")]
    pub reference: PathBuf,

    /// Output VCF file path (use '-' for stdout)
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Sequencing assay methodology for ITD quantification
    #[arg(short = 'm', long = "method", value_enum, default_value_t = AssayMethod::Hc)]
    pub method: AssayMethod,

    /// Reference genome build for built-in classic FLT3 locus
    #[arg(short = 'g', long = "genome", value_enum, default_value_t = GenomeBuild::Auto)]
    pub genome: GenomeBuild,

    /// User-provided target region FASTA file (e.g. FLT3_dna_e14e15.fa) to override/calibrate target locus
    #[arg(long = "target-fasta")]
    pub target_fasta: Option<PathBuf>,

    /// Custom target genomic region(s) in chr:start-end (1-based inclusive), repeatable
    #[arg(long = "region")]
    pub regions: Vec<String>,

    /// Custom BED file specifying target region(s)
    #[arg(long = "bed")]
    pub bed: Option<PathBuf>,

    /// BAM read extraction mode
    #[arg(long = "read-mode", value_enum, default_value_t = ReadMode::Targeted)]
    pub read_mode: ReadMode,

    /// Sample name in VCF output (defaults to BAM @RG SM tag or BAM filename stem)
    #[arg(long = "sample")]
    pub sample: Option<String>,

    /// Disable built-in 3' Illumina adapter trimming
    #[arg(long = "no-trim-adapters", default_value_t = false)]
    pub no_trim_adapters: bool,

    /// Exclude BAM reads flagged as PCR/optical duplicates (0x400); keep disabled for amplicon data
    #[arg(long = "exclude-duplicates", default_value_t = false)]
    pub exclude_duplicates: bool,

    /// Minimum effective supporting mutant read/fragment depth for PASS
    #[arg(long = "min-reads", default_value_t = 2.0)]
    pub min_reads: f64,

    /// Minimum variant allele frequency (AF) for PASS
    #[arg(long = "min-vaf", default_value_t = 0.0)]
    pub min_vaf: f64,

    /// Minimum net ITD insertion length (bp)
    #[arg(long = "min-itd-len", default_value_t = 3)]
    pub min_itd_len: usize,

    /// Maximum net ITD insertion length (bp)
    #[arg(long = "max-itd-len", default_value_t = 500)]
    pub max_itd_len: usize,

    /// Minimum soft-clip length (bp) to attempt split-read clip anchoring
    #[arg(long = "min-clip-len", default_value_t = 9)]
    pub min_clip_len: usize,

    /// Minimum bp extending past a breakpoint required to count WT/Mut support
    #[arg(long = "buffer", default_value_t = 10)]
    pub buffer: usize,

    /// Maximum edit distance (NM) allowed for stringent WT or Mut read assignment
    #[arg(long = "max-edit-dist", default_value_t = 5)]
    pub max_edit_dist: usize,

    /// Maximum soft-clip length allowed for stringent WT or Mut read assignment
    #[arg(long = "max-clip", default_value_t = 3)]
    pub max_clip: usize,

    /// Maximum total insertion bases allowed for stringent WT or Mut read assignment
    #[arg(long = "max-ins", default_value_t = 2)]
    pub max_ins: usize,

    /// Maximum total deletion bases allowed for stringent WT or Mut read assignment
    #[arg(long = "max-del", default_value_t = 2)]
    pub max_del: usize,

    /// Maximum edit distance for greedy candidate allele clustering
    #[arg(long = "cluster-edit-dist", default_value_t = 5)]
    pub cluster_edit_dist: usize,

    /// Only report in-frame ITDs (SVLEN % 3 == 0)
    #[arg(long = "in-frame-only", default_value_t = false)]
    pub in_frame_only: bool,

    /// Include non-PASS filtered candidates in the output VCF
    #[arg(long = "include-filtered", default_value_t = false)]
    pub include_filtered: bool,

    /// Number of worker threads (0 = auto)
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    pub threads: usize,
}
