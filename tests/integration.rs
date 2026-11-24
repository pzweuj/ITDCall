use std::fs::{self, File};
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;

use noodles::bam;
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record_buf::{Cigar as CigarBuf, QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::map::{self, header::tag::SORT_ORDER};
use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
use tempfile::TempDir;

use itdcall::cli::{AssayMethod, Cli, GenomeBuild, ReadMode};
use itdcall::pipeline::run_cli;
use itdcall::seq::revcomp;

/// Exact 1500 bp minus-strand sequence from FLT3_ITD_ext (`FLT3_dna_e14e15.fa`).
const FLT3_E14E15_MINUS: &[u8] = "TGAATTTTTTTAGAGACGGGATTTCACTGTGTTGCCTGGGCTGGTGTCTAATTCTTGGACTCAAGTGATCCTCCCATCTTGGCCTCCCAAAGTGCTGGGACTACAGGCGTGAGCCATCGTGCCTGGACCAGATGTAATTCTTGTAAATTTACCCTAGTAAATGGTTTTAAAGAATGGATTGAGAATAATCCGCAATTTTCTAGGGAGGAAATGTATCCTTTACAGGGACATTGCCTGATTGTCTGTGGGGAGTAGAGTATATGTAGAGTGGTTGTTAGGACTGAAAATGATTATTACTGAAACAGGATGTGAGAGATTATAATGAGTTGTCCACTATTTATAATGTCACACAGGAATTCTGTTTCATCGCTGAGTGACACTCTTTTGTTGCAGGCCCCTTCCCTTTCATCCAAGACAACATCTCATTCTATGCAACAATTGGTGTTTGTCTCCTCTTCATTGTCGTTTTAACCCTGCTAATTTGTCACAAGTACAAAAAGGTAAAAGCAAAGGTAAAAATTCATTATTCTTTCCTCTATCTGCAGAACTGCCTATTCCTAACTGACTCATCATTTCATCTCTGAAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAGGTACAGTATAGTGGAAGGACAGCAACAAAGATGCACAAAAATGGGAGGCACAGTTTCCCACCCATGCCTTCTTCTCTTTTCCATCCTTTTAATGGTTACTGTTTGCCATGTTTCAAGGCTAAAAATGGAGTGGATTGGGGTGTCAACCCAGTCATGAATACAAAATTCAAGTCATAAATACAAAACTCCACCTCTTCTTCCACTTCTCTCTTTGTTTATTTTTATTTTCTTTTATTTTTTTGAGACAGGGTCTCACTCTGTCACCTAGGCTGGAGTACAGTAGTGTGATCATGGCTCACTGCAGCCTCCACCCCCTGGGCTCAGGCTCCTCAGCAGCTGGGGCTACAGGCGTGCACCACCACGCTCAGCTAAGTTTTTTATTTTTGGTAGAGACGGGGTTTCACCATGTTGCCTAGGCTGATCTGGGACTCCTGAGCTCAAACTATCTGCCCGCCTCGGCCTCCCAAAGTGTTGGGATTACAGGCGTGAGCCACTGTGCCAGGCCTATTTAGTTTTTATAGGCTTTGATATTTGGGGATACAGTTGCCACAGGCAGAAGTTCTTTCTTACATTTGTGTTGTTTTAT".as_bytes();

/// Write a single-contig reference FASTA and matching .fai index.
fn write_ref_fasta(dir: &Path, contig: &str, seq: &[u8]) -> std::path::PathBuf {
    let fa_path = dir.join("ref.fa");
    let mut f = File::create(&fa_path).unwrap();
    writeln!(f, ">{contig}").unwrap();
    let line_bases = 80usize;
    let header_len = contig.len() + 2; // ">contig\n"
    for chunk in seq.chunks(line_bases) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    drop(f);

    let fai_path = dir.join("ref.fa.fai");
    let mut fai = File::create(&fai_path).unwrap();
    writeln!(
        fai,
        "{contig}\t{}\t{header_len}\t{line_bases}\t{}",
        seq.len(),
        line_bases + 1
    )
    .unwrap();
    fa_path
}

struct SimRead {
    name: String,
    pos_1based: usize, // 1-based alignment start on contig (approximate, realigned by itdcall)
    seq_fwd_on_ref: Vec<u8>,
    is_paired: bool,
    is_r1: bool,
    is_rev: bool,
}

/// Write a synthetic BAM file containing simulated reads.
fn write_synthetic_bam(
    dir: &Path,
    filename: &str,
    contig: &str,
    contig_len: usize,
    reads: &[SimRead],
) -> std::path::PathBuf {
    let bam_path = dir.join(filename);
    let mut header_builder = Header::builder().set_header(
        Map::<map::Header>::builder()
            .insert(SORT_ORDER, "coordinate")
            .build()
            .unwrap(),
    );
    header_builder = header_builder.add_reference_sequence(
        contig,
        Map::<ReferenceSequence>::new(NonZeroUsize::new(contig_len).unwrap()),
    );
    let header = header_builder.build();

    let file = File::create(&bam_path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    let mut sorted_reads: Vec<_> = reads.iter().collect();
    sorted_reads.sort_by_key(|r| r.pos_1based);
    for r in sorted_reads {
        let mut flags = Flags::empty();
        if r.is_paired {
            flags |= Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
            if r.is_r1 {
                flags |= Flags::FIRST_SEGMENT;
            } else {
                flags |= Flags::LAST_SEGMENT;
            }
        }
        if r.is_rev {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        let seq_len = r.seq_fwd_on_ref.len();
        let rec = RecordBuf::builder()
            .set_name(r.name.as_bytes())
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(r.pos_1based.max(1)).unwrap())
            .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(CigarBuf::from(vec![Op::new(Kind::Match, seq_len)]))
            .set_sequence(Sequence::from(r.seq_fwd_on_ref.clone()))
            .set_quality_scores(QualityScores::from(vec![35u8; seq_len]))
            .build();
        writer.write_alignment_record(&header, &rec).unwrap();
    }
    writer.try_finish().unwrap();
    bam_path
}

fn index_bam(path: &Path) {
    let bam::Index::Bai(index) = bam::fs::index(path).unwrap() else {
        panic!("expected a BAI index");
    };
    bam::bai::fs::write(path.with_extension("bam.bai"), &index).unwrap();
}

fn single_end_insertion_reads(wt: &[u8], insertion: &[u8]) -> Vec<SimRead> {
    let mut seq = wt[700..850].to_vec();
    seq.splice(59..59, insertion.iter().copied());
    (0..6)
        .map(|i| SimRead {
            name: format!("insertion_{i}"),
            pos_1based: 701,
            seq_fwd_on_ref: seq.clone(),
            is_paired: false,
            is_r1: false,
            is_rev: false,
        })
        .collect()
}

/// Generate paired-end reads tiled uniformly across a template (hybrid capture model).
fn simulate_hc_pairs(
    prefix: &str,
    template: &[u8],
    contig_offset: usize,
    read_len: usize,
    insert_size: usize,
    step: usize,
    copies: usize,
) -> Vec<SimRead> {
    let mut reads = Vec::new();
    let mut idx = 0usize;
    let max_start = template.len().saturating_sub(insert_size);
    for _ in 0..copies {
        let mut start = 10usize;
        while start <= max_start.saturating_sub(10) {
            let r1_seq = template[start..start + read_len].to_vec();
            let r2_start = start + insert_size - read_len;
            let r2_seq = template[r2_start..r2_start + read_len].to_vec();
            let qname = format!("{prefix}_{idx}");
            reads.push(SimRead {
                name: qname.clone(),
                pos_1based: contig_offset + start + 1,
                seq_fwd_on_ref: r1_seq,
                is_paired: true,
                is_r1: true,
                is_rev: false,
            });
            reads.push(SimRead {
                name: qname,
                pos_1based: contig_offset + r2_start + 1,
                seq_fwd_on_ref: r2_seq,
                is_paired: true,
                is_r1: false,
                is_rev: true,
            });
            idx += 1;
            start += step;
        }
    }
    reads
}

/// Generate fixed-endpoint paired-end amplicon reads.
fn simulate_amplicon_pairs(
    prefix: &str,
    template: &[u8],
    contig_offset: usize,
    amp_start: usize,
    amp_end: usize,
    read_len: usize,
    count: usize,
) -> Vec<SimRead> {
    let mut reads = Vec::new();
    let r1_seq = template[amp_start..amp_start + read_len].to_vec();
    let r2_start = amp_end - read_len;
    let r2_seq = template[r2_start..amp_end].to_vec();
    for idx in 0..count {
        let qname = format!("{prefix}_{idx}");
        reads.push(SimRead {
            name: qname.clone(),
            pos_1based: contig_offset + amp_start + 1,
            seq_fwd_on_ref: r1_seq.clone(),
            is_paired: true,
            is_r1: true,
            is_rev: false,
        });
        reads.push(SimRead {
            name: qname,
            pos_1based: contig_offset + r2_start + 1,
            seq_fwd_on_ref: r2_seq.clone(),
            is_paired: true,
            is_r1: false,
            is_rev: true,
        });
    }
    reads
}

#[derive(Debug)]
struct ParsedVcfRow {
    chrom: String,
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    filter: String,
    info: String,
}

fn parse_vcf_rows(vcf_path: &Path) -> Vec<ParsedVcfRow> {
    let text = fs::read_to_string(vcf_path).unwrap();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 10, "VCF line must have 10 columns: {line}");
        rows.push(ParsedVcfRow {
            chrom: cols[0].to_string(),
            pos: cols[1].parse().unwrap(),
            ref_allele: cols[3].to_string(),
            alt_allele: cols[4].to_string(),
            filter: cols[6].to_string(),
            info: cols[7].to_string(),
        });
    }
    rows
}

fn get_info_val<'a>(info: &'a str, key: &str) -> &'a str {
    for kv in info.split(';') {
        if let Some((k, v)) = kv.split_once('=')
            && k == key
        {
            return v;
        }
    }
    panic!("Missing INFO key {key} in {info}");
}

fn default_cli(
    bam: std::path::PathBuf,
    reference: std::path::PathBuf,
    output: std::path::PathBuf,
) -> Cli {
    Cli {
        bam,
        reference,
        output,
        method: AssayMethod::Hc,
        genome: GenomeBuild::Auto,
        target_fasta: None,
        regions: Vec::new(),
        bed: None,
        read_mode: ReadMode::Targeted,
        sample: Some("TEST_SAMPLE".to_string()),
        no_trim_adapters: false,
        exclude_duplicates: false,
        min_reads: 2.0,
        min_vaf: 0.0,
        min_itd_len: 3,
        max_itd_len: 500,
        min_clip_len: 9,
        buffer: 10,
        max_edit_dist: 5,
        max_clip: 3,
        max_ins: 2,
        max_del: 2,
        cluster_edit_dist: 5,
        in_frame_only: false,
        include_filtered: false,
        threads: 2,
    }
}

#[test]
fn test_hc_short_and_long_itd_and_user_minus_strand_target_fasta() {
    let tmp = TempDir::new().unwrap();
    // Plus-strand reference contig containing 500 bp padding + FLT3 (+) + 500 bp padding
    let flt3_plus = revcomp(FLT3_E14E15_MINUS);
    let pad_left = vec![b'A'; 500];
    let pad_right = vec![b'C'; 500];
    let mut contig_seq = Vec::new();
    contig_seq.extend_from_slice(&pad_left);
    contig_seq.extend_from_slice(&flt3_plus);
    contig_seq.extend_from_slice(&pad_right);

    let ref_fa = write_ref_fasta(tmp.path(), "chr13", &contig_seq);

    // User-provided calibration target FASTA on MINUS strand (just like FLT3_dna_e14e15.fa!)
    let target_fa = tmp.path().join("FLT3_dna_e14e15.fa");
    {
        let mut f = File::create(&target_fa).unwrap();
        writeln!(f, ">FLT3_dna_e14e15").unwrap();
        f.write_all(FLT3_E14E15_MINUS).unwrap();
        writeln!(f).unwrap();
    }

    // Create a LONG 180 bp ITD (longer than 120 bp read length!): duplicate flt3_plus[600..780]
    let mut mut_long = Vec::new();
    mut_long.extend_from_slice(&flt3_plus[..780]);
    mut_long.extend_from_slice(&flt3_plus[600..780]);
    mut_long.extend_from_slice(&flt3_plus[780..]);

    // Simulate 50% VAF (equal copies of WT and Mut) under Hybrid Capture
    let mut reads = simulate_hc_pairs("wt", &flt3_plus, 500, 120, 220, 15, 2);
    reads.extend(simulate_hc_pairs("mut", &mut_long, 500, 120, 220, 15, 2));

    let bam = write_synthetic_bam(tmp.path(), "hc_long.bam", "chr13", contig_seq.len(), &reads);
    let vcf = tmp.path().join("hc_long.vcf");

    let mut cli = default_cli(bam, ref_fa, vcf.clone());
    cli.method = AssayMethod::Hc;
    cli.target_fasta = Some(target_fa);
    run_cli(&cli).unwrap();

    let rows = parse_vcf_rows(&vcf);
    assert_eq!(rows.len(), 1, "Expected exactly 1 ITD call, got {rows:?}");
    let r = &rows[0];
    assert_eq!(r.chrom, "chr13");
    assert_eq!(r.filter, "PASS");
    assert_eq!(get_info_val(&r.info, "SVLEN"), "180");
    assert_eq!(get_info_val(&r.info, "METHOD"), "hc");
    assert_eq!(r.alt_allele.len() - r.ref_allele.len(), 180);

    // Verify genomic allele reconstruction equals mut_long byte-for-byte!
    let local_pos = r.pos - 501; // 0-based inside flt3_plus
    let mut reconstructed = Vec::new();
    reconstructed.extend_from_slice(&flt3_plus[..local_pos]);
    reconstructed.extend_from_slice(r.alt_allele.as_bytes());
    reconstructed.extend_from_slice(&flt3_plus[local_pos + r.ref_allele.len()..]);
    assert_eq!(reconstructed, mut_long);

    // Verify HC unbiased VAF is near 0.50 despite 180 bp ITD > 120 bp read length!
    let af: f64 = get_info_val(&r.info, "AF").parse().unwrap();
    assert!(
        (af - 0.50).abs() < 0.12,
        "Expected unbiased HC AF ~0.50 for 180 bp ITD, got {af}"
    );
    // Verify NO annotation fields (CDS, AA, GENE) are emitted
    assert!(!r.info.contains("CDS="));
    assert!(!r.info.contains("AA="));
    assert!(!r.info.contains("GENE="));
}

#[test]
fn test_amplicon_complex_itd_exact_vaf_and_builtin_auto_flt3() {
    let tmp = TempDir::new().unwrap();
    let flt3_plus = revcomp(FLT3_E14E15_MINUS);
    // Built-in auto-detection should find FLT3 even on a cropped contig without --target-fasta!
    let ref_fa = write_ref_fasta(tmp.path(), "chr13", &flt3_plus);

    // Complex ITD: 36 bp dup of flt3_plus[780..816] + 3 bp non-templated insertion "CCG" -> SVLEN = 39
    // Internal junction sits at 816..819.
    let mut mut_complex = Vec::new();
    mut_complex.extend_from_slice(&flt3_plus[..816]);
    mut_complex.extend_from_slice(b"CCG");
    mut_complex.extend_from_slice(&flt3_plus[780..]);

    // Realistic clinical amplicon spanning [700..940] on WT (240 bp) -> [700..979] on Mut (279 bp),
    // sequenced with 2x180 bp PE reads so both R1 and R2 span across the 816..819 junction.
    // 60 WT pairs + 40 Mut pairs => true AF = 0.40, true AR = 40/60 = 0.6667
    let mut reads = simulate_amplicon_pairs("wt", &flt3_plus, 0, 700, 940, 180, 60);
    reads.extend(simulate_amplicon_pairs(
        "mut",
        &mut_complex,
        0,
        700,
        979,
        180,
        40,
    ));

    let bam = write_synthetic_bam(tmp.path(), "amp.bam", "chr13", flt3_plus.len(), &reads);
    let vcf = tmp.path().join("amp.vcf");

    let mut cli = default_cli(bam, ref_fa, vcf.clone());
    cli.method = AssayMethod::Amplicon;
    cli.genome = GenomeBuild::Auto;
    run_cli(&cli).unwrap();

    let rows = parse_vcf_rows(&vcf);
    assert_eq!(rows.len(), 1, "Expected 1 complex ITD call, got {rows:?}");
    let r = &rows[0];
    assert_eq!(get_info_val(&r.info, "SVLEN"), "39");
    assert_eq!(get_info_val(&r.info, "DUPLEN"), "36");
    assert_eq!(get_info_val(&r.info, "INSLEN"), "3");
    assert_eq!(get_info_val(&r.info, "METHOD"), "amplicon");

    let af: f64 = get_info_val(&r.info, "AF").parse().unwrap();
    let ar: f64 = get_info_val(&r.info, "AR").parse().unwrap();
    assert!(
        (af - 0.40).abs() < 1e-3,
        "Expected exact amplicon AF=0.40, got {af}"
    );
    assert!(
        (ar - 0.6667).abs() < 1e-2,
        "Expected exact amplicon AR=0.6667, got {ar}"
    );
}

#[test]
fn test_multiclone_and_single_end_and_wt_negative_control() {
    let tmp = TempDir::new().unwrap();
    let flt3_plus = revcomp(FLT3_E14E15_MINUS);
    let ref_fa = write_ref_fasta(tmp.path(), "13", &flt3_plus);

    // Clone A: 24 bp dup [760..784]; Clone B: 60 bp dup [790..850]
    let mut mut_a = Vec::new();
    mut_a.extend_from_slice(&flt3_plus[..784]);
    mut_a.extend_from_slice(&flt3_plus[760..784]);
    mut_a.extend_from_slice(&flt3_plus[784..]);

    let mut mut_b = Vec::new();
    mut_b.extend_from_slice(&flt3_plus[..850]);
    mut_b.extend_from_slice(&flt3_plus[790..850]);
    mut_b.extend_from_slice(&flt3_plus[850..]);

    // 50 WT + 25 Clone A + 25 Clone B (amplicon [680..940], 2x220 bp MiSeq PE reads)
    let mut reads = simulate_amplicon_pairs("wt", &flt3_plus, 0, 680, 940, 220, 50);
    reads.extend(simulate_amplicon_pairs(
        "cloneA", &mut_a, 0, 680, 964, 220, 25,
    ));
    reads.extend(simulate_amplicon_pairs(
        "cloneB", &mut_b, 0, 680, 1000, 220, 25,
    ));

    let bam = write_synthetic_bam(tmp.path(), "multi.bam", "13", flt3_plus.len(), &reads);
    let vcf = tmp.path().join("multi.vcf");

    let mut cli = default_cli(bam, ref_fa.clone(), vcf.clone());
    cli.method = AssayMethod::Amplicon;
    run_cli(&cli).unwrap();

    let rows = parse_vcf_rows(&vcf);
    assert_eq!(
        rows.len(),
        2,
        "Expected 2 multi-clone ITD calls, got {rows:?}"
    );
    let svlens: Vec<&str> = rows
        .iter()
        .map(|r| get_info_val(&r.info, "SVLEN"))
        .collect();
    assert!(svlens.contains(&"24"));
    assert!(svlens.contains(&"60"));
    for r in &rows {
        let af: f64 = get_info_val(&r.info, "AF").parse().unwrap();
        let ar: f64 = get_info_val(&r.info, "AR").parse().unwrap();
        assert!((af - 0.25).abs() < 1e-2, "Expected AF=0.25, got {af}");
        assert!(
            (ar - 0.50).abs() < 1e-2,
            "Expected multi-clone AR=0.50, got {ar}"
        );
    }

    reads.retain(|read| !read.name.starts_with("wt_"));
    cli.bam = write_synthetic_bam(tmp.path(), "multi_no_wt.bam", "13", flt3_plus.len(), &reads);
    run_cli(&cli).unwrap();
    let rows = parse_vcf_rows(&vcf);
    assert_eq!(rows.len(), 2);
    for row in &rows {
        assert_eq!(get_info_val(&row.info, "AF").parse::<f64>().unwrap(), 0.5);
        assert_eq!(get_info_val(&row.info, "AR"), "Inf");
    }

    // Single-end reads carrying a 9 bp short ITD (tests inline insertion + single-end path)
    let mut mut_short = Vec::new();
    mut_short.extend_from_slice(&flt3_plus[..759]);
    mut_short.extend_from_slice(&flt3_plus[750..759]);
    mut_short.extend_from_slice(&flt3_plus[759..]);
    let mut se_reads = Vec::new();
    for i in 0..30 {
        se_reads.push(SimRead {
            name: format!("se_wt_{i}"),
            pos_1based: 701,
            seq_fwd_on_ref: flt3_plus[700..850].to_vec(),
            is_paired: false,
            is_r1: false,
            is_rev: false,
        });
    }
    for i in 0..10 {
        se_reads.push(SimRead {
            name: format!("se_mut_{i}"),
            pos_1based: 701,
            seq_fwd_on_ref: mut_short[700..859].to_vec(),
            is_paired: false,
            is_r1: false,
            is_rev: false,
        });
    }
    let bam_se = write_synthetic_bam(tmp.path(), "se.bam", "13", flt3_plus.len(), &se_reads);
    let vcf_se = tmp.path().join("se.vcf");
    let mut cli_se = default_cli(bam_se, ref_fa.clone(), vcf_se.clone());
    cli_se.method = AssayMethod::Amplicon;
    cli_se.regions = vec!["13:1-1500".to_string()];
    run_cli(&cli_se).unwrap();
    let se_rows = parse_vcf_rows(&vcf_se);
    assert_eq!(
        se_rows.len(),
        1,
        "Expected 1 short single-end ITD call, got {se_rows:?}"
    );
    assert_eq!(get_info_val(&se_rows[0].info, "SVLEN"), "9");
    let se_af: f64 = get_info_val(&se_rows[0].info, "AF").parse().unwrap();
    assert!(
        (se_af - 0.25).abs() < 1e-3,
        "Expected single-end AF=0.25, got {se_af}"
    );

    // Wild-type negative control
    let wt_only = simulate_hc_pairs("wt", &flt3_plus, 0, 120, 220, 20, 1);
    let bam_wt = write_synthetic_bam(tmp.path(), "wt.bam", "13", flt3_plus.len(), &wt_only);
    let vcf_wt = tmp.path().join("wt.vcf");
    let cli_wt = default_cli(bam_wt, ref_fa, vcf_wt.clone());
    run_cli(&cli_wt).unwrap();
    let wt_rows = parse_vcf_rows(&vcf_wt);
    assert!(
        wt_rows.is_empty(),
        "Expected 0 calls on WT negative control, got {wt_rows:?}"
    );
}

#[test]
fn test_rejects_non_templated_insertion() {
    let tmp = TempDir::new().unwrap();
    let wt = revcomp(FLT3_E14E15_MINUS);
    let reference = write_ref_fasta(tmp.path(), "13", &wt);
    let vcf = tmp.path().join("out.vcf");
    let mut cli = default_cli(tmp.path().join("novel.bam"), reference, vcf.clone());
    cli.method = AssayMethod::Amplicon;
    cli.include_filtered = true;
    // The second insertion shares a single base with the adjacent reference.
    for insertion in [b"CTGCGGGCGGGA", b"TTGCGCGCTTGC"] {
        let reads = single_end_insertion_reads(&wt, insertion);
        cli.bam = write_synthetic_bam(tmp.path(), "novel.bam", "13", wt.len(), &reads);
        run_cli(&cli).unwrap();
        assert!(parse_vcf_rows(&vcf).is_empty());
    }
}

#[test]
fn test_variant_ref_uses_genomic_fasta_with_custom_target_snp() {
    let tmp = TempDir::new().unwrap();
    let wt = revcomp(FLT3_E14E15_MINUS);
    let reference = write_ref_fasta(tmp.path(), "13", &wt);
    let mut target = wt.clone();
    target[749] = b'A';
    let target_path = tmp.path().join("target.fa");
    fs::write(
        &target_path,
        format!(
            ">target 13:1-1500\n{}\n",
            std::str::from_utf8(&target).unwrap()
        ),
    )
    .unwrap();
    let reads = single_end_insertion_reads(&target, &target[750..759]);
    let bam = write_synthetic_bam(tmp.path(), "target.bam", "13", wt.len(), &reads);
    let vcf = tmp.path().join("out.vcf");
    let mut cli = default_cli(bam, reference, vcf.clone());
    cli.method = AssayMethod::Amplicon;
    cli.target_fasta = Some(target_path);
    run_cli(&cli).unwrap();
    let rows = parse_vcf_rows(&vcf);
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(get_info_val(&row.info, "SVLEN"), "9");
    assert_eq!(
        row.ref_allele.as_bytes(),
        &wt[row.pos - 1..][..row.ref_allele.len()]
    );
}

#[test]
fn test_targeted_bam_index_preserves_calls_and_missing_contig_error() {
    let tmp = TempDir::new().unwrap();
    let wt = revcomp(FLT3_E14E15_MINUS);
    let reference = write_ref_fasta(tmp.path(), "13", &wt);
    let mut reads = single_end_insertion_reads(&wt, &wt[750..759]);
    let mates: Vec<_> = reads
        .iter_mut()
        .map(|read| {
            read.is_paired = true;
            read.is_r1 = true;
            SimRead {
                name: read.name.clone(),
                pos_1based: 951,
                seq_fwd_on_ref: wt[950..1030].to_vec(),
                is_paired: true,
                is_r1: false,
                is_rev: true,
            }
        })
        .collect();
    reads.extend(mates);
    let bam = write_synthetic_bam(tmp.path(), "paired.bam", "13", wt.len(), &reads);
    let vcf = tmp.path().join("out.vcf");
    let mut cli = default_cli(bam.clone(), reference, vcf.clone());
    cli.method = AssayMethod::Amplicon;
    cli.regions = vec!["13:700-900".to_string()];
    run_cli(&cli).unwrap();
    let unindexed_vcf = fs::read(&vcf).unwrap();
    assert_eq!(parse_vcf_rows(&vcf).len(), 1);
    index_bam(&bam);
    run_cli(&cli).unwrap();
    assert_eq!(fs::read(&vcf).unwrap(), unindexed_vcf);
    fs::rename(bam.with_extension("bam.bai"), bam.with_extension("bai")).unwrap();
    run_cli(&cli).unwrap();
    assert_eq!(fs::read(&vcf).unwrap(), unindexed_vcf);

    let wrong_bam = write_synthetic_bam(tmp.path(), "wrong.bam", "14", wt.len(), &reads);
    cli.bam = wrong_bam.clone();
    let unindexed_error = run_cli(&cli).unwrap_err().to_string();
    index_bam(&wrong_bam);
    let indexed_error = run_cli(&cli).unwrap_err().to_string();
    assert!(unindexed_error.contains("not found in BAM header"));
    assert_eq!(unindexed_error, indexed_error);
}

#[test]
fn test_overlapping_regions_and_unindexed_fasta_with_description() {
    let tmp = TempDir::new().unwrap();
    let wt = revcomp(FLT3_E14E15_MINUS);
    let reference = tmp.path().join("ref.fa");
    fs::write(
        &reference,
        format!(
            ">13 example description\n{}\n",
            std::str::from_utf8(&wt).unwrap()
        ),
    )
    .unwrap();
    let reads = single_end_insertion_reads(&wt, &wt[750..759]);
    let bam = write_synthetic_bam(tmp.path(), "single.bam", "13", wt.len(), &reads);
    let vcf = tmp.path().join("out.vcf");
    let mut cli = default_cli(bam, reference, vcf.clone());
    cli.method = AssayMethod::Amplicon;
    cli.regions = vec!["13:1-1500".to_string(), "13:700-920".to_string()];
    run_cli(&cli).unwrap();
    let rows = parse_vcf_rows(&vcf);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].chrom, "13");
    assert_eq!(get_info_val(&rows[0].info, "SVLEN"), "9");
    assert!(
        fs::read_to_string(vcf)
            .unwrap()
            .contains("##contig=<ID=13,length=1500>")
    );
}
