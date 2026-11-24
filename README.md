# itdcall

[简体中文](README.zh-CN.md)

`itdcall` detects internal tandem duplications (ITDs) in next-generation sequencing data. It reads an aligned BAM and a reference FASTA and writes a VCF. It runs as a single Rust executable without external tools such as BWA or SAMtools.

The caller uses local alignment, read extension, allele clustering, and realignment for quantification. It supports paired-end, single-end, and orphan reads, trims Illumina 3′ adapters by default, and quantifies samples with multiple ITD clones. The approach was inspired by [FLT3_ITD_ext](https://github.com/ht50/FLT3_ITD_ext).

## Build and install

GitHub Releases provides a statically linked Linux x86-64 executable. Extract the archive and run `./itdcall --help`.

To build from source with Rust and Cargo:

```sh
cargo build --release --locked
```

The executable is at `target/release/itdcall` on Linux and macOS, or `target/release/itdcall.exe` on Windows.

## Usage

Hybrid capture with the built-in FLT3 target:

```sh
itdcall -b sample.bam -r hg38.fa -o sample.itd.vcf -m hc
```

Amplicon sequencing with an explicit reference build:

```sh
itdcall -b amplicon.bam -r hg19.fa -o amplicon.itd.vcf -m amplicon --genome hg19
```

Custom targets:

```sh
itdcall -b sample.bam -r hg19.fa -o out.vcf --target-fasta FLT3_dna_e14e15.fa
itdcall -b sample.bam -r hg38.fa -o out.vcf --region chr13:28033301-28034800
itdcall -b sample.bam -r hg38.fa -o out.vcf --bed targets.bed
```

- The built-in target covers FLT3 exons 14–15 at `chr13:28607438-28608937` on hg19 and `chr13:28033301-28034800` on hg38. The default `--genome auto` detects the build and accepts both `chr13` and `13`.
- `--target-fasta` locates the supplied sequence in the reference on either strand and uses the corresponding genomic sequence for calling. With a target FASTA, `--region` can provide a search hint.
- `--region` accepts repeatable 1-based inclusive intervals. BED intervals are 0-based and half-open.
- A `.bai` accelerates BAM extraction in `targeted` mode. Without it, extraction uses the same overlap rule in a streaming scan. `loose` also scans for unmapped reads, and `all` considers the entire BAM. A `.fai` enables random access to the reference FASTA.

## Common options

| Option | Default | Description |
| --- | --- | --- |
| `-m, --method` | `hc` | Hybrid capture (`hc` or `capture`) or amplicon sequencing (`amplicon`). |
| `-g, --genome` | `auto` | Built-in FLT3 reference build, `auto`, `hg19`, or `hg38`. |
| `--read-mode` | `targeted` | BAM extraction mode, `targeted`, `loose`, or `all`. |
| `--min-reads` | `2` | Minimum effective variant depth (VD) for PASS. |
| `--min-vaf` | `0.0` | Minimum variant allele fraction (AF) for PASS. |
| `--min-itd-len` / `--max-itd-len` | `3` / `500` | Net insertion length range in bp. The minimum duplicated segment is 3 bp by default and follows the minimum ITD length when set to 1 or 2. |
| `--in-frame-only` | Off | Require the net insertion length to be a multiple of three for PASS. |
| `--include-filtered` | Off | Include filtered calls in the VCF. |
| `--no-trim-adapters` | Off | Disable Illumina adapter trimming. |
| `--exclude-duplicates` | Off | Exclude BAM reads marked as duplicates (`0x400`). Leave off for amplicon data. |
| `-t, --threads` | `0` (automatic) | Number of worker threads. |

## Quantification and output

- **Hybrid capture** averages support across the left and right breakpoints. VD is `(MutLeft + MutRight) / 2`. DP also uses boundary counts from wild-type and mutant fragments.
- **Amplicon** counts fragments. VD is the number supporting the current allele. DP includes wild-type fragments and fragments supporting any detected clone at the breakpoints.
- **AF** is `VD / DP`. **AR** is `VD / WT` and is reported as `Inf` when WT depth is zero.

The VCFv4.2 output uses left-aligned alleles in genomic plus-strand coordinates. It includes `SVTYPE=ITD`, `SVLEN` (net insertion length), `DUPSTART`, `DUPEND`, `DUPLEN`, `INSLEN` (non-templated insertion length), `METHOD`, and `LOCUS`. The sample fields are `GT:DP:VD:AF:AR:CR`. `CR` reports coverage radius in six bins (0, 1–5, 6–10, 11–25, 26–50, and over 50 bp). Gene, transcript, and HGVS annotations are outside the output.

If overlapping target intervals yield the same `CHROM/POS/REF/ALT`, the caller keeps one record, prioritizing PASS, higher VD, and then higher DP.

## Tests and releases

```sh
cargo test --locked
```

To release, push a commit whose entire subject is a version such as `v0.0.1`:

```sh
git commit -m "v0.0.1"
git push
```

The [release workflow](.github/workflows/release.yml) finds version commits in a branch push, runs the tests, and publishes a Linux x86-64 musl archive to GitHub Releases. It tags the commit with its subject. Set the Cargo package version and update `Cargo.lock` before committing if the executable's `--version` should match the release tag.

## License

Licensed under the [MIT License](LICENSE). The license preserves the copyright notice from [FLT3_ITD_ext](https://github.com/ht50/FLT3_ITD_ext/blob/master/LICENSE).
