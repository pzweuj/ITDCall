# itdcall

[English (primary README)](README.md)

纯 Rust 编写、零外部依赖的二代测序（NGS）内部串联重复（Internal Tandem Duplication, ITD）检测软件。  
核心算法参考 [FLT3_ITD_ext](https://github.com/ht50/FLT3_ITD_ext) 的局部比对—序列延伸—聚类—重比对定量思路，采用纯 Rust 实现（内置 BAM/FASTA 解析、k-mer 种子 + 带状仿射空位 Smith-Waterman 局部比对与软剪切锚定、贪心编辑距离聚类及方法学自适应定量），无需 `bwa` / `samtools` / `sumaclust` / `bedtools` 等外部工具即可单二进制跨平台运行。

## 核心特性

- **双方法学原生适配**：
  - **探针捕获（`--method hc` / `capture`）**：针对随机打断片段特点，采用左右断点边界平均法（`VD = (MutLeft + MutRight) / 2`，`DP = (WTLeft + WTRight + MutWTLeft + MutWTRight + OtherMut) / 2`），即使对于超过读长的超长 ITD（如 180–500 bp）也能给出无偏的 VAF 与 Allelic Ratio（AR）。
  - **扩增子测序（`--method amplicon`）**：针对固定引物起止位置特点，采用片段级分子定量（`VD = Mut`，`DP = WT + Mut + OtherMut`），同时默认保留 `0x400` 标记以免扩增子读段被误过滤。
- **输入输出简洁规范**：
  - 输入：BAM 文件 + 参考基因组 FASTA。`targeted` 模式可用 `.bai` 加速，索引查询与流式扫描均提取比对区间与目标重叠的读段。`loose` 模式通过流式扫描额外保留未比对读段。参考 FASTA 有 `.fai` 时可随机访问，无索引时自动回退。
  - 输出：标准 VCFv4.2（不含基因/转录本/HGVS 注释，仅输出左对齐标准化后的基因组 `+` 链 `CHROM`/`POS`/`REF`/`ALT` 及 `SVLEN`/`DUPSTART`/`DUPEND`/`DUPLEN`/`INSLEN`/`DP`/`VD`/`AF`/`AR`/`CR` 等定量与质控指标）。
- **内置经典 FLT3 区域 + 自定义区域校正**：
  - 内置 `hg19`（`chr13:28607438-28608937`）与 `hg38`（`chr13:28033301-28034800`）经典 FLT3 外显子 14–15 区域，默认 `--genome auto` 自动识别版本与 `chr13`/`13` 命名。
  - 支持通过 `--target-fasta` 传入用户自定义的局部目标区域 FASTA 文件（类似 `FLT3_dna_e14e15.fa`）。程序在输入参考基因组中定位目标并识别方向，从 `--reference` 提取对应的 `+` 链序列用于比对和生成 VCF。
  - 同时支持 `--region chr:start-end`（可重复指定）或 `--bed regions.bed` 扩展至任意基因/区域。
- **测序数据兼容**：以双端测序（Paired-end）为主，同时兼容单端（Single-end）及孤立（Orphan）reads；内置 Illumina 3' 接头裁剪与多克隆（Multi-clone）共存在定量。

## 构建

```bash
cargo build --release
```

编译产物位于 `target/release/itdcall`（Windows 为 `target\release\itdcall.exe`）。

## 使用示例

### 1. 探针捕获（Hybrid Capture），使用内置 FLT3 区域（自动识别 hg19 / hg38）

```bash
itdcall \
  -b sample.bam \
  -r hg38.fa \
  -o sample.itd.vcf \
  -m hc
```

### 2. 扩增子测序（Amplicon），使用内置 FLT3 区域

```bash
itdcall \
  -b amplicon.bam \
  -r hg19.fa \
  -o amplicon.itd.vcf \
  -m amplicon \
  --genome hg19
```

### 3. 用户提供局部目标 FASTA 进行区域校正（如 `FLT3_dna_e14e15.fa`）

```bash
itdcall \
  -b sample.bam \
  -r hg19.fa \
  -o sample.itd.vcf \
  -m hc \
  --target-fasta FLT3_dna_e14e15.fa
```

若 FASTA 位于特殊自定义区域，也可附加 `--region chr13:28607438-28608937` 加速定位。

### 4. 自定义基因组区域或 BED 多区域检测

```bash
# 单区域/多区域
itdcall -b sample.bam -r hg38.fa -o out.vcf -m hc --region chr13:28033301-28034800

# BED 批量区域
itdcall -b sample.bam -r hg38.fa -o out.vcf -m hc --bed targets.bed
```

### 5. 常用调优参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-m, --method` | `hc` | 定量方法学：`hc`（探针捕获）或 `amplicon`（扩增子） |
| `-g, --genome` | `auto` | 内置 FLT3 参考版本：`auto` / `hg19` / `hg38` |
| `--target-fasta` | 无 | 用户自定义局部参考 FASTA（类似 `FLT3_dna_e14e15.fa`） |
| `--region` / `--bed` | 无 | 自定义基因组坐标区间（1-based 闭区间 / BED） |
| `--read-mode` | `targeted` | BAM 读段提取模式：`targeted` / `loose`（含 unmapped）/ `all` |
| `--min-reads` | `2` | `PASS` 最小有效支持深度（VD） |
| `--min-vaf` | `0.0` | `PASS` 最小变异等位基因频率（AF） |
| `--min-itd-len` / `--max-itd-len` | `3` / `500` | 检测 ITD 净插入长度范围（bp）。候选默认还需至少 3 bp 的重复片段，若下限设为 1 或 2，则采用该值作为最短重复长度 |
| `--in-frame-only` | `false` | 仅保留 3 bp 整数倍（in-frame）ITD |
| `--include-filtered` | `false` | 在 VCF 中输出非 `PASS` 过滤候选 |
| `--no-trim-adapters` | `false` | 关闭内置 Illumina 3' 接头裁剪 |
| `--exclude-duplicates` | `false` | 排除 BAM `0x400` 重复标记（扩增子请保持关闭） |
| `-t, --threads` | `0`（自动） | 并行线程数 |

## 输出 VCF 字段说明（无注释）

- `INFO`:
  - `SVTYPE=ITD`：变异类型
  - `SVLEN`：ITD 净插入长度（bp）
  - `DUPSTART` / `DUPEND` / `DUPLEN`：重复片段在参考基因组 `+` 链上的起止坐标及长度
  - `INSLEN`：断点处非模板插入序列长度（bp）
  - `METHOD`：定量方法学（`hc` / `amplicon`）
  - `LOCUS`：目标区域标识
  - `DP` / `VD` / `AF` / `AR`：有效总深度、变异支持深度、等位基因频率、等位基因比值。`AR` 用 `VD` 除以实际 WT 深度（探针捕获取两侧 WT 深度均值，扩增子取 WT 片段数），无 WT 支持时输出 `Inf`。
- `FORMAT`: `GT:DP:VD:AF:AR:CR`（`CR` 为重复区间内部覆盖半径分布直方图 `0,1-5,6-10,11-25,26-50,>50`）

重叠目标检出同一 `CHROM/POS/REF/ALT` 时，仅保留一条记录，依次优先选择 `PASS`、较高 `VD`、较高 `DP` 的记录。

## 测试

```bash
cargo test
```

## 发布

将提交标题写成独立的版本号（如 `v0.0.1`）并推送到 GitHub 分支。[发布工作流](.github/workflows/release.yml)会自动测试、编译 Linux x86-64 musl 版本，并为匹配的提交创建同名标签和 GitHub Release。若需让二进制的 `--version` 与标签一致，请在提交前同步更新 Cargo 包版本及 `Cargo.lock`。

## 许可证

本项目采用 [MIT 许可证](LICENSE)，其中保留了 [FLT3_ITD_ext 的原始版权声明](https://github.com/ht50/FLT3_ITD_ext/blob/master/LICENSE)。
