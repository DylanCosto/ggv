# Development

## Crates


| Crate | Role |
|---|---|
| `ggv-model` | Loci, regions, contigs. No I/O, so it builds for WASM. |
| `ggv-io` | Indexed BAM, CRAM (with its reference), FASTA, and VCF/BCF (TBI/CSI) through pasta; locus resolution (`chr7` ⇄ `7`, `chrM` ⇄ `MT`). |
| `ggv-query` | `pileup` (per position) and `summarize` (per region). |
| `ggv-annot` | GTF, GFF3 and BED gene annotations (plain or gzipped): transcripts with an interval index and gene-name search. |
| `ggv-layout` | Decodes reads (blocks, indels, mismatches vs reference), packs rows, bins coverage. |
| `ggv-render` | Builds a backend-neutral `Scene`; rasterizes it to PNG (tiny-skia, built-in 5×7 font) or writes SVG. |
| `ggv-cli` | The `ggv` binary. |
| `ggv-desktop` | The desktop app (egui + wgpu). |

The MCP server is `ggv mcp` (see [Agents (MCP)](mcp.md)). Planned next: remote files (HTTPS).

## Correctness

`tests/data/hg002-revio` is a 4 MB slice of real PacBio Revio HiFi reads of HG002 kept in the
repository, so the long-read checks run everywhere, CI included
(`crates/ggv-layout/tests/long_reads.rs`): a heterozygous 6.8 kb deletion is carried only by HP 1
reads and spanned intact only by HP 2 reads; every supplementary record lists its other parts;
GGV's allele phasing agrees with PacBio's haplotagging on every read both phase; at least 97% of
5mC calls land on CpGs on both strands; and at the imprinted GNAS A/B region methylation splits
by haplotype (about 17% against 88%).

`scripts/compare_mpileup.py` checks `ggv query pileup` against `samtools mpileup -A -B -Q 0 -q 0
-d 0 -x --reverse-del` at every position: depth, per-strand base counts, per-strand deletions and
insertions, ref skips, and reference base. On the HG002 chr20 fixtures (Illumina 60×, HiFi 48×)
all 22,913 and 374,967 positions match exactly.

`crates/ggv-layout/tests/coverage_parity.rs` checks that coverage derived from decoded reads
(what `ggv render` and the app draw) equals the pileup engine at every position: depth,
deletions, allele counts, alt fraction. It runs when `GGV_TEST_BAM`, `GGV_TEST_REGIONS` (and
optionally `GGV_TEST_REF`) are set. It passes on the HG002 fixtures and on 27,008 positions of
deep NA12878 BRCA1 data. With `GGV_TEST_WGS_BAM` and `GGV_TEST_WGS_REGIONS` it also checks the
index depth estimate against counted depth (mean within 25%, bins correlated at r > 0.8).

Note `-A`: by default samtools mpileup silently drops anomalous pairs (paired, not properly
paired). GGV, like IGV, shows them.

## Building on Windows

The checkout path here is long enough to hit the 260-character limit when linking, so build with
a short target dir if you see `LNK1104`:

```sh
CARGO_TARGET_DIR=C:/t cargo build --release
```

pasta (0.9.0) comes from its public GitHub repo at a pinned commit; `.cargo/config.toml` has
cargo fetch it with the git CLI.
