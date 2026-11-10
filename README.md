# GGV: Graphic Genome Viewer

**A fast genome viewer for people and AI agents.** Look at reads, variants and genes in a
desktop app, ask the same questions from the command line, or let an agent ask them. Every
route gives the same numbers.

[![Latest release](https://img.shields.io/github/v/release/DylanCosto/ggv)](https://github.com/DylanCosto/ggv/releases/latest)
[![CI](https://github.com/DylanCosto/ggv/actions/workflows/ci.yml/badge.svg)](https://github.com/DylanCosto/ggv/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/github/license/DylanCosto/ggv)](https://github.com/DylanCosto/ggv/blob/main/LICENSE)

[Install](#install) · [Try the demo](#try-it-in-a-minute) · [Desktop app](#desktop-app) ·
[Command line](#command-line) · [AI agents](#ai-agents) · [Documentation](#documentation)

![GGV showing a heterozygous 6.8 kb deletion in HG002 long reads, grouped and coloured by haplotype: the reads on haplotype 1 carry the deletion and the reads on haplotype 2 span it intact](https://raw.githubusercontent.com/DylanCosto/ggv/main/docs/images/ggv-desktop.png)

*A heterozygous 6.8 kb deletion in PacBio HiFi reads of HG002, grouped by haplotype: the reads
that carry it are all on haplotype 1. This is the first stop of the built-in demo.*

## What you get

| Part | What it does |
|---|---|
| **Desktop app** | `ggv-desktop` opens BAM/CRAM, VCF/BCF, gene annotations and signal tracks. Frames take about a millisecond, so panning and zooming stay smooth on deep data. |
| **Command line** | `ggv` answers "what's at this base?" as TSV or JSON, draws PNG and SVG pictures without a display, and writes a report for every site in a list. |
| **AI agents** | `ggv mcp` serves the same answers over the Model Context Protocol. Attached to the desktop app, an agent sees what you are looking at. |

GGV 0.9 is a **beta for local files**. It shows and counts what is in your files; it is not a
variant caller.

## Try it in a minute

1. [Install GGV](#install).
2. Run `ggv-desktop --demo`, or start the app and press **Open the demo**.

The demo downloads 4.3 MB of real long-read data once (PacBio HiFi reads of GIAB HG002) and
walks you through six stops: a deletion on one haplotype, exact counts at a base, methylation,
variant review, the command palette and agents, then how to open your own files.

## Install

Download the package for your system from
[Releases](https://github.com/DylanCosto/ggv/releases/latest):

| System | Package | Needs |
|---|---|---|
| Windows 10/11 (x64) | `ggv-<version>-windows-x64.zip` | DirectX 12 or Vulkan |
| Linux (x86_64) | `ggv-<version>-linux-x86_64.tar.gz` | glibc 2.35+ (Ubuntu 22.04, Debian 12, Fedora 36 or newer), Vulkan or OpenGL |
| macOS 11+ (Apple silicon and Intel) | `ggv-<version>-macos-universal.zip` | nothing else |

Unpack it, then run the installer in that folder. Nothing needs administrator rights.

```sh
powershell -ExecutionPolicy Bypass -File install.ps1   # Windows
./install.sh                                            # Linux and macOS
```

The builds are unsigned, so there is a warning the first time. On Windows choose
**More info → Run anyway**; on macOS `install.sh` clears the quarantine for you.
[docs/install.md](https://github.com/DylanCosto/ggv/blob/main/docs/install.md) covers
uninstalling and checking downloads.

## Desktop app

```sh
ggv-desktop sample.bam chr17:43,044,806 --genome hg38     # built-in reference and genes
ggv-desktop tumor.bam normal.bam calls.vcf.gz genes.gtf.gz BRCA1 -f GRCh38.fa
```

Give it files in any order; a last argument that isn't a file is a locus or a gene name. Or
start it empty and drag files onto the window.

**Exact counts at any base.** Click a base and the Evidence panel counts every read there by
allele and strand, with the VCF record and what the change does to the gene.

![The Evidence panel at a heterozygous variant in GNAS: 10 reads with C and 15 with T, by strand, with the HGVS consequence c.393C>T p.Ile131=, and the gene track showing the amino acids](https://raw.githubusercontent.com/DylanCosto/ggv/main/docs/images/evidence.png)

**Methylation on each haplotype.** Colour by base modification and group by haplotype to see
imprinting directly: here one copy of GNAS is mostly unmethylated (blue) and the other mostly
methylated (red).

![Reads at the imprinted GNAS locus coloured by base modification and grouped by haplotype: haplotype 1 is mostly blue (unmethylated), haplotype 2 mostly red (methylated)](https://raw.githubusercontent.com/DylanCosto/ggv/main/docs/images/methylation.png)

**And more:**

- **Reads your way:** sort, group and colour by strand, haplotype, read group or any tag;
  filter with expressions such as `mapq >= 20 && hp == 1 && !dup`.
- **Review variants:** step through a VCF and mark each site real, artifact or unsure. Your
  decisions are saved beside the file as you go.
- **Two places at once:** split the window to see both ends of a rearrangement, or two genes,
  side by side.
- **Find anything:** Ctrl+K searches every command, gene, bookmark and recent file. `?` lists
  the keyboard shortcuts.
- **Whole chromosomes:** zoomed all the way out, depth is estimated from the index without
  reading the reads.

[docs/desktop.md](https://github.com/DylanCosto/ggv/blob/main/docs/desktop.md) has the rest:
splicing, structural variants, signal tracks, sessions and reports.

## Command line

```sh
ggv query pileup  sample.bam chr7:55,191,822 -f ref.fa         # one base
ggv query summary sample.bam chr7:55191000-55192000 -f ref.fa  # a region
ggv render sample.bam chr7:55,191,822 -f ref.fa -o snv.png     # a picture, no display needed
ggv batch sites.bed -b tumor.bam -b normal.bam -f ref.fa -v calls.vcf.gz -o review/
ggv index sample.bam                                           # build a missing index
```

`query pileup` at one base gives depth and alleles, here 10 reads with C and 15 with T (some
columns left out):

```text
#contig  pos       ref  depth  A  C   G  T   N  del  ins  ...  alt_fraction  top_alt
chr20    58903752  C    25     0  10  0  15  0  0    0    ...  0.6           T:15
```

`render --sort base` draws the same site as the app would:

![A headless render of the same heterozygous site: coverage with the alternate allele in colour, and 25 reads sorted so the 15 carrying T are on top](https://raw.githubusercontent.com/DylanCosto/ggv/main/docs/images/render.png)

`batch` writes a picture and an HTML evidence report for every site in a BED, VCF or list of
loci, with a contact sheet linking them. Add `--format json` to any query for scripts.
[docs/cli.md](https://github.com/DylanCosto/ggv/blob/main/docs/cli.md) has every option.

## AI agents

Add GGV to an MCP client with your files as the defaults:

```sh
claude mcp add ggv -- ggv mcp -b tumor.bam -b normal.bam -f GRCh38.fa -v calls.vcf.gz -g genes.gtf.gz
```

Or attach it to the desktop app, so calls use whatever is on screen and "is this variant
real?" needs no coordinates:

```sh
claude mcp add ggv -- ggv mcp --attach
```

For clients configured with JSON:

```json
{"mcpServers": {"ggv": {"command": "ggv", "args": ["mcp", "--attach"]}}}
```

Agents get numbers, not pictures. This is the `evidence` answer for the site above, abridged:

```json
{
  "locus": "chr20:58903752", "ref": "C", "depth": 25,
  "alleles": [
    {"allele": "C", "reads": 10, "forward": 4, "reverse": 6,  "fraction": 0.4, "reference": true},
    {"allele": "T", "reads": 15, "forward": 2, "reverse": 13, "fraction": 0.6}
  ],
  "genes": [{"gene": "GNAS", "transcript": "NM_000516.7", "hgvs_c": "c.393C>T",
             "hgvs_p": "p.Ile131=", "effect": "synonymous", "where": "exon 5/13"}],
  "warnings": []
}
```

| Tool | Answers |
|---|---|
| `evidence` | reads per allele and strand at a base, artifact warnings, the VCF record, the HGVS consequence |
| `reads`, `summarize`, `coverage` | reads grouped by allele, strand or haplotype; region statistics; depth |
| `variants`, `genes`, `consequence` | VCF records, transcripts, and what any change does to a gene |
| `phase`, `svs`, `repeats`, `junctions` | haplotypes from the reads, structural variants, repeat lengths, splicing |
| `compare`, `signal`, `render` | where two files differ (tumor against normal), signal tracks, and a picture when one helps |

Attached to the app, an agent can also move the view, change the display, add bookmarks and
overlays, and wait for you to click a base.
[docs/mcp.md](https://github.com/DylanCosto/ggv/blob/main/docs/mcp.md) lists every tool.

## What it reads

| Kind | Formats |
|---|---|
| Alignments | BAM and CRAM, indexed |
| Variants | VCF and BCF, bgzipped and indexed |
| Annotations | GTF, GFF3 and BED, plain or gzipped |
| References | indexed FASTA (plain or bgzipped) and 2bit |
| Signals | bigWig, bedGraph and SEG |

A file without an index is no obstacle: the app offers to build it, and `ggv index` does the
same from a terminal.

**Built-in genomes** save you finding a reference. `--genome hg38` downloads the sequence and
RefSeq genes once (877 MB for hg38) and keeps them. Also available: hg19, T2T (hs1), mm39,
mm10 and sacCer3, with ClinVar and other tracks for hg38 and hg19.
See [docs/genomes.md](https://github.com/DylanCosto/ggv/blob/main/docs/genomes.md).

**Not yet:** remote files (HTTPS, S3), cancelling a slow load, CRAM with a 2bit reference.

## Documentation

| Guide | Covers |
|---|---|
| [Desktop app](https://github.com/DylanCosto/ggv/blob/main/docs/desktop.md) | reading alignments, variants, review, splicing, signal tracks, keyboard shortcuts |
| [Command line](https://github.com/DylanCosto/ggv/blob/main/docs/cli.md) | queries, rendering and batch reports |
| [Agents (MCP)](https://github.com/DylanCosto/ggv/blob/main/docs/mcp.md) | every tool, and working attached to the desktop app |
| [Genomes](https://github.com/DylanCosto/ggv/blob/main/docs/genomes.md) | built-in references and tracks |
| [Install](https://github.com/DylanCosto/ggv/blob/main/docs/install.md) | packages, installers and first run |
| [Development](https://github.com/DylanCosto/ggv/blob/main/docs/development.md) | crates, correctness checks and building |

## Build from source

Rust 1.91 or newer:

```sh
cargo build --release -p ggv-cli -p ggv-desktop
```

The binaries are `target/release/ggv` and `target/release/ggv-desktop`. To package as a
release does, see `scripts/windows`, `scripts/linux` and `scripts/macos`. GGV reads files
through [pasta](https://github.com/DylanCosto/pasta), a pure-Rust library for
SAM/BAM/CRAM/VCF/FASTA.

## Feedback

Found a bug? [Report it](https://github.com/DylanCosto/ggv/issues/new?template=bug_report.yml)
with the diagnostics from **Commands → About GGV**.

## License

MIT; see [LICENSE](https://github.com/DylanCosto/ggv/blob/main/LICENSE). Bundled open-source
components are listed in
[THIRD_PARTY_NOTICES.md](https://github.com/DylanCosto/ggv/blob/main/THIRD_PARTY_NOTICES.md),
including the MANE transcript list (NCBI/EMBL-EBI, public domain).
