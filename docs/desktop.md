# The desktop app

```sh
ggv-desktop sample.bam chr17:41,244,435 -f ref.fa              # one track
ggv-desktop tumor.bam normal.bam chr17:43,044,806 -f ref.fa    # one track per BAM
ggv-desktop *.bam --locus chr17:43,090,000-43,100,000 -f ref.fa
ggv-desktop a.bam b.bam BRCA1 -g Homo_sapiens.GRCh38.112.gtf.gz -f ref.fa   # jump to a gene
ggv-desktop -g genes.gtf.gz TP53                                           # genes only, no BAM
ggv-desktop child.bam trio.vcf.gz giab.vcf.gz genes.gtf.gz chr20:1,895,892 -f ref.fa
```

CRAM works anywhere a BAM does (`.crai` index next to it); pass the reference it was compressed
against with `-f`. On the HG002 fixtures, `ggv query pileup` over CRAM is byte-identical to the
BAM across all of chr20.

Positional files are sorted by extension: `.bam`/`.cram` are alignment tracks, `.vcf(.gz)`/`.bcf`
variant tracks, `.gtf`/`.gff3`/`.bed` (optionally `.gz`) gene tracks; a trailing argument that is
none of these is the locus.

Tracks stack below one shared ruler and reference row and share the locus and zoom. Each has
its own coverage, read rows and vertical scroll. Add tracks with **Add BAM…** (several at once)
and close one with the × in its header. BAMs may name contigs differently (`chr17` vs `17`);
each track maps the view's contig through the same aliases as the CLI, and says so in its
header when a file lacks the contig.

**Reading alignments** (toolbar, or right-click a base for IGV's context menu):

- **Sort** the rows at a base by base (mismatches A/C/G/T/N first, then deletions, then
  reference), strand, MAPQ, start, or insert size. As in IGV, rows are reordered by the read
  over that base, so reads stay packed. `--sort base` sorts at the locus on startup.
- **Group** by strand, haplotype (`HP` tag, for phased long reads), read group (`RG`) or **any
  tag** (type it in the Group or Colour menu: `CB` cell barcodes, `PS` phase sets, `UB`, `MI`, an
  aligner's own), each group in its own labelled band. **Colour by tag** gives each value its
  own colour, the same in the app and in `ggv render` (`--group tag --tag CB`, `--color-by tag`;
  MCP `group_by` / `color_by` `tag:CB`). Reads are reloaded with the tag when it changes.
- **Link split reads** (checkbox, `--link`, MCP `link`): each split read's primary and
  supplementary records go on one row, joined by a thin blue line, in the primary's group (as
  IGV's "link supplementary alignments"); mates of a pair are never joined.
- **Colour** by strand, haplotype, or insert size and pair orientation: inserts above the
  99.5th / below the 0.5th percentile of the loaded pairs red / blue, mate on another contig
  orange, outward (RL, for mates that don't overlap) green, same-strand LL teal / RR dark blue.
- **Base modifications** (`MM`/`ML` tags, e.g. 5mC/5hmC/6mA from ONT or PacBio): IGV's
  two-colour scale, red where a base is likely modified and blue where it is likely not,
  stronger the more certain (codes at one base are summed). Hover shows the likelihood at the
  base and each read's calls per code. Calls are decoded with pasta's port of HTSlib's MM/ML
  parser; on the ONT fixture they match `samtools mpileup --output-mods` exactly on all
  106,437 (position, code) entries over 100 kb (counts and summed likelihoods; calls on
  inserted bases, which samtools also prints, aren't drawn).
- **Soft clips**: show clipped bases where they would align, paler, with their mismatches.
- Copy a read name or position from the context menu. Read tooltips add HP, RG, insert size
  and mate position.
- Startup flags: `--color` (`strand`, `haplotype`, `insert-size`, `base-mods`), `--group`,
  `--sort`, `--soft-clips`.
- Right-click a track for **Move track up / down** and **Close track**.

**Sessions** (**Session** menu): **Save as…** writes the open files, view, outlined base,
display options and every split pane (with whether they are linked) to a small JSON file;
**Open…** restores one; the window is also saved every 20 seconds while it changes and on exit,
and comes back with **Restore last session** or `ggv-desktop --restore` (`--session FILE` opens a
saved one). Saving never leaves a half-written file: each save goes to a temporary file that
replaces the old one only once complete, and the previous version is kept as `<name>.bak`, which
GGV opens instead (and says so) if the file is ever damaged. Files inside the session's folder
are stored relative to it (with `/`), so a project folder holding the session and its data can
be moved, copied to another machine or shared. A file that isn't where the session says opens a
**Files not found** window: **Locate…** one and GGV looks for the others in the same folder (and
the folders in it), or **Look in a folder…**; the rest of the session opens meanwhile. Options
are stored by name, unknown fields are ignored, so session files stay readable and editable by
hand:

```json
{
  "version": 1,
  "bams": ["hifi.bam", "illumina.bam"],
  "vcfs": ["giab.vcf.gz"],
  "genes": ["gencode.gtf"],
  "reference": "GRCh38.fa",
  "locus": "chr20:1,895,842-1,895,942",
  "highlight": "chr20:1895892",
  "color": "strand",
  "group": "strand",
  "sort": {"pos": 1895892, "by": "base"}
}
```

Sorting and grouping happen on the loader thread with each load, and re-pack from memory
(no disk) when changed; frame times are unchanged (three exome tracks grouped and coloured by
insert size: 0.5–0.7 ms median, under 4 ms worst).

**Variant tracks** (`-v`/`--vcf`, positional `.vcf.gz`/`.bcf`, or **VCF…**) read bgzipped VCF
with a TBI or CSI index, or BCF with CSI, through pasta's HTSlib-compatible engine. They sit at
the top, as in IGV: a site row (red for the alt allele fraction, from INFO/AF or the genotypes;
faded when FILTER is not PASS) and a genotype row per sample (grey hom-ref, blue het, cyan
hom-alt, pale no-call). Hover a site for position, type, REF/ALT, QUAL, FILTER, AF, genotype
counts and INFO; hover a sample row for its GT, DP, GQ and AD. gVCF reference blocks
(`<*>`/`<NON_REF>`) are skipped. **< Variant** / **Variant >** (or Ctrl+B / Ctrl+F, over the
track to use) jump to the previous/next site from the outlined base, searching beyond what is
loaded and on into the next (or previous) chromosomes, outlining the site and zooming to read
level if the view is wide.

**Gene tracks** (`-g`/`--genes`, or **Genes…**) read GTF, GFF3 or BED files, gzipped or not. They
sit at the bottom, as in IGV: thick CDS, thin UTRs, introns with strand chevrons, and names
underneath (skipped where they would collide). One transcript per gene is shown (MANE Select,
then MANE Plus Clinical, then Ensembl canonical, then the curated transcript with the longest CDS,
then predicted models); MANE comes from the file's tags or, for files without them such as UCSC's
ncbiRefSeq, from the MANE GRCh38 v1.5 list (NCBI/EMBL-EBI, public domain) built into GGV. Tick **All transcripts** to see every one. Hover
for the transcript id, biotype, exon count, coding length, tags, and which exon (numbered in the
direction of transcription) is under the pointer. Double-click a gene to fit it in the view. Type a gene name, gene id or transcript id
(`BRCA1`, `ENSG00000012048`, `ENST00000357654.9`) in the locus box to jump to it. With no BAM
open, an annotation alone defines the genome.

Annotations over 4 MB are parsed once and cached (in `~/Library/Caches/ggv/annotations`,
`$XDG_CACHE_HOME/ggv/annotations` or `%LOCALAPPDATA%\ggv\annotations`; `GGV_ANNOTATION_CACHE` names
another directory, or `off`). Later opens, in the app, `ggv batch` and `ggv mcp`, read the cache:
GENCODE v50 (644,292 transcripts) opens in about 10 ms and 24 MB instead of 3.6 s and 670 MB,
and a contig's transcripts are decoded the first time it is viewed (chr1 in about 17 ms). A cache
is keyed by the file's path, size and modification time, so an edited file is parsed again; delete
the directory to reclaim the space.

Ensembl 112's full GRCh38 GTF (55 MB gzipped, 1.5 GB of text; 254,129 transcripts, 63,140
genes) parses in about 1.9 s on a background thread: decompression runs on its own thread, and
exon/CDS lines of known transcripts take a fast path that reads only the transcript id. Frame
times with it loaded: about 1 ms across all of chr17, 0.4 ms at 2 Mb, and unchanged (about
1 ms) with three deep exome BAMs at 10 kb.

**Working faster:**

- **Open files** by dropping them on the window (BAM/CRAM, VCF/BCF, annotations, a FASTA as the
  reference, or a session `.json`), from **Session → Recent files**, or the toolbar buttons. A
  file without an index gets an offer to build one (BAI, falling back to CSI for long contigs;
  CRAI; TBI; FAI) on a background thread, after which it opens; rebuilt indexes give results
  identical to samtools-made ones.
- **Back / Forward** (`<` `>` or Alt+Left / Alt+Right) step through the places you jumped to
  (locus box, gene names, variant and mate jumps).
- **Filters** (hide duplicates, secondary, supplementary, QC-fail; minimum MAPQ) apply live from
  memory to reads, coverage and counts, and match what loading with the same filters gives.
- **Coverage** scale: autoscale to the view, to the loaded data, or log; optionally the same
  scale for every track.
- **Session → Export view as PNG…** saves the tracks exactly as drawn. **Dark** switches theme.
- **Group by base here** (right-click a base, or `--group base`): reads in bands by the base they
  carry there (A/C/G/T, deletion, not covering), for eyeballing which reads support an allele.
- **Evidence** panel: exact counts at the outlined base (click a base to outline it, click it
  again or Esc to clear; once it is panned off screen the panel follows the view's centre again)
  from the same pileup engine as `ggv query pileup`: per-track depth, bases by strand, deletions,
  insertions, MAPQ and base-quality summaries, soft clips, the top alt with a warning when it is on
  one strand only while the other reads are on both (at a capture edge every read can be on one
  strand, and the panel says strand bias can't be judged there); **Copy** buttons for the locus
  and for each file's counts (as a small table), and a click on a file's name copies it with its
  full sample name; plus the genes (exon number, coding or not) and every VCF record there with
  genotypes, DP and AD.
- **Variant triage** (**Review…**, or `--review`): steps through a VCF's sites. Press 1 / 2 / 3
  to mark the outlined site real / artifact / unsure (with an optional note) and move on, N / P
  to skip. Decisions are saved on every change to `<vcf>.ggv-review.tsv` next to the VCF, so a
  review resumes where it stopped; reviewed sites get a coloured tick on the track. When the
  VCF's folder can't be written (a read-only share), decisions go to GGV's own folder instead and
  the review says where; **Save reviews in…** picks a folder for every review and queue. Saves
  are atomic with a backup, like sessions, and a save that fails stays on screen until one works.
- **Assembly checks**: the status bar names the open genome's build and whether a file declares
  it (a VCF's `##reference`, a file name like `hg38.ncbiRefSeq.gtf.gz`) or GGV inferred it from
  chromosome 1's length, which is a hint rather than a verification. GGV compares the contig
  lengths every file states (reference, BAM/CRAM headers, VCF `##contig` lines, bigWig) and
  checks VCF REF alleles in view against the reference; when files disagree (an hg19 BAM with
  hg38 genes, a VCF called on another reference) the label becomes **Assembly mismatch**, and
  **View → Assembly and reference…** explains each file's build and every disagreement before
  mismatch colours or consequence labels are trusted. Custom assemblies are simply "unknown".
- **Mates and split reads**: right-click a read for **Go to mate** and **Go to split part**
  (from `SA`), which jump there (across chromosomes) and outline it; tooltips list both.
- **Region statistics**: Shift-drag across the tracks to select a region; a window gives each
  track's read counts by filter, mean/median/min/max depth, breadth at 1/10/20/30×, and the
  variable sites, matching `ggv query summary`, with **Copy as JSON**.
  `--select chr20:1,890,000-1,900,000` opens with one.
- **Split view**: **Split** opens another pane beside the current one, with the same files and
  display options but its own position, zoom, history and outlined base (up to 4). Click a pane
  to make it active: the toolbar, locus box and keys act on it. Files opened and options changed
  in any pane apply to all; **Close** in a pane's header removes it. `--split LOCUS` (repeatable)
  opens panes on startup, e.g. both ends of a translocation:
  `ggv-desktop tumor.bam chr9:130,700,000-130,720,000 --split chr22:23,280,000-23,300,000`.
- **Whole chromosomes**: views wider than 5 Mb show depth estimated from the index (BAI/CSI
  offsets, or CRAI slice sizes) instead of nothing, in a lighter blue-grey and labelled as an
  estimate. The index gives read density per 16 kb window without decoding a read; each part of
  the view (up to 64) is then calibrated in three density strata (its sparsest, typical and
  densest windows, a third of its data each), each on its own sample counted exactly with the
  read filters and sized by data (about 8 MB of alignments in all), so dense windows full of
  repeats or duplicates don't inherit the typical window's ratio. On two whole-genome BAMs, in
  1 Mb bins: unfiltered, the median bin is within 1.1–1.7% of counted depth and nine in ten
  within 5.5%, in about 0.1 s for 20–40 Mb. With a MAPQ filter the median bin is within 1.6%
  and nine in ten within 10%, but single bins in repeats can still be off by a fifth to a half:
  the index can't see MAPQ. Zoom under 5 Mb (or use `ggv query summary`) for counted depth.

**Beyond IGV** — everything below is also in the **command palette** (Ctrl+K), and `?` lists
every shortcut:

- **What a change does**: with a gene annotation and a reference, the evidence panel, variant
  tooltips and reports give each change's consequence for every transcript there, in HGVS:
  `CDC25B ENST00000245960.10  c.42C>T  p.Leu14=  synonymous  (exon 1/16)`, missense, stop
  gained/lost, start lost, frameshift, inframe insertion/deletion, splice donor/acceptor/region,
  intron (`c.436+143A>G`), UTR (`c.-14`, `c.*20`). Indels are moved to their most 3′ position in
  a repeat, as HGVS requires. On the HG002 fixtures all seven coding variants checked against
  Ensembl VEP agree in `c.`, `p.` and consequence, on both strands, and every complete coding
  sequence in GENCODE v46 and Ensembl 112 translates from start to stop with no internal stop.
  Hovering a gene names the base's `c.` position, codon and amino acid.
- **Amino acids** on the gene track when zoomed in (a codon at least 9 points wide): one letter
  per codon, alternate codons shaded, the start codon green and stops red.
- **Group by alleles** (Group menu, or `--group alleles`): splits reads into two haplotypes by the
  alleles they carry at the heterozygous sites in the loaded span, for data without `HP` tags.
  Sites are found from the reads (one alternate base at 25–75%, depth ≥ 6), oriented left to
  right through reads that span two or more, and each read goes to the haplotype its alleles
  vote for; the group labels give the site and block counts. Against phase from transmission in
  the HG002 trio (which uses no reads) it agrees on 100% of reads both place (HiFi and Illumina),
  in under 1 ms per load. Needs a reference (sites come from mismatches).
- **Command palette** (Ctrl+K, or **Commands**): one fuzzy box for genes, loci, contigs,
  bookmarks, recent files, saved filters and every action. Matching all 63,140 Ensembl genes
  takes about 1.5 ms per keystroke; results are re-ranked only when the query changes.
- **Review queue** (**Tools → Review queue**, or `--queue sites.bed`): step through any list of
  places (BED, VCF, `chr1:100-200` lines, or text pasted into **Review queue: paste loci…**).
  N / P step, 1 / 2 / 3 mark pass / fail / unsure (with a note) and move to the next undecided
  item, U jumps to the next undecided one. Decisions save on every change to
  `<list>.ggv-queue.tsv` (pasted lists: GGV's config folder), so a review resumes where it
  stopped, and export as TSV or BED.
- **Bookmarks** (Ctrl+D bookmarks the selection, the outlined base or the view; Ctrl+Shift+B
  lists them): notes, jump, delete, **Export BED…**. They show as purple bands on the ruler,
  are shared by split panes and saved in sessions.
- **Linked panes**: tick **Link panes** (or `--link-panes`) and every pane pans and zooms with the
  active one, keeping its offset (both ends of a translocation move together).
- **Reads load ahead**: each track's loader keeps its recent read-level spans and, while idle,
  decodes the spans on either side of the view (the direction of travel first). In the fling
  benchmark at 10 kb views, 51 of 79 loads were served from memory instead of the disk; at
  1 kb, all of them. Frames are unchanged (0.3–0.4 ms median).
- **Files reload when they change** (on by default; **Tools**): a background thread checks each
  open file and its index once a second and reloads one in place, keeping the view, once it has
  stopped changing (a pipeline re-run). Idle, it costs a few `stat` calls a second.
- **Filter expressions** (Ctrl+E, or `--filter`): `mapq >= 20 && hp == 1 && !dup`,
  `flag & 0x800`, `name ~ "ccs"`, `insert > 1000 || mate_elsewhere`, `tag(HP) == 2`. Fields cover
  MAPQ, flags, positions, strand, pairing, HP/RG, name, mismatch/indel/clip counts, base mods
  and split reads; the window lists them and checks the expression as you type. Apply to every
  track, or right-click a track → **Filter this track…** for one file (remembered for that file
  the next time it opens). Save filters for reuse; they appear in the palette. Expressions apply
  where reads are loaded (views within the read window) and run on the loader thread with no allocation.
- **Find a read by name** (Ctrl+Shift+F, or `--find-read`): loaded reads first (instant), then
  every open file, contig by contig on all cores, stopping at the first hit; it jumps there and
  outlines the read, with buttons for its mate and split parts and **Copy as SAM**.
- **Find a sequence** (Ctrl+Shift+S, or `--find-seq PRIMER/2`): a primer or oligo (4–64 bases,
  IUPAC codes allowed) on both strands with up to 3 mismatches, in the view, the contig or the
  whole reference, offline. Matching is bit-parallel and split across cores: a 500 Mb genome in
  0.2 s exact and 0.3 s with 3 mismatches (from the OS cache). Hits are marked above the
  reference row; click one to go there.
- **Compare two tracks** (tumor/normal; **Tools → Compare**, or `--compare`): positions where any
  allele's fraction (A/C/G/T/N/deletion) differs by more than a threshold (default 20%; at
  100%, a full 100% difference) with enough depth in both (default 10) are marked orange through both tracks, from what is already
  loaded (no disk). `]` / `[` step to the next / previous one; beyond what is loaded a
  background pileup scan of both files finds it, across chromosomes.
- **Measure** (Ctrl+drag): length, GC content, and each track's mean depth and read count over
  the span, live while dragging, from loaded data.
- **Copy as**: right-click to copy a read's sequence or SAM line, the selection (or view) as
  FASTA, the selected reads as SAM, or the view as an image (Ctrl+Shift+C) straight to the
  clipboard. Copies that need the disk run on a background thread.
- **Evidence report** (Ctrl+Shift+R, or `--report out.html`): one self-contained HTML file for
  the outlined base: the view as an embedded PNG, the VCF records with genotypes, each track's
  allele counts by strand with MAPQ and base-quality summaries (and a warning for one-strand
  support), the selected region's summary, and review, queue and bookmark notes. Light and dark
  aware; prints to PDF.
- **Deep links**: **Session → Copy link to this view** gives
  `ggv://open?file=…&file=…&ref=…&locus=…&highlight=…`; `ggv-desktop "ggv://open?…"` opens it (each
  file is its own `file=` value, so paths may contain commas). `register-file-types.ps1` and
  `install.sh` register the `ggv://` scheme so links open GGV from a browser, chat or notebook.
  **Copy command that reopens this view** gives the equivalent `ggv-desktop …` line.

| Do | How |
|---|---|
| Anything | Ctrl+K, type a few letters |
| Go to a locus | Type in the locus box (Ctrl+L): `chr17`, `chr17:41,244,435` (a base, outlined) or `chr17:41,240,000-41,250,000` |
| Zoom | Mouse wheel (around the pointer; over rows that don't all fit, the wheel scrolls them and Ctrl+wheel zooms), `+`/`−`, double-click to zoom in (on a gene: fit the gene) |
| Step through variants | **< Variant** / **Variant >**, Ctrl+B / Ctrl+F |
| Pan | Drag (horizontal: genome; vertical: read rows), arrow keys, trackpad sideways scroll |
| Scroll rows | Drag up/down in a track, or `↑`/`↓`, PageUp/PageDown over it |
| Arrange tracks | Drag a track's header bar to move it (a blue line shows where it lands); drag the edge between two tracks to give one's height to the other, or into the room below the last track to move the track under it (the tracks further away stay put; a height you set is kept when zooming out hides the reads; the last track's bottom edge stops at the window's bottom), double-click an edge for the automatic height. Heights saved on a bigger screen shrink to fit. Saved with the session |
| Inspect | Hover a read (name, span, strand, MAPQ, flags, mismatches, the base under the pointer) or a coverage column (depth and allele counts) |
| Exact counts | Tick **Evidence**, click a base to outline it (again, or Esc, to clear) |
| Region stats | Shift-drag across the tracks |
| Two places at once | **Split**, then click a pane to steer it; **Link panes** to move them together |
| Work through a list | `--queue sites.bed`, then N / P and 1 / 2 / 3 |
| Measure | Ctrl+drag |

How it stays fast:

- **One loader thread per track**, so tracks load in parallel; the slowest track sets the
  wait, not the sum.
- **One GPU draw per track.** Reads, coverage, mismatches and indels are uploaded once per
  load as 32-byte instances whose x is in bases. Panning and zooming change only a uniform, so
  frame cost doesn't grow with the number of reads. Level of detail (sequence bars vs letters,
  hiding 1 bp insertions zoomed out) is decided in the shader, per instance.
- **Loading never blocks drawing.** A background thread loads about 3 view widths (at least
  20 kb), decodes reads once, derives coverage from them, and builds the instances. The UI
  uploads them when they arrive. While you keep moving, any result newer than what's on screen
  is shown, so data streams in mid-fling.
- **Only visible text is laid out per frame**: the ruler, and letters and labels for the rows
  and bases in view.
- **Wide loads decode on several threads**: spans of 500 kb or more are cut into chunks decoded
  in parallel and joined in order (the same reads, counts and coverage as one pass).
- **Zoomed out, reads carry only what can show**: at four or more bases a pixel a read's blocks
  join into bars across gaps under a pixel, and mismatches are drawn only in columns with a
  single-base variant in 35% of reads, not every sequencing error. A few megabases of nanopore
  reads went from ten million instances to two; zooming in rebuilds everything.
- **Nothing scans every loaded base per frame**: the deepest coverage in view comes from a
  per-kilobase index, splice junctions and VCF REF checks are counted once per load.
- Wider than 5 Mb: depth estimated from the index (above), in milliseconds.
- **Reads up to 1 Mb** by default (IGV shows 30 kb), set on the toolbar ("Reads up to"), with
  `--read-window`, or by an agent (MCP `display`), from 10 kb to 50 Mb and remembered. Cost
  grows with the reads in view; on 35x short reads (HG003 NovaSeq, about 265,000 reads per Mb,
  one track with the reference):

  | Window | Reads | Load | Frame median / 95th pct while panning | Worst pause | Peak memory |
  |---|---|---|---|---|---|
  | 100 kb | 26k | 0.12 s | 0.6 / 1.0 ms | 5 ms | 0.44 GB |
  | 500 kb | 133k | 0.34 s | 1.4 / 6.4 ms | 19 ms | 0.62 GB |
  | 1 Mb | 265k | 0.60 s | 2.5 / 6.3 ms | 37 ms | 0.97 GB |
  | 2 Mb | 529k | 1.0 s | 5.0 / 10.6 ms | 97 ms | 1.51 GB |
  | 5 Mb | 1.34M | 2.8 s | 11.9 / 29.5 ms | 230 ms | 3.45 GB |

  Sparse data (long reads, exomes, low-pass) goes much wider at the same cost. Tracks with more
  drawing instances than one GPU buffer holds are split across several, so no size is refused.
- Wider than the read window: binned depth only (up to 5 Mb), counted straight from each record's
  encoded CIGAR with no decoding or allocation (a 3 Mb exome span: ~0.25 s → ~0.15 s end to
  end, identical bins; `coverage_parity.rs` checks them against per-base coverage).
- The reference row's instances are built on the loader thread too; building them on the UI
  thread on every load was the main source of the occasional 4–5 ms frames (worst frames now
  2.4–3.3 ms).

Measured with `--bench` (200 frames each of panning, zooming, and fast flinging with background
loads; vsync off) on NA12878 exome BRCA1 data (13–15k reads, depth up to 709, 60–80k
instances):

| View | Frame CPU p50 / p95 | Wall p50 (fps) | Loads during bench |
|---|---|---|---|
| 101 bp, letters | 0.49 / 0.75 ms | 0.75 ms (~1,300) | 0 |
| 10 kb | 0.33–0.53 / 0.48–1.14 ms | 0.6–0.96 ms (~1,000–1,700) | 33–42, no frame spike |
| 1 Mb, binned depth | 0.47–0.59 / 0.61–1.04 ms | 0.73–0.93 ms (~1,100–1,350) | 4 |

Three GRCh38 exome tracks over BRCA1 (40k reads, 144k instances), 10 kb view: frame CPU p50
0.72–1.04 ms, p95 1.0–2.1 ms, worst 4.1 ms (650–930 fps), with 125 background loads during the
run. One of those tracks alone: p50 0.5–0.7 ms.

Loads: a 20 kb read span in about 19 ms, a 3 Mb binned span in about 240 ms. Launch to first
frame with data takes about 0.8 s.

Automation flags, for CI or checking a change: `--screenshot out.png` saves the real GPU frame
once data has loaded (in every pane), then exits; `--bench` prints the table above, then exits;
`--zoom-cycle` zooms all the way out and back in and reports each track's load state.
With split panes, scripted runs print where each pane ended up. They never overwrite the saved
last session.

## Splicing (RNA-seq)

`ggv query junctions` counts the splice junctions (introns skipped by `N` in the CIGAR) of one or
more files over a region, with the same read filters as everything else, and lines them up:

```sh
ggv query junctions -f GRCh38.fa -g gencode.gtf.gz affected.bam unaffected.bam chr17:7,668,421-7,687,490
```

Each junction gets its splice-site motif (GT/AG, GC/AG, AT/AC or non-canonical) and the strand it
implies (falling back to the aligner's `XS` tag, then the annotation), whether the annotation
knows it, and per file its reads and its **share** of every junction leaving its donor and
entering its acceptor: a shift in those shares between samples is aberrant splicing. On two
DRAGEN RNA BAMs, TP53 (both files) takes 0.3 s. `--format json` for scripts; MCP `junctions` gives
the same, defaulting to the window's files and view when attached.

`ggv render --sashimi` (MCP `render` with `sashimi`) adds a sashimi band under each file's
coverage: an arc per intron with its read count, + strand above the line in red, − below in blue
(`--min-junction-reads`, default 2). In the app, **Junctions** on the toolbar (or `display` with
`junctions`) draws the arcs over each track's coverage, counted from the reads that pass the
track's filters, so they agree with `ggv query junctions`.

## Signal tracks

bigWig (`.bw`), bedGraph (`.bedGraph`, `.bg`, gzipped too) and SEG (`.seg`, copy number) open
like any other file (command line, drag and drop, or MCP `open`). Signals draw as a bar chart
scaled to the values in view (negative values below the zero line); SEG files draw one heatmap
row per sample, blue for loss and red for gain (log2 ratios, full colour at ±1.5). Hover for the
value. bigWig files are read through their index, from the coarsest zoom level at most an eighth
of a pixel column for wide views, so a whole chromosome draws in milliseconds; bedGraph and SEG
are read whole on a background thread.

```sh
ggv-desktop sample.bam chip.bw atac.bedGraph.gz tumor.seg --genome hg38
ggv query signal chip.bw chr7:55,000,000-56,000,000 --bins 100 --format json
```

MCP has a `signal` tool (mean, min and max per bin, per file and SEG sample; attached, the
window's signal files and view). On UCSC's yeast GC bigWig, every raw value on chrI and chrM
matches UCSC's own API, and whole-chromosome bins from zoom levels are within 3% of the range
of means computed from raw values.
