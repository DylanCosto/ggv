# Agents (MCP)

`ggv mcp` serves the same data to AI agents over the Model Context Protocol (stdio), headless.
Agents don't get pictures; they get the numbers GGV draws from, computed by the same engines as
the desktop app and `ggv query`, so an agent and a person looking at the same base see the same
counts. Files given on the command line are the defaults for every call:

```sh
claude mcp add ggv -- ggv mcp -b tumor.bam -b normal.bam -f GRCh38.fa -v calls.vcf.gz -g genes.gtf.gz
```

or, for clients configured with JSON:

```json
{"mcpServers": {"ggv": {"command": "ggv", "args": ["mcp", "-b", "/data/tumor.bam", "-b", "/data/normal.bam",
  "-f", "/ref/GRCh38.fa", "-v", "/data/calls.vcf.gz", "-g", "/ref/genes.gtf.gz"]}}}
```

| Tool | Answers |
|---|---|
| `info`, `open` | the default files, their samples, read groups and main contigs; add or replace defaults |
| `evidence` | one base (or up to 50): per file, reads per allele by strand, and per allele the median MAPQ, median distance from the read end and soft-clipped share; quality means; **warnings** for one-strand support, alt near read ends, alt reads mapping worse, low depth, MAPQ 0, soft-clip pile-ups; the VCF records there; each change's HGVS consequence |
| `reads` | reads at a base or in a region (≤ 100 kb) grouped by allele, strand, HP tag, allele phase or not at all, with group statistics; raw reads only when `list` is given (≤ 50 per group) |
| `summarize` | `ggv query summary`: reads used and why others were dropped, depth statistics, breadth, the most variable positions |
| `coverage` | binned depth, counted up to 5 Mb and estimated from the index beyond (a whole chromosome in about 25 ms) |
| `variants`, `genes` | VCF records with genotypes; transcripts by region or by gene name/ID |
| `consequence` | HGVS `c.`/`p.` and effect of any REF→ALT (substitutions, and indels as `del`/`dup`/`ins`/`delins`/`fs` at the protein level too), with a warning if REF doesn't match the reference |
| `phase` | haplotypes from the reads' own alleles: sites, blocks, reads per haplotype, which allele each carries |
| `svs` | structural variants (≥ 50 bp) in a region (≤ 2 Mb) against the reads: each VCF record (sequence-resolved or symbolic) with, per file, supporting reads (a CIGAR deletion/insertion of about its size near the breakpoint, or a split alignment across it), reference reads, reads with another allele there, clipped reads, allele fraction, a genotype from the reads, and all of it by HP tag; plus clusters of large events in the reads that no record explains. On HG002 ONT against the T2T Q100 truth, the reads' genotype matched 34 of 54 SVs in ordinary regions, most of the rest being one SV listed once per haplotype (flagged as `same_site`) |
| `repeats` | tandem repeat lengths in the reads (e.g. HTT CAG, FMR1 CGG, VNTRs): each read spanning the repeat with flank gives its length; alleles by haplotype when reads are phased, else split at the largest gap, with range, reads, HP split, change from the reference and motif copies. HG002 ONT agrees with TRGT-style calls at 35 of 35 chr17 loci tested and at HTT (87 and 108 bp) |
| `junctions` | RNA-seq splice junctions in a region (≤ 10 Mb) or gene, every file side by side: intron, length, strand and motif (from the reference), known or novel (novel donor, novel acceptor, novel combination of known sites) with its genes, and per file reads, unique reads and the junction's share of its donor's and acceptor's reads |
| `compare` | positions where an allele's fraction differs between two files (tumor/normal) by more than the threshold |
| `render` | a PNG, for clients that can see images or to show a person |
| `signal` | binned values of bigWig / bedGraph / SEG tracks (mean, min, max; per SEG sample) |

Every tool takes `filter` (the desktop app's expression language: `"mapq >= 20 && !dup"`),
`min_mapq`, `include_duplicates`, `exclude_flags` and `min_base_quality`, applied to every count (bases below
the minimum quality are not observed: not in depth, alleles, phasing or pictures; `reads` lists
them in their own group), and every answer carries a `query` block with the full paths, filters
and region that produced it; each result names its file by path as well as by short name. Answers
are compact JSON, summaries first: `evidence` for two BAMs is about 500 tokens. `evidence` and
`consequence` give one transcript per gene (MANE Select, else MANE Plus Clinical, else canonical,
else the curated one with the longest CDS; MANE from the file's tags or GGV's built-in MANE v1.5
list; where the gene's main transcript doesn't reach the change, its best one there, curated before
predicted) and say how many others there were; `transcripts: "all"` lists every one, and `genes: []`
leaves annotations out. Open files are
kept between calls and annotations parsed once; on the HG002 fixtures every call takes 0–26 ms.
Errors come back as tool results the model can read (`contig chr99 is not in this file (it has
e.g. chr20)`), not protocol failures.

## Attached to the app

`ggv mcp --attach` makes the tools about what the person is looking at. Each open ggv-desktop
window listens on a loopback port (127.0.0.1 only, behind a random token kept in the user's own
config folder), and before every call `ggv mcp` reads the window's state, so "is this variant
real?" needs no coordinates or file names:

```sh
claude mcp add ggv -- ggv mcp --attach
```

Whatever a call leaves out comes from the active pane, and every answer's `desktop` block says
what was used and where it came from:

| Left out | Taken from the window |
|---|---|
| files | the pane's BAMs, reference, VCFs and annotations |
| the base (`evidence`, `consequence`) | the outlined base while it's on screen, else the view centre (what the Evidence panel shows); `evidence` takes every base of a selection of up to 50 bp made after it |
| `reads` | whichever the person set last: the selection, or the outlined base |
| the region (`summarize`, `coverage`, `variants`, `genes`, `phase`, `compare`) | the selection, else the view |
| `render` | the view: a picture of what the person sees |
| read filters | that file's filter in the app (flags, MAPQ, expression), with a note if the files' filters differ |

Three more tools appear: `view` (each pane's view, outlined base, selection, files, filters and
display options, and the other GGV windows), `navigate` (go to a locus, region or gene, or
select a region, which opens its statistics in the Region window) and `display` (change the
window's colouring, grouping, sort, soft clips, squished rows, transcripts and read filters, for
every file or one file's own expression, as the toolbar would; every value is checked before
any is applied, and the reply's `before` restores the previous settings). `open` opens files in the
window. Five more manage the window: `select` (select a region, or `clear`), `bookmarks` (list,
add with a note, remove, clear, or export as BED), `session` (save to or load from a `.ggv.json`
file), `tracks` (list the open tracks, close one, or move it up or down) and `review` (the
real/artifact/unsure decisions made in variant triage, per VCF).

Six more do in the window what the person would, through the same code as its buttons and
keys, each returning what changed or an error that says why nothing did:

| Tool | What it does |
|---|---|
| `review` | Variant triage: `start` a VCF, `decide` the current site real / artifact / unsure with a `note`, `note`, `next`, `prev`, `stop`; without an action, every VCF's decisions and the current site |
| `queue` | A review queue: `open` a list file or `paste` loci, `decide` pass / fail / unsure with a note, `next`, `prev`, `next_undecided`, `go` to an item, `export` as BED or TSV, `stop` |
| `panes` | Split view: `split` at a locus, `close`, `activate` a pane, `link` panes to scroll together |
| `action` | Any of the window's own actions by name (`back`, `forward`, `next_variant`, `next_diff`, `toggle_evidence`, `clear_selection`, `reload_all`, `copy_fasta`, `color:strand`…); without a name, the list. Actions that need a file dialog are refused in favour of the tools that take a path |
| `rows` | Scroll an alignment track's read rows (as the person does with the wheel or scrollbar over rows that don't fit, by dragging, or with ↑/↓, PageUp/PageDown): to a `row`, `by` some rows, `to` top or bottom, to the row holding a named `read` (outlined), or to a `group`; without one, each track's row count, the rows on screen and its groups |
| `search` | The window's sequence search: `start` (sequence, mismatches, scope), progress and hits, `go` to a hit, `stop` |
| `show_compare` | The visible tumor/normal comparison: tracks, threshold and minimum depth, `next` / `prev` difference, `close` |

`view` reports `api` (`version` and the commands this window understands) so a client can check
before relying on one, and each pane's `assembly` (build, how it is known, disagreements). A
review an agent finishes on its own, then hands back:

```text
review {action: "start", vcf: "calls.vcf.gz"}        → current.site: chr20:31157163 AGAATG>A
review {action: "decide", decision: "artifact", note: "strand bias"}   (saved; moves on)
review {}                                            → current.site: chr20:31157392 …
review {action: "decide", decision: "real"}
session {action: "save", path: "case.ggv.json"}       (files, panes, review file beside the VCF)
review {action: "stop"}                              (the window is the person's again)
```

Moving to the next site takes a moment: while `current.site` is null and `finding_next` is true,
ask again.

Two more make the agent part of the session rather than a bystander. `overlays` marks regions on
the view for the person: a tinted band across every track with a short label in the ruler
(`add` with a region, label and one of seven colours; `list`, `remove`, `clear`), shared by the
window's panes, saved with the session, and cleared from the toolbar's **Clear overlays**.
`events` reports what the person does as they do it — outlining a base, selecting or clearing a
region (once the drag ends), adding a bookmark, deciding a triage site — waiting up to 300 s for
the next (`wait_s`), from a `cursor` (`since`) and for some `kinds` only. An agent that waits for
`outline` and then calls `evidence` answers "is this real?" for whatever base the person
clicks. Only the person's actions are events: the agent's own commands, and the window loading,
are not. With several windows the focused one is used, else the last used; `--attach PID` or a
call's `window` names one, and `pane` picks a pane. Files given with `-b`/`-f`/... are used
when no window is open. Reading the window takes well under a millisecond, and arguments a call
does give always win.
