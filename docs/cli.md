# The command line

```sh
ggv info sample.bam
ggv query pileup  sample.bam chr7:55,191,822 -f ref.fa --format json
ggv query pileup  sample.bam chr7:55191000-55192000 --min-alt-fraction 0.2
ggv query summary sample.bam chr7:55191000-55192000 -f ref.fa --format json
```

```sh
ggv render sample.bam chr7:55,191,822 -f ref.fa -o snv.png            # ±50 bp, base outlined
ggv render sample.bam chr7:55191000-55241000 -f ref.fa -o view.svg    # SVG, real text
ggv render sample.bam chr7:55,191,822 -f ref.fa --scale 2 --color-by strand -o snv@2x.png
ggv render hifi.bam chr20:58888000-58890000 -f ref.fa --color-by base-mods -o methylation.png
```

`--color-by` takes `strand`, `haplotype` (HP tag), `insert-size` (pair orientation and insert
size) or `base-mods` (5mC and other MM/ML calls, red likely modified, blue not), with the app's
colours. `--group` puts reads in labelled bands as the app does: `strand`, `haplotype`,
`read-group`, `alleles` (phased from the reads' own alleles) or, for a single base, `base`;
`--sort` orders them by what they show at a single base. `--color-by base-mods --group
haplotype` shows allele-specific methylation (imprinting) at a glance. `ggv batch` takes
`--color-by` and `--group` too, and MCP `render` takes `color_by` and `group_by` (attached, it
uses the window's colouring and grouping).

Renders follow IGV's conventions: coverage with allele-coloured bars where ≥20% of reads
disagree with the reference, reads packed into rows (squished past 40 rows), mismatches
coloured A/C/G/T and shaded by base quality, letters when zoomed in, black lines for deletions
(labelled from 2 bp), purple bars for insertions (≥5 bp when zoomed out), and strand arrows.
Regions wider than `--read-window` (1 Mb, as in the app) show coverage only.

**Batch snapshots**: a picture and an evidence report for every place in a list, with no
window, on all cores:

```sh
ggv batch sites.bed -b tumor.bam -b normal.bam -f ref.fa -v calls.vcf.gz -g genes.gtf.gz -o review/
```

The list is BED, VCF, or one locus per line (`chr1:100`, `chr1:100-200`, text after it is a
label), gzipped or not. For each place GGV writes `NNN-contig-pos.png` (every BAM stacked, the
place outlined; single bases get ±50 bp) and `NNN-contig-pos.html` (the image embedded, the VCF
records there, the consequences of each change, allele counts by strand at a single base or a
region summary with its variable sites). `index.html` is a contact sheet of every place linking
to its report, and `summary.tsv` has one row per place (consequence, depth and top allele per
BAM). A place that fails (an unknown contig) is reported and skipped. Nine places over two
BAMs take about 0.3 s. The exit status is 0 when every place was written, 3 when some failed
(`summary.tsv` says which and why) and 1 when none was.

Each position reports depth, A/C/G/T/N by strand, deletions, insertions after the base,
reference skips, mean MAPQ and MAPQ-0 count, mean base quality, soft clips starting/ending
there, and the top alt allele with its fraction. `summary` adds read-filter counts, mean/median
depth, breadth at 1/10/20/30×, and the most variable sites.

Defaults count what IGV shows: unmapped, secondary, QC-fail and duplicate reads are dropped
(`--exclude-flags 0x704`), no MAPQ or base-quality floor, overlapping mates both count.
