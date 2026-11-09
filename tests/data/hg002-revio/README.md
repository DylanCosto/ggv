# HG002 PacBio Revio long-read fixture (GRCh38)

Two slices of real PacBio Revio HiFi data for GIAB HG002, small enough to keep in the
repository (about 4 MB), with everything long reads add over short ones: haplotype tags,
methylation calls, split reads and a structural variant on one haplotype. Built by
`scripts/make_revio_fixture.py` (its docstring has every command).

| Slice | What it tests |
|---|---|
| chr20:58,800,000-58,950,000 (GNAS) | 5mC calls (MM/ML) on CpGs; imprinting: at GNAS A/B (58,888,000-58,890,000) one haplotype is about 17% methylated and the other 88%; 12 phased heterozygous SNPs; a 1 kb insertion (pbsv) near 58,830,532 |
| chr20:31,150,000-31,200,000 | heterozygous 6,802 bp deletion pbsv.DEL.567 (31,170,827-31,177,629): every read carrying it is HP 1, every intact read spanning it HP 2; split reads (SA, supplementary records) |

Files:

- `HG002.revio.slices.bam` (+ `.bai`): 535 records (419 primary, 116 supplementary), median read
  length 15.7 kb, all tags as published (HP/PS on 253 records, MM/ML on all, SA on 159).
- `HG002.revio.deepvariant.phased.vcf.gz`, `HG002.revio.pbsv.vcf.gz` (+ `.tbi`): the same run's
  DeepVariant phased small variants and pbsv structural variants, cut to the slices.
- `GRCh38.chr20.slices.fa.gz` (+ `.fai`, `.gzi`): UCSC hg38 chr20 at full length, with only
  31,110,000-31,240,000 and 58,760,000-58,990,000 kept (the slices plus the longest read) and N
  elsewhere, BGZF-compressed.
- `ncbiRefSeq.slices.gtf`: the UCSC hg38 RefSeq transcripts overlapping the kept windows.

Sources (public):

- Reads and calls: PacBio Revio 2022Q4 HG002 rep1,
  `https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.*`
  (pbmm2 alignment, DeepVariant, WhatsHap haplotagging, pbsv).
- Reference: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz`.
- Genes: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz`.

Used by `crates/ggv-layout/tests/long_reads.rs`. To look at it:

```sh
ggv-desktop tests/data/hg002-revio/HG002.revio.slices.bam -f tests/data/hg002-revio/GRCh38.chr20.slices.fa.gz \
  -v tests/data/hg002-revio/HG002.revio.pbsv.vcf.gz -g tests/data/hg002-revio/ncbiRefSeq.slices.gtf \
  chr20:31,165,000-31,184,000 --group haplotype --color base-mods
```
