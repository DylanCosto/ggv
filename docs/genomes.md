# Genomes and built-in tracks

Instead of finding a reference and a gene annotation, name a built-in genome; GGV downloads
UCSC's `.2bit` and NCBI RefSeq genes once and keeps them:

```sh
ggv-desktop sample.bam --genome hg38          # or Reference… → a genome in the app
ggv query pileup sample.bam chr7:55,191,822 --genome hg38
ggv genomes                                    # what's built in, and what's downloaded
ggv genomes --fetch T2T                        # download ahead of time
```

| Genome | Also | Download |
|---|---|---|
| `hg38` | GRCh38, human | 877 MB |
| `hg19` | GRCh37, b37 | 836 MB |
| `hs1` | T2T, CHM13 | 852 MB |
| `mm39` | GRCm39, mouse | 745 MB |
| `mm10` | GRCm38 | 743 MB |
| `sacCer3` | yeast, R64 | 3 MB |

`--genome` works for the app, `ggv query`, `ggv render`, `ggv batch` and `ggv mcp` (reference and,
where used, genes); `-f`/`-g` still win. Files are kept in `%LOCALAPPDATA%\ggv\genomes`,
`~/Library/Caches/ggv/genomes` or `$XDG_CACHE_HOME/ggv/genomes` (`GGV_GENOME_DIR` overrides);
downloads use the system's `curl` and resume if interrupted. `.2bit` references can be read at
any position without unpacking (UCSC's FASTA downloads are plain gzip), and any `.2bit` works as
`-f` for BAM; CRAM decoding still needs a FASTA. On the long-read fixture, `--genome hg38` gives
byte-identical pileups, summaries and pictures to the chr20 FASTA, at 30 ms a query. The app
shows the reference row (and codons in genes) from the genome alone, without a BAM.

## Built-in tracks

The app's **Tracks…** menu lists annotation, variant and pathogenicity tracks for the open
genome's build (found from chromosome 1's length, whatever reference is open), and `--track`
names them; each downloads once into `<genome cache>/<build>/tracks`:

```sh
ggv-desktop sample.bam -f GRCh38.fa --track mane --track clinvar --track alphamissense
```

| Track | Builds | What | Download |
|---|---|---|---|
| `gencode` | hg38, hg19 (lifted), mm39, mm10 | GENCODE comprehensive | 30–160 MB |
| `gencode-basic` | hg38 | GENCODE basic | 80 MB |
| `mane` | hg38 | MANE Select / Plus Clinical | 9 MB |
| `refseq` | all | NCBI RefSeq, as UCSC packages it | 0.4–42 MB |
| `clinvar` | hg38, hg19 | ClinVar VCF, newest dated weekly release | 200 MB |
| `alphamissense` | hg38, hg19 | AlphaMissense, one bigWig per alternate base | 1.8 GB |
| `revel` | hg38, hg19 | REVEL, one bigWig per alternate base | 1.9 GB |

The scores open as four signal tracks (the score of a change to A, C, G or T at each base).
ClinVar's sites are coloured by clinical significance (pathogenic red, uncertain dark blue,
benign green, conflicting light blue, other grey, with a legend in the track header); any VCF
can be coloured by an INFO field instead of allele frequency from the track's right-click menu
or MCP `display` (`variant_color`). Clicking a variant opens its details: ClinVar's summary
(significance, review stars, condition, gene, consequence), the consequences in the open gene
tracks, genotypes, and every INFO field with the header's description.
CADD (39 GB) and phyloP (10 GB) are left out until GGV reads remote files.
