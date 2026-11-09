"""Builds tests/data/hg002-revio from slices of PacBio's public Revio HG002 data.

Step 1 (needs samtools and bcftools; they read only the slices over HTTPS):

    V=https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38
    samtools view -b -o u.bam "$V.bam" chr20:58800000-58950000 chr20:31150000-31200000
    samtools sort -o HG002.revio.slices.bam u.bam && samtools index HG002.revio.slices.bam
    for k in deepvariant.phased pbsv; do
      bcftools view -r chr20:58800000-58950000,chr20:31150000-31200000 -Oz -o HG002.revio.$k.vcf.gz "$V.$k.vcf.gz"
      bcftools index -t HG002.revio.$k.vcf.gz
    done
    curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz && gunzip chr20.fa.gz
    curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

Step 2 (this script, Python 3 only):

    python scripts/make_revio_fixture.py SLICES_DIR tests/data/hg002-revio

copies the reads and calls, writes chr20 at full length with only the windows the reads reach
kept (N elsewhere) as BGZF with .fai and .gzi, and keeps the RefSeq transcripts that overlap
the windows.
"""

import gzip, shutil, struct, sys, zlib
from pathlib import Path

# The slices, widened by the longest read (about 32 kb) so every aligned base has its reference.
KEEP = [(31_110_000, 31_240_000), (58_760_000, 58_990_000)]
LINE = 60


def bgzf_block(data: bytes) -> bytes:
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    cdata = c.compress(data) + c.flush()
    bsize = 18 + len(cdata) + 8
    header = b'\x1f\x8b\x08\x04' + b'\0\0\0\0' + b'\0\xff' + struct.pack('<H', 6) + b'BC' + struct.pack('<HH', 2, bsize - 1)
    return header + cdata + struct.pack('<II', zlib.crc32(data) & 0xFFFFFFFF, len(data))


def write_bgzf(path: Path, data: bytes) -> list[tuple[int, int]]:
    """Writes `data` as BGZF; returns (compressed, uncompressed) offsets of every block start
    after the first, as the .gzi lists them."""
    offsets, coff = [], 0
    with open(path, 'wb') as f:
        for uoff in range(0, len(data), 0xFF00):
            if uoff:
                offsets.append((coff, uoff))
            block = bgzf_block(data[uoff:uoff + 0xFF00])
            f.write(block)
            coff += len(block)
        f.write(bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000'))  # EOF block
    return offsets


def masked_reference(src: Path, out: Path):
    name, seq = None, []
    with open(src) as f:
        for line in f:
            if line.startswith('>'):
                assert name is None, 'one sequence expected'
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    s = ''.join(seq)
    masked = bytearray(b'N' * len(s))
    for a, b in KEEP:
        masked[a:b] = s[a:b].encode()
    header = f'>{name}\n'.encode()
    body = b'\n'.join(bytes(masked[i:i + LINE]) for i in range(0, len(masked), LINE)) + b'\n'
    offsets = write_bgzf(out, header + body)
    Path(str(out) + '.fai').write_text(f'{name}\t{len(s)}\t{len(header)}\t{LINE}\t{LINE + 1}\n', newline='\n')
    with open(str(out) + '.gzi', 'wb') as f:
        f.write(struct.pack('<Q', len(offsets)))
        for c, u in offsets:
            f.write(struct.pack('<QQ', c, u))
    print(f'{out.name}: {len(s):,} bp, {sum(b - a for a, b in KEEP):,} kept, {out.stat().st_size:,} bytes')


def gene_subset(src: Path, out: Path):
    # Keep whole transcripts that overlap a window (every line of an overlapping transcript).
    keep_tx = set()
    with gzip.open(src, 'rt') as f:
        for line in f:
            p = line.split('\t')
            if len(p) > 8 and p[0] == 'chr20' and p[2] == 'transcript' and any(int(p[3]) <= b and int(p[4]) >= a for a, b in KEEP):
                keep_tx.add(p[8].split('transcript_id "')[1].split('"')[0])
    n = 0
    with gzip.open(src, 'rt') as f, open(out, 'w', newline='\n') as w:
        for line in f:
            p = line.split('\t')
            if len(p) > 8 and p[0] == 'chr20' and 'transcript_id "' in p[8] and p[8].split('transcript_id "')[1].split('"')[0] in keep_tx:
                w.write(line)
                n += 1
    print(f'{out.name}: {len(keep_tx)} transcripts, {n} lines')


def main():
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    for name in ['HG002.revio.slices.bam', 'HG002.revio.slices.bam.bai', 'HG002.revio.deepvariant.phased.vcf.gz',
                 'HG002.revio.deepvariant.phased.vcf.gz.tbi', 'HG002.revio.pbsv.vcf.gz', 'HG002.revio.pbsv.vcf.gz.tbi']:
        shutil.copyfile(src / name, out / name)
    masked_reference(src / 'GRCh38.chr20.fa', out / 'GRCh38.chr20.slices.fa.gz')
    gene_subset(src / 'hg38.ncbiRefSeq.gtf.gz', out / 'ncbiRefSeq.slices.gtf')


if __name__ == '__main__':
    main()
