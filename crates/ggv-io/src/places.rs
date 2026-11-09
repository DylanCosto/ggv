//! Lists of places: BED, VCF, `chr1:100-200` lines or pasted text, for review queues and
//! batch snapshots.

use std::{io::Read, path::Path};

use ggv_model::{Contig, Locus, Region};

/// A place from a list. `end: None` means to the end of the contig (a bare contig name).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Place {
    pub contig: String,
    /// 0-based.
    pub start: u64,
    /// 0-based exclusive.
    pub end: Option<u64>,
    pub label: String,
    /// `REF>ALT` for a VCF row, so variants at one site stay apart.
    pub alleles: Option<String>,
    /// The locus as written, for a place given as `chr1:100-200` or a bare name: a contig name
    /// may itself contain `:` (`HLA-A*01:01:01:01`), which only a file's contigs can settle.
    pub written: Option<String>,
    /// A zero-length BED interval: an insertion site, right after the base at `start`.
    pub insertion: bool,
    /// The other reading of a row that is both a locus with a label and a BED/VCF row whose
    /// contig has a `:` (`chr20:58909366<TAB>30<TAB>99`, `HLA-A*01:01:01:01<TAB>9<TAB>29`):
    /// [`Place::settle`] picks one with a file's contigs.
    pub other: Option<Box<Place>>,
}

impl Place {
    /// A stable key: `chr1:101-200` (1-based, inclusive).
    pub fn key(&self) -> String {
        match self.end {
            Some(e) if e == self.start + 1 => format!("{}:{}", self.contig, self.start + 1),
            Some(e) => format!("{}:{}-{}", self.contig, self.start + 1, e),
            None => self.contig.clone(),
        }
    }

    /// Which reading of an ambiguous row (see `other`) is meant in a file with these contigs:
    /// the BED/VCF one when its first column is exactly a contig there (or an alias of one),
    /// else the locus. Without a file, [`parse_lines`] put the likelier reading first.
    pub fn settle(&self, contigs: &[Contig]) -> &Place {
        let Some(other) = &self.other else { return self };
        let (row, locus) = if self.written.is_none() { (self, &**other) } else { (&**other, self) };
        let known = |n: &str| contigs.iter().any(|c| c.name == n) || ggv_model::contig_aliases(n).iter().any(|a| contigs.iter().any(|c| &c.name == a));
        if known(&row.contig) { row } else { locus }
    }

    /// The place in a file with these contigs (trying `chr` aliases), [settled](Self::settle)
    /// first. A locus written as exactly one of its contig names is that whole contig, as `GGV
    /// query` reads the same text.
    pub fn resolve(&self, contigs: &[Contig]) -> crate::Result<Region> {
        let p = self.settle(contigs);
        if !std::ptr::eq(p, self) {
            return p.resolve(contigs);
        }
        match &self.written {
            Some(text) => crate::resolve_locus_str(contigs, text),
            None => crate::resolve_locus(contigs, &Locus { contig: self.contig.clone(), start: Some(self.start + 1), end: self.end }),
        }
    }
}

/// One line of a list that is neither a place nor ignorable (blank, comment, header).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unreadable {
    /// 1-based.
    pub line: usize,
    pub text: String,
}

fn number(s: &str) -> Option<u64> {
    let d: String = s.chars().filter(|c| *c != ',' && *c != '_').collect();
    if d.is_empty() || !d.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    d.parse().ok()
}

/// REF and ALT columns as a VCF row has them: REF bases, ALT bases, symbolic alleles,
/// breakends, `*` or `.`, comma-separated. A BED row's fifth column is a score, never these.
fn vcf_alleles(reference: &str, alt: &str) -> bool {
    let bases = |s: &str| !s.is_empty() && s.bytes().all(|b| b"ACGTNacgtn".contains(&b));
    bases(reference)
        && alt.split(',').all(|a| {
            bases(a) || a == "*" || a == "." || (a.starts_with('<') && a.ends_with('>')) || a.contains('[') || a.contains(']')
        })
}

/// What kind of list a text is, when its file name (or the caller) says.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// Decide per line (pasted text, `.txt`, `.tsv`, ...).
    Auto,
    Bed,
    Vcf,
}

impl Format {
    /// From a file name: `.bed` and `.vcf`, also gzipped; anything else is `Auto`.
    pub fn of_path(path: &Path) -> Self {
        let name = path.file_name().map(|n| n.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
        let name = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".bgz")).unwrap_or(&name);
        if name.ends_with(".vcf") {
            Self::Vcf
        } else if name.ends_with(".bed") {
            Self::Bed
        } else {
            Self::Auto
        }
    }
}

/// Whether a headerless row with a numeric third column is VCF rather than BED: the eight
/// fixed VCF columns, REF/ALT alleles, and a FILTER where BED12 has a number (thickStart). A
/// BED row with allele-like name and score (`chr1 109 120 AC . +`) stays BED.
fn vcf_row(f: &[&str]) -> bool {
    f.len() >= 8 && vcf_alleles(f[3], f[4]) && number(f[6]).is_none()
}

/// A symbolic allele's INFO/END (1-based, inclusive), when valid for a record at `pos`.
fn symbolic_end(f: &[&str], pos: u64) -> Option<u64> {
    if !f[4].split(',').any(|a| a.starts_with('<')) {
        return None;
    }
    let end = f.get(7)?.split(';').find_map(|kv| kv.strip_prefix("END="))?;
    number(end).filter(|&e| e >= pos)
}

/// Reads places from text: one per line, as a locus (`chr1:100-200 note`), BED (`chr1 99 200
/// name`), VCF (`chr1 100 id REF ALT ...`) or `contig position`. Blank lines, `#` comments and
/// `track`/`browser` lines are ignored; other unreadable lines are counted.
///
/// A zero-length BED interval (`chr1 105 105`, an insertion site) is the base before it, where
/// VCF anchors an insertion and a pileup counts it (`INS after`).
pub fn parse(text: &str) -> (Vec<Place>, usize) {
    parse_as(text, Format::Auto)
}

/// [`parse`] for a list of known `format` ([`Format::of_path`]).
///
/// A VCF row is recognised by the format, the file's `##fileformat=VCF` or `#CHROM` header, or
/// else by its columns (a non-numeric ID, or a whole VCF row: see `vcf_row`), so a numeric ID
/// is never taken for a BED end coordinate nor a BED row for VCF. A symbolic allele's place
/// runs to its INFO/END; a missing or invalid END (before POS) leaves it at REF's length.
pub fn parse_as(text: &str, format: Format) -> (Vec<Place>, usize) {
    let mut items = Vec::new();
    let mut skipped = 0;
    for e in parse_lines(text, format) {
        match e {
            Ok(p) => items.push(p),
            Err(_) => skipped += 1,
        }
    }
    (items, skipped)
}

/// [`parse_as`], keeping each unreadable line (with its number) in list order, for callers that
/// report them.
pub fn parse_lines(text: &str, format: Format) -> Vec<Result<Place, Unreadable>> {
    let mut out = Vec::new();
    let mut vcf = format == Format::Vcf;
    // Notepad, Excel and Windows PowerShell start UTF-8 text with a byte-order mark, which
    // `trim` keeps: it would become part of the first contig name.
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.starts_with("##fileformat=VCF") || line.starts_with("#CHROM") {
            vcf = true;
        }
        if line.is_empty() || line.starts_with('#') || line.starts_with("track ") || line.starts_with("browser ") {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        // A whole BED or VCF row is one even when the contig has a `:` (`HLA-A*01:01:01:01 9
        // 29`). Anything less after a `chr1:105` is its label (`2 of 3 callers`, `30 25`), and
        // BED and VCF are tab-separated, so a space after the locus starts a label too.
        let bed = f.len() >= 3 && number(f[1]).zip(number(f[2])).is_some_and(|(s, e)| e >= s);
        let vcf_like = f.len() >= 5 && number(f[1]).is_some_and(|p| p > 0) && (vcf || vcf_alleles(f[3], f[4]));
        let tabular = (bed || vcf_like) && !line[f[0].len()..].starts_with(' ');
        let item = if f.len() == 1 || (f[0].contains(':') && !tabular) {
            locus_place(&f)
        } else if f[0].contains(':') && format == Format::Auto && !vcf {
            // `chr20:58909366<TAB>30<TAB>99` (a TSV of locus, depth, GQ) is also a BED row on a
            // contig named `chr20:58909366`. Only a file's contigs can say ([`Place::settle`]);
            // until then a locus whose own contig has a `:` (an HLA allele name, never a
            // `chr1`-style name) reads as that contig's row, any other as a labelled locus.
            match (locus_place(&f).filter(|l| l.end.is_some()), row_place(&f, vcf, format)) {
                (Some(l), Some(r)) => Some(if l.contig.contains(':') {
                    Place { other: Some(Box::new(l)), ..r }
                } else {
                    Place { other: Some(Box::new(r)), ..l }
                }),
                (l, r) => r.or(l),
            }
        } else {
            row_place(&f, vcf, format)
        };
        out.push(item.ok_or_else(|| Unreadable { line: n + 1, text: line.to_owned() }));
    }
    out
}

/// The words after the first `from`, as a label.
fn rest(f: &[&str], from: usize) -> String {
    f.get(from..).map(|r| r.join(" ")).unwrap_or_default()
}

/// A line read as a locus (`chr1:100-200`, a bare contig name) and its label.
fn locus_place(f: &[&str]) -> Option<Place> {
    f[0].parse::<Locus>().ok().map(|l| Place {
        contig: l.contig,
        start: l.start.map_or(0, |s| s - 1),
        end: match (l.start, l.end) {
            (_, Some(e)) => Some(e),
            (Some(s), None) => Some(s),
            (None, None) => None,
        },
        label: rest(f, 1),
        alleles: None,
        written: Some(f[0].to_owned()),
        insertion: false,
        other: None,
    })
}

/// A line read as a VCF, BED or `contig position` row.
fn row_place(f: &[&str], vcf: bool, format: Format) -> Option<Place> {
    if f.len() >= 5 && number(f[1]).is_some() && (vcf || (format == Format::Auto && vcf_alleles(f[3], f[4]) && (number(f[2]).is_none() || vcf_row(f)))) {
        // VCF: CHROM POS ID REF ALT.
        number(f[1]).filter(|&p| p > 0).map(|p| Place {
            contig: f[0].to_owned(),
            start: p - 1,
            end: Some(symbolic_end(f, p).unwrap_or(p - 1 + f[3].len().max(1) as u64)),
            label: if f[2] == "." { format!("{}>{}", f[3], f[4]) } else { f[2].to_owned() },
            alleles: Some(format!("{}>{}", f[3], f[4])),
            written: None,
            insertion: false,
            other: None,
        })
    } else if f.len() >= 3 && number(f[1]).is_some() && number(f[2]).is_some() {
        let (s, e) = (number(f[1]).unwrap(), number(f[2]).unwrap());
        // Zero length is an insertion site: the base before it (the first base at 0).
        let insertion = e == s;
        let (s, e) = if insertion { (s.saturating_sub(1), s.max(1)) } else { (s, e) };
        (e > s).then(|| Place { contig: f[0].to_owned(), start: s, end: Some(e), label: rest(f, 3), alleles: None, written: None, insertion, other: None })
    } else if f.len() >= 2 && number(f[1]).is_some_and(|p| p > 0) {
        let p = number(f[1]).unwrap();
        Some(Place { contig: f[0].to_owned(), start: p - 1, end: Some(p), label: rest(f, 2), alleles: None, written: None, insertion: false, other: None })
    } else {
        None
    }
}

/// Reads a list file, gunzipping `.gz`/`.bgz`; UTF-8, or UTF-16 with a byte-order mark (what
/// Windows PowerShell 5.1's `>` and Excel's "Unicode Text" write).
pub fn read_file(path: &Path) -> std::io::Result<String> {
    let name = path.to_string_lossy().to_ascii_lowercase();
    let file = std::fs::File::open(path)?;
    let mut bytes = Vec::new();
    if name.ends_with(".gz") || name.ends_with(".bgz") {
        flate2::read::MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else {
        std::io::BufReader::new(file).read_to_end(&mut bytes)?;
    }
    decode(bytes)
}

/// Text from bytes: UTF-16 LE/BE after its byte-order mark, else UTF-8 (a UTF-8 mark is left
/// for [`parse_lines`]).
fn decode(bytes: Vec<u8>) -> std::io::Result<String> {
    let utf16 = |rest: &[u8], unit: fn([u8; 2]) -> u16| {
        let units: Vec<u16> = rest.chunks(2).map(|c| unit([c[0], *c.get(1).unwrap_or(&0)])).collect();
        String::from_utf16(&units).map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid UTF-16 text"))
    };
    match bytes.as_slice() {
        [0xFF, 0xFE, rest @ ..] => utf16(rest, u16::from_le_bytes),
        [0xFE, 0xFF, rest @ ..] => utf16(rest, u16::from_be_bytes),
        _ => String::from_utf8(bytes).map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData, "stream did not contain valid UTF-8")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_every_list_format() {
        let text = "\
track name=x
# a comment
chr1\t99\t200\tpeak 1
chr2:1,000-2,000 odd coverage
chr3:500
chrX
chr4\t1234\trs1\tA\tG\t50\tPASS\t.
chr5\t1234\t.\tAC\tA\t.\t.\t.
chr6 777 het
not a place at all
";
        let (items, skipped) = parse(text);
        assert_eq!(skipped, 1);
        let keys: Vec<String> = items.iter().map(Place::key).collect();
        assert_eq!(keys, ["chr1:100-200", "chr2:1000-2000", "chr3:500", "chrX", "chr4:1234", "chr5:1234-1235", "chr6:777"]);
        assert_eq!(items[0].label, "peak 1");
        assert_eq!(items[1].label, "odd coverage");
        assert_eq!(items[4].label, "rs1");
        assert_eq!(items[5].label, "AC>A");
        assert_eq!(items[3].end, None);
        assert_eq!(items[4].alleles.as_deref(), Some("A>G"));
        assert_eq!(items[0].alleles, None);
    }

    #[test]
    fn numeric_vcf_ids_are_not_bed_ends() {
        // A VCF whose IDs are numbers is still a VCF, not a BED with an end column: with a
        // header, and without one.
        let with_header = "##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr20	30006342	30006343	G	A	50	PASS	.
chr20	30006342	12345	G	A	50	PASS	.
";
        let bare = "chr20	30006342	30006343	G	A	50	PASS	.
chr20	100	5	GT	<DEL>	.	.	.
";
        for text in [with_header, bare] {
            let (items, skipped) = parse(text);
            assert_eq!(skipped, 0);
            assert_eq!(items[0].key(), "chr20:30006342");
            assert_eq!(items[0].label, "30006343");
        }
        assert_eq!(parse(bare).0[1].key(), "chr20:100-101");
        // BED with score and strand columns stays BED.
        let (items, _) = parse("chr1	99	200	ACGT	500	+
");
        assert_eq!(items[0].key(), "chr1:100-200");
    }

    #[test]
    fn symbolic_vcf_end_is_the_place() {
        // A symbolic deletion and the same deletion as sequence give one place.
        let symbolic = "chrQA\t120\tdel\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=170\n";
        let sequence = format!("chrQA\t120\tdel\tA{}\tA\t.\tPASS\t.\n", "C".repeat(50));
        assert_eq!(parse(symbolic).0[0].key(), "chrQA:120-170");
        assert_eq!(parse(&sequence).0[0].key(), "chrQA:120-170");
        // A malformed or backwards END falls back to REF.
        assert_eq!(parse("chrQA\t120\tdel\tA\t<DEL>\t.\tPASS\tEND=x1\n").0[0].key(), "chrQA:120");
        assert_eq!(parse("chrQA\t120\tdel\tA\t<DEL>\t.\tPASS\tEND=100\n").0[0].key(), "chrQA:120");
        // END on a sequence allele is not used.
        assert_eq!(parse("chrQA\t120\tx\tA\tG\t.\tPASS\tEND=170\n").0[0].key(), "chrQA:120");
    }

    #[test]
    fn bed_rows_with_allele_like_fields_stay_bed() {
        let row = "chrQA\t109\t120\tAC\t.\t+\n";
        for (items, _) in [parse(row), parse_as(row, Format::Bed)] {
            assert_eq!(items[0].key(), "chrQA:110-120");
            assert_eq!(items[0].label, "AC . +");
        }
        // BED12 with allele-like name and a dot score.
        assert_eq!(parse("chr1\t99\t200\tAC\t.\t+\t99\t200\t0\t1\t101,\t0,\n").0[0].key(), "chr1:100-200");
        // A .vcf file is VCF even without a header.
        assert_eq!(parse_as("chr1\t100\t5\tGT\tG\n", Format::Vcf).0[0].key(), "chr1:100-101");
        assert_eq!(Format::of_path(Path::new("x/sites.BED.gz")), Format::Bed);
        assert_eq!(Format::of_path(Path::new("calls.vcf.bgz")), Format::Vcf);
        assert_eq!(Format::of_path(Path::new("list.txt")), Format::Auto);
    }

    #[test]
    fn a_byte_order_mark_is_not_part_of_the_first_place() {
        // Windows PowerShell `>` and Notepad write one; the first contig was `\u{feff}chr20`.
        let (items, skipped) = parse("\u{feff}chr20:58888000\r\nchr20:58889000\r\n");
        assert_eq!((items[0].key().as_str(), skipped), ("chr20:58888000", 0));
        let (items, _) = parse_as("\u{feff}chr20\t58887999\t58888000\tsite1\n", Format::Bed);
        assert_eq!(items[0].key(), "chr20:58888000");
        // A VCF header after the mark is still a header.
        let (items, skipped) = parse("\u{feff}##fileformat=VCFv4.2\nchr1\t100\t5\tG\tA\t.\t.\t.\n");
        assert_eq!((items[0].key().as_str(), skipped), ("chr1:100", 0));
    }

    #[test]
    fn contigs_with_colons() {
        let contigs = [Contig { name: "HLA-A*01:01:01:01".into(), length: 100 }, Contig { name: "HLA-A*01:01:01".into(), length: 50 }];
        // A BED row on an HLA contig is BED, not the locus `HLA-A*01:01:01:1` labelled "9 29".
        let (items, _) = parse("HLA-A*01:01:01:01\t9\t29\tbedrow\n");
        assert_eq!((items[0].contig.as_str(), items[0].start, items[0].end), ("HLA-A*01:01:01:01", 9, Some(29)));
        assert_eq!(items[0].label, "bedrow");
        // So is a VCF row.
        assert_eq!(parse("HLA-A*01:01:01:01\t15\trs1\tA\tG\t.\t.\t.\n").0[0].key(), "HLA-A*01:01:01:01:15");
        // A bare name that is a contig is that whole contig, even when a shorter one exists.
        let (items, _) = parse("HLA-A*01:01:01:01\nHLA-A*01:01:01:01:15 locusline\nchr1:100-200 5\n");
        let r = items[0].resolve(&contigs).unwrap();
        assert_eq!((r.contig.as_str(), r.start, r.end), ("HLA-A*01:01:01:01", 0, 100));
        let r = items[1].resolve(&contigs).unwrap();
        assert_eq!((r.contig.as_str(), r.start, r.end), ("HLA-A*01:01:01:01", 14, 15));
        assert_eq!(items[1].label, "locusline");
        // A locus with a numeric label stays a locus.
        assert_eq!((items[2].key().as_str(), items[2].label.as_str()), ("chr1:100-200", "5"));
    }

    #[test]
    fn a_locus_label_may_start_with_numbers() {
        // bughunt3: `2 of 3 ...` read as VCF (POS 2, ID `of`), `30 25` as a backwards BED row.
        let text = "chrQ:105 2 of 3 callers agree\nchrQ:105 2 callers\nchrQ:120 30 25\nchrQ:120 30 99\nchrQ:130\t2\tof 3\n";
        let (items, skipped) = parse(text);
        assert_eq!(skipped, 0);
        let got: Vec<(String, &str)> = items.iter().map(|p| (p.key(), p.label.as_str())).collect();
        assert_eq!(got, [("chrQ:105".into(), "2 of 3 callers agree"), ("chrQ:105".into(), "2 callers"), ("chrQ:120".into(), "30 25"), ("chrQ:120".into(), "30 99"), ("chrQ:130".into(), "2 of 3")]);
        // Without a `:`, REF/ALT that aren't alleles don't make a VCF row either.
        let (items, _) = parse("chr6 777 of 3 callers\n");
        assert_eq!((items[0].key().as_str(), items[0].label.as_str()), ("chr6:777", "of 3 callers"));
    }

    #[test]
    fn a_tsv_locus_with_numbers_is_a_locus_unless_it_is_a_contig() {
        // bughunt4: `30<TAB>99` after a locus made a BED row on contig `chr20:58909366`.
        let (items, skipped) = parse("chr20:58909366	30	99
chr20:58909366	30	25
");
        assert_eq!(skipped, 0);
        for p in &items {
            assert_eq!((p.key().as_str(), p.written.as_deref()), ("chr20:58909366", Some("chr20:58909366")));
        }
        assert_eq!(items[0].label, "30 99");
        // With a file's contigs: the locus, unless the first column is itself a contig there.
        let chr20 = [Contig { name: "chr20".into(), length: 64_444_167 }];
        let r = items[0].resolve(&chr20).unwrap();
        assert_eq!((r.contig.as_str(), r.start, r.end), ("chr20", 58_909_365, 58_909_366));
        assert_eq!(items[0].settle(&chr20).label, "30 99");
        let odd = [Contig { name: "chr20:58909366".into(), length: 1000 }];
        let p = items[0].settle(&odd);
        assert_eq!((p.contig.as_str(), p.start, p.end, p.label.as_str()), ("chr20:58909366", 30, Some(99), ""));
        assert_eq!(items[0].resolve(&odd).unwrap().start, 30);
        // An HLA-shaped first column reads as a BED row until contigs say otherwise.
        let (items, _) = parse("HLA-A*01:01:01:01	9	29	bedrow
");
        let hla = [Contig { name: "HLA-A*01:01:01".into(), length: 100 }];
        let p = items[0].settle(&hla);
        assert_eq!((p.contig.as_str(), p.start, p.label.as_str()), ("HLA-A*01:01:01", 0, "9 29 bedrow"));
        // A .bed file or a VCF is never a locus list.
        assert_eq!(parse_as("chr20:5	30	99
", Format::Bed).0[0].other, None);
    }

    #[test]
    fn utf16_lists_are_decoded() {
        // Windows PowerShell 5.1 `>` / `Out-File` writes UTF-16 LE with a byte-order mark.
        let text = "chrQ:15\r\nchrQ:100-110\r\n";
        let le: Vec<u8> = [0xFF, 0xFE].into_iter().chain(text.encode_utf16().flat_map(u16::to_le_bytes)).collect();
        let be: Vec<u8> = [0xFE, 0xFF].into_iter().chain(text.encode_utf16().flat_map(u16::to_be_bytes)).collect();
        let dir = std::env::temp_dir().join(format!("ggv-places-utf16-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for (name, bytes) in [("le.txt", le), ("be.txt", be), ("u8.txt", [&[0xEF, 0xBB, 0xBF][..], text.as_bytes()].concat())] {
            let p = dir.join(name);
            std::fs::write(&p, bytes).unwrap();
            let keys: Vec<String> = parse(&read_file(&p).unwrap()).0.iter().map(Place::key).collect();
            assert_eq!(keys, ["chrQ:15", "chrQ:100-110"], "{name}");
        }
        std::fs::write(dir.join("bad.txt"), [0x63, 0xFF, 0x0A]).unwrap();
        assert!(read_file(&dir.join("bad.txt")).unwrap_err().to_string().contains("UTF-8"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn zero_length_bed_is_an_insertion_site_and_bad_lines_are_kept() {
        let text = "chrQ:200\nchrQ\t105\t105\tinsertion-site\nchrQ\t0\t0\nchrQ 10x\nchrQ\t50\t40\n";
        let (items, skipped) = parse(text);
        let keys: Vec<String> = items.iter().map(Place::key).collect();
        assert_eq!(keys, ["chrQ:200", "chrQ:105", "chrQ:1"]);
        assert_eq!(items.iter().map(|p| p.insertion).collect::<Vec<_>>(), [false, true, true]);
        assert_eq!(skipped, 2);
        let lines = parse_lines(text, Format::Auto);
        let bad: Vec<(usize, &str)> = lines.iter().filter_map(|l| l.as_ref().err()).map(|u| (u.line, u.text.as_str())).collect();
        assert_eq!(bad, [(4, "chrQ 10x"), (5, "chrQ\t50\t40")]);
    }
}
