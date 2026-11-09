//! Deep links and reopen commands for a view.
//!
//! `ggv://open?file=/data/a.bam&file=/data/calls.vcf.gz&ref=/data/hg38.fa&locus=chr17:43,044,295-43,125,483`
//! opens those files at that place (once the scheme is registered, from a browser, chat or
//! notebook too). Each file is its own `file=` value, so paths may contain commas.

use std::path::PathBuf;

use crate::session::Session;

/// What a link asks for.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Link {
    pub files: Vec<PathBuf>,
    pub reference: Option<PathBuf>,
    pub locus: Option<String>,
    pub highlight: Option<String>,
}

fn encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b"-._~/:,".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

fn decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'%' if i + 2 < b.len() => {
                let hex = |c: u8| (c as char).to_digit(16);
                match (hex(b[i + 1]), hex(b[i + 2])) {
                    (Some(h), Some(l)) => {
                        out.push((h * 16 + l) as u8);
                        i += 3;
                        continue;
                    }
                    _ => out.push(b'%'),
                }
            }
            b'+' => out.push(b' '),
            c => out.push(c),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// A link to `s`'s files, view and outlined base.
pub fn to_link(s: &Session) -> String {
    let mut parts = Vec::new();
    // Every kind of track (signal tracks were left out).
    for f in s.bams.iter().chain(&s.vcfs).chain(&s.genes).chain(&s.signals) {
        parts.push(format!("file={}", encode(&f.to_string_lossy())));
    }
    if let Some(r) = &s.reference {
        parts.push(format!("ref={}", encode(&r.to_string_lossy())));
    }
    if let Some(l) = &s.locus {
        parts.push(format!("locus={}", encode(l)));
    }
    if let Some(h) = &s.highlight {
        parts.push(format!("highlight={}", encode(h)));
    }
    format!("ggv://open?{}", parts.join("&"))
}

/// The schemes links are read in: `ggv://`, and `gv://` from before GGV was renamed (links
/// copied then still open).
const SCHEMES: [&str; 2] = ["ggv://", "gv://"];

/// What follows the scheme, if `s` is a link.
fn after_scheme(s: &str) -> Option<&str> {
    SCHEMES.iter().find_map(|scheme| s.get(..scheme.len()).filter(|p| p.eq_ignore_ascii_case(scheme)).map(|_| &s[scheme.len()..]))
}

pub fn is_link(s: &str) -> bool {
    after_scheme(s).is_some()
}

pub fn parse(link: &str) -> Result<Link, String> {
    let Some(rest) = after_scheme(link) else {
        return Err("not an ggv:// link".into());
    };
    let (action, query) = rest.split_once('?').unwrap_or((rest, ""));
    if !action.trim_end_matches('/').eq_ignore_ascii_case("open") {
        return Err(format!("unknown ggv:// action '{action}' (expected open)"));
    }
    let mut out = Link::default();
    for pair in query.split('&').filter(|p| !p.is_empty()) {
        let (k, v) = pair.split_once('=').unwrap_or((pair, ""));
        let v = decode(v);
        match k {
            "file" => out.files.push(v.into()),
            "ref" | "reference" => out.reference = Some(v.into()),
            "locus" => out.locus = Some(v),
            "highlight" => out.highlight = Some(v),
            _ => {} // newer links may carry more; ignore what we don't know
        }
    }
    Ok(out)
}

/// Quotes an argument for the platform's shell.
fn quote(arg: &str) -> String {
    let safe = |c: char| c.is_ascii_alphanumeric() || "-._/:,=+@".contains(c) || (cfg!(windows) && c == '\\');
    if !arg.is_empty() && arg.chars().all(safe) {
        arg.to_owned()
    } else if cfg!(windows) {
        format!("\"{}\"", arg.replace('"', "\\\""))
    } else {
        format!("'{}'", arg.replace('\'', "'\\''"))
    }
}

/// The `ggv-desktop` command that reopens `s`.
pub fn to_command(s: &Session) -> String {
    let mut args = vec!["ggv-desktop".to_owned()];
    // Alignments and signal files are positional (sorted by extension on startup).
    args.extend(s.bams.iter().chain(&s.signals).map(|p| quote(&p.to_string_lossy())));
    for v in &s.vcfs {
        args.push("-v".into());
        args.push(quote(&v.to_string_lossy()));
    }
    for g in &s.genes {
        args.push("-g".into());
        args.push(quote(&g.to_string_lossy()));
    }
    if let Some(r) = &s.reference {
        args.push("-f".into());
        args.push(quote(&r.to_string_lossy()));
    }
    // The outlined base, if any, is the more useful place to come back to.
    if let Some(l) = s.highlight.as_ref().or(s.locus.as_ref()) {
        args.push("--locus".into());
        args.push(quote(l));
    }
    args.join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_round_trip_awkward_paths() {
        let s = Session {
            bams: vec!["C:\\data\\tumor, run 2.bam".into(), "/d/normal.bam".into()],
            vcfs: vec!["/d/calls&more.vcf.gz".into()],
            reference: Some("/ref/hg38 #1.fa".into()),
            locus: Some("chr17:43,044,295-43,125,483".into()),
            highlight: Some("chr17:43092919".into()),
            ..Default::default()
        };
        let link = to_link(&s);
        assert!(link.starts_with("ggv://open?file=C:%5Cdata%5Ctumor,%20run%202.bam&"), "{link}");
        let back = parse(&link).unwrap();
        assert_eq!(back.files, vec![PathBuf::from("C:\\data\\tumor, run 2.bam"), "/d/normal.bam".into(), "/d/calls&more.vcf.gz".into()]);
        assert_eq!(back.reference, Some("/ref/hg38 #1.fa".into()));
        assert_eq!(back.locus.as_deref(), Some("chr17:43,044,295-43,125,483"));
        assert_eq!(back.highlight.as_deref(), Some("chr17:43092919"));
    }

    #[test]
    fn every_track_kind_survives_a_link_and_a_command() {
        let s = Session {
            bams: vec!["/d/a.bam".into()],
            vcfs: vec!["/d/c.vcf.gz".into()],
            genes: vec!["/d/g.gtf".into()],
            signals: vec!["/d/s.bw".into(), "/d/t.bedGraph".into()],
            reference: Some("/d/r.fa".into()),
            ..Default::default()
        };
        let sorted = crate::files::sort(parse(&to_link(&s)).unwrap().files);
        assert_eq!((sorted.bams, sorted.vcfs, sorted.genes, sorted.signals), (s.bams.clone(), s.vcfs.clone(), s.genes.clone(), s.signals.clone()));
        let c = to_command(&s);
        assert!(c.contains("/d/s.bw") && c.contains("/d/t.bedGraph"), "{c}");
    }

    #[test]
    fn rejects_other_links_and_ignores_unknown_keys() {
        assert!(parse("http://x").is_err());
        assert!(parse("ggv://delete?file=x").is_err());
        let l = parse("GGV://open/?file=a.bam&future=1&locus=chr1%3A5").unwrap();
        assert_eq!(l.files, vec![PathBuf::from("a.bam")]);
        assert_eq!(l.locus.as_deref(), Some("chr1:5"));
        // Links from before the rename still open.
        assert!(is_link("gv://open?file=a.bam") && !is_link("ggv:/open") && !is_link("g"));
        assert_eq!(parse("GV://open?file=a.bam&locus=chr2%3A7").unwrap().locus.as_deref(), Some("chr2:7"));
        assert_eq!(decode("100%"), "100%");
    }

    #[test]
    fn command_quotes_what_needs_it() {
        let s = Session { bams: vec!["/d/a b.bam".into()], genes: vec!["/g.gtf".into()], locus: Some("chr1:1-100".into()), ..Default::default() };
        let c = to_command(&s);
        if cfg!(windows) {
            assert_eq!(c, "ggv-desktop \"/d/a b.bam\" -g /g.gtf --locus chr1:1-100");
        } else {
            assert_eq!(c, "ggv-desktop '/d/a b.bam' -g /g.gtf --locus chr1:1-100");
        }
    }
}
