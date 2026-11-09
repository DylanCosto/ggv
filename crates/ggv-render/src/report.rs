//! Evidence reports: one self-contained HTML file for a variant or region, with the view as an
//! embedded PNG, allele counts, strand and mapping-quality statistics, and review notes. It
//! opens in any browser and can be attached or printed to PDF.

#[derive(Default)]
pub struct Report {
    pub title: String,
    pub subtitle: String,
    pub png: Option<Vec<u8>>,
    pub sections: Vec<Section>,
}

#[derive(Default)]
pub struct Section {
    pub heading: String,
    pub lines: Vec<String>,
    /// Header row, then rows.
    pub table: Vec<Vec<String>>,
    /// A line to highlight as a warning (e.g. one-strand support).
    pub warning: Option<String>,
}

pub fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            c => out.push(c),
        }
    }
    out
}

/// A relative file name as a URL path: every byte but unreserved ones and `/` percent-encoded,
/// so a `#`, `%` or `?` in a name is part of the path, not a fragment or query.
pub fn url_path(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b"-._~/".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

pub fn base64(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for c in data.chunks(3) {
        let n = (u32::from(c[0]) << 16) | (u32::from(*c.get(1).unwrap_or(&0)) << 8) | u32::from(*c.get(2).unwrap_or(&0));
        out.push(T[(n >> 18) as usize & 63] as char);
        out.push(T[(n >> 12) as usize & 63] as char);
        out.push(if c.len() > 1 { T[(n >> 6) as usize & 63] as char } else { '=' });
        out.push(if c.len() > 2 { T[n as usize & 63] as char } else { '=' });
    }
    out
}

/// `YYYY-MM-DD HH:MM UTC` for a Unix time (civil-from-days, no time zone database needed).
pub fn utc_stamp(secs: u64) -> String {
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02} {:02}:{:02} UTC", rem / 3600, rem / 60 % 60)
}

impl Report {
    pub fn html(&self) -> String {
        let mut h = String::new();
        h.push_str("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">");
        h.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">");
        h.push_str(&format!("<title>{}</title>", escape(&self.title)));
        h.push_str(
            "<style>
:root{--bg:#fff;--fg:#1d1d1f;--muted:#6b6b70;--rule:#e2e2e6;--warn:#b3261e;--head:#f5f5f7}
@media (prefers-color-scheme:dark){:root{--bg:#18181a;--fg:#e8e8ea;--muted:#a0a0a6;--rule:#333338;--warn:#ff8a80;--head:#232326}}
body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:0 auto;max-width:1100px;padding:24px 16px}
h1{font-size:22px;margin:0 0 2px}h2{font-size:16px;margin:28px 0 6px;border-bottom:1px solid var(--rule);padding-bottom:4px}
.sub{color:var(--muted);margin-bottom:16px}img{max-width:100%;border:1px solid var(--rule);border-radius:4px}
table{border-collapse:collapse;font:13px ui-monospace,Consolas,monospace;margin:6px 0}
th,td{padding:3px 10px;text-align:right;border-bottom:1px solid var(--rule)}th{background:var(--head);font-weight:600}
td:first-child,th:first-child{text-align:left}p{margin:4px 0}.warn{color:var(--warn);font-weight:600}
footer{color:var(--muted);font-size:12px;margin-top:32px}
</style></head><body>\n",
        );
        h.push_str(&format!("<h1>{}</h1><div class=\"sub\">{}</div>\n", escape(&self.title), escape(&self.subtitle)));
        if let Some(png) = &self.png {
            h.push_str(&format!("<img alt=\"GGV view of {}\" src=\"data:image/png;base64,{}\">\n", escape(&self.title), base64(png)));
        }
        for s in &self.sections {
            h.push_str(&format!("<h2>{}</h2>\n", escape(&s.heading)));
            for l in &s.lines {
                h.push_str(&format!("<p>{}</p>\n", escape(l)));
            }
            if let Some(w) = &s.warning {
                h.push_str(&format!("<p class=\"warn\">{}</p>\n", escape(w)));
            }
            if let Some((head, rows)) = s.table.split_first() {
                h.push_str("<table><tr>");
                for c in head {
                    h.push_str(&format!("<th>{}</th>", escape(c)));
                }
                h.push_str("</tr>\n");
                for r in rows {
                    h.push_str("<tr>");
                    for c in r {
                        h.push_str(&format!("<td>{}</td>", escape(c)));
                    }
                    h.push_str("</tr>\n");
                }
                h.push_str("</table>\n");
            }
        }
        h.push_str("<footer>Made with GGV</footer></body></html>\n");
        h
    }
}

/// One place on a contact sheet: a thumbnail linking to its report.
pub struct Card {
    /// Image file, relative to the sheet.
    pub image: String,
    /// Report file, relative to the sheet.
    pub link: String,
    pub title: String,
    pub lines: Vec<String>,
}

/// An `index.html` of cards: every place at a glance, each opening its report.
pub fn contact_sheet(title: &str, subtitle: &str, cards: &[Card]) -> String {
    let mut h = String::new();
    h.push_str("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">");
    h.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">");
    h.push_str(&format!("<title>{}</title>", escape(title)));
    h.push_str(
        "<style>
:root{--bg:#fff;--fg:#1d1d1f;--muted:#6b6b70;--rule:#e2e2e6;--card:#fafafb}
@media (prefers-color-scheme:dark){:root{--bg:#18181a;--fg:#e8e8ea;--muted:#a0a0a6;--rule:#333338;--card:#202023}}
body{background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif;margin:0 auto;max-width:1400px;padding:24px 16px}
h1{font-size:22px;margin:0 0 2px}.sub{color:var(--muted);margin-bottom:18px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:14px}
.card{display:block;background:var(--card);border:1px solid var(--rule);border-radius:6px;padding:8px;color:inherit;text-decoration:none}
a.card:hover{border-color:var(--muted)}.card img{width:100%;height:170px;object-fit:cover;object-position:top;border-radius:3px;display:block}
.t{font-weight:600;margin:6px 0 2px;font-family:ui-monospace,Consolas,monospace}.l{color:var(--muted);font-size:13px;overflow-wrap:anywhere}
</style></head><body>\n",
    );
    h.push_str(&format!("<h1>{}</h1><div class=\"sub\">{}</div>\n<div class=\"grid\">\n", escape(title), escape(subtitle)));
    for c in cards {
        // A failed place has no report: its card is no link (`href=""` reopened the sheet).
        let tag = if c.link.is_empty() { "div" } else { "a" };
        if c.link.is_empty() {
            h.push_str("<div class=\"card\">");
        } else {
            h.push_str(&format!("<a class=\"card\" href=\"{}\">", escape(&url_path(&c.link))));
        }
        if !c.image.is_empty() {
            h.push_str(&format!("<img loading=\"lazy\" alt=\"{}\" src=\"{}\">", escape(&c.title), escape(&url_path(&c.image))));
        }
        h.push_str(&format!("<div class=\"t\">{}</div>", escape(&c.title)));
        for l in &c.lines {
            h.push_str(&format!("<div class=\"l\">{}</div>", escape(l)));
        }
        h.push_str(&format!("</{tag}>\n"));
    }
    h.push_str("</div></body></html>\n");
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_matches_rfc_4648() {
        assert_eq!(base64(b""), "");
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
        assert_eq!(base64(b"foobar"), "Zm9vYmFy");
    }

    #[test]
    fn dates_and_escaping() {
        assert_eq!(utc_stamp(0), "1970-01-01 00:00 UTC");
        assert_eq!(utc_stamp(1_793_750_400 + 3_600 * 13 + 60 * 7), "2026-11-04 13:07 UTC");
        assert_eq!(utc_stamp(951_782_400), "2000-02-29 00:00 UTC");
        let r = Report {
            title: "chr1:100 <A>".into(),
            sections: vec![Section { heading: "x".into(), table: vec![vec!["a".into()], vec!["1&2".into()]], ..Default::default() }],
            ..Default::default()
        };
        let html = r.html();
        assert!(html.contains("<h1>chr1:100 &lt;A&gt;</h1>") && html.contains("<td>1&amp;2</td>"));
    }

    #[test]
    fn contact_sheet_links_are_urls_and_failed_cards_are_not_links() {
        // A PanSN contig's `#` started the URL fragment; a failed card linked to `href=""`.
        let card = |image: &str, link: &str| Card { image: image.into(), link: link.into(), title: "t".into(), lines: vec![] };
        let html = contact_sheet("x", "y", &[card("001-HG002#1#chrQ-15.png", "001-HG002#1#chrQ-15.html"), card("", "")]);
        assert!(html.contains("href=\"001-HG002%231%23chrQ-15.html\"") && html.contains("src=\"001-HG002%231%23chrQ-15.png\""), "{html}");
        assert!(!html.contains("href=\"\"") && html.contains("<div class=\"card\"><div class=\"t\">t</div></div>"), "{html}");
        assert_eq!(url_path("a b%?.html"), "a%20b%25%3F.html");
    }
}
