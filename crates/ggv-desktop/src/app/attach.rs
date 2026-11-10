//! What a pane tells an attached agent (`ggv mcp --attach`), and the commands it takes.

use std::path::{Path, PathBuf};

use eframe::egui;
use serde_json::{Value, json};

use super::{App, Hit};

fn paths(v: &[PathBuf]) -> Vec<String> {
    v.iter().map(|p| p.display().to_string()).collect()
}

/// Index extensions GGV reads next to a data file (`x.bam.bai` or `x.bai`).
const INDEXES: &[&str] = &["bai", "csi", "crai", "tbi", "fai", "gzi"];

/// Checks an export destination before anything is written: it must end in one of `exts`, and
/// must not be (under any spelling or link) one of `inputs` or an index of one.
pub(super) fn check_export(path: &Path, exts: &[&str], inputs: &[PathBuf]) -> Result<(), String> {
    let name = path.file_name().map(|n| n.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
    if !exts.iter().any(|e| name.ends_with(&format!(".{e}"))) {
        let list = exts.iter().map(|e| format!(".{e}")).collect::<Vec<_>>().join(" or ");
        return Err(format!("{}: the file must end in {list}", path.display()));
    }
    // Only an existing file can be an open one.
    let Ok(dest) = std::fs::canonicalize(path) else { return Ok(()) };
    for input in inputs {
        let mut protected = vec![input.clone()];
        for e in INDEXES {
            protected.push(input.with_added_extension(e));
            protected.push(input.with_extension(e));
        }
        if protected.iter().any(|p| std::fs::canonicalize(p).is_ok_and(|c| c == dest)) {
            return Err(format!("{} is open in GGV (or is an index of {}): choose another file", path.display(), input.display()));
        }
    }
    Ok(())
}

impl App {
    /// This pane's state: view, outlined base, selection, files, each file's read filter and
    /// the display options. Positions are 1-based, regions inclusive (as typed in the app).
    pub fn attach_state(&self) -> Value {
        let s = self.session();
        let name = |c: usize| self.contigs.get(c).map(|x| x.name.clone()).unwrap_or_default();
        let view = self.contigs.get(self.contig).map(|c| {
            let (vs, ve) = self.visible();
            let (start, end) = (vs.max(0.0) as u64 + 1, (ve.min(self.contig_len()) as u64).max(vs.max(0.0) as u64 + 1));
            json!({"region": format!("{}:{start}-{end}", c.name), "contig": c.name, "start": start, "end": end, "centre": self.center.max(0.0) as u64 + 1})
        });
        let filters: Vec<Value> = s
            .bams
            .iter()
            .map(|p| {
                let f = self.track_style(p).filter;
                json!({
                    "path": p.display().to_string(),
                    "exclude_flags": format!("0x{:x}", f.exclude_flags()),
                    "min_mapq": f.min_mapq,
                    "expression": f.expr.map(crate::expr::text),
                })
            })
            .collect();
        json!({
            "view": view,
            "highlight": self.highlight.map(|(c, p)| format!("{}:{}", name(c), p + 1)),
            "selection": self.selection.map(|(c, a, b)| format!("{}:{}-{}", name(c), a + 1, b)),
            "files": {"bams": paths(&s.bams), "reference": s.reference.as_ref().map(|p| p.display().to_string()), "vcfs": paths(&s.vcfs), "genes": paths(&s.genes), "signals": paths(&s.signals)},
            "assembly": self.assembly_state(),
            "error": self.error,
            "warning": self.warning(),
            "notice": self.notice.as_ref().map(|(n, _)| n.clone()),
            "track_errors": self.tracks.iter().filter_map(|t| t.error.as_ref().map(|e| json!({"file": t.loader.bam_path.display().to_string(), "error": e}))).chain(self.variants.iter().filter_map(|v| v.error.as_ref().map(|e| json!({"file": v.path.display().to_string(), "error": e})))).collect::<Vec<_>>(),
            "filters": filters,
            "display": {
                "color": s.color, "group": s.group, "group_pos": s.group_pos, "sort": s.sort,
                "soft_clips": s.soft_clips, "squished": s.squished, "junctions": s.junctions, "all_transcripts": s.all_transcripts, "link_splits": s.link_splits,
                "read_window": super::read_window(),
                "variant_colors": self.variants.iter().map(|v| (v.file_name(), json!(v.color.name()))).collect::<serde_json::Map<_, _>>(),
                "filter": {
                    "min_mapq": self.style.filter.min_mapq,
                    "hide_duplicates": self.style.filter.hide_duplicates,
                    "hide_secondary": self.style.filter.hide_secondary,
                    "hide_supplementary": self.style.filter.hide_supplementary,
                    "hide_qc_fail": self.style.filter.hide_qc_fail,
                    "expression": self.style.filter.expr.map(crate::expr::text).unwrap_or_default(),
                },
            },
            "ready": self.is_ready(),
        })
    }

    /// Moves the view to a locus, region or gene (a single base is outlined, as typed in the
    /// locus box).
    pub fn attach_goto(&mut self, text: &str) -> Result<(), String> {
        if self.contigs.is_empty() {
            return Err("no files are open in this pane".into());
        }
        match self.resolve(text) {
            Err(Some(e)) => Err(e),
            // Genes still loading: the pane goes there once they have.
            Err(None) | Ok(_) => {
                self.goto(text);
                // A selection the new view doesn't show is no longer what the person looks at.
                let (vs, ve) = self.visible();
                if self.selection.is_some_and(|(c, a, b)| c != self.contig || (b as f64) < vs || (a as f64) > ve) {
                    self.selection = None;
                }
                Ok(())
            }
        }
    }

    /// Selects a region (shown, with its statistics in the Region window).
    pub fn attach_select(&mut self, text: &str) -> Result<(), String> {
        let region = self.resolve(text).map_err(|e| e.unwrap_or_else(|| format!("`{text}`: genes are still loading; try again")))?;
        if region.len() > 10_000_000 {
            return Err("select at most 10 Mb".into());
        }
        let c = self.contigs.iter().position(|x| x.name == region.contig).ok_or("unknown contig")?;
        let (vs, ve) = self.visible();
        if c != self.contig || (region.start as f64) < vs || (region.end as f64) > ve {
            let pad = (region.len() / 10).max(10);
            let shown = format!("{}:{}-{}", region.contig, region.start.saturating_sub(pad) + 1, region.end + pad);
            self.goto(&shown);
        }
        self.selection = Some((c, region.start, region.end));
        Ok(())
    }

    /// Changes how this pane shows reads: colour, grouping, sort, soft clips, squished rows,
    /// transcripts and read filters (for every file, or one file's filter expression). Every
    /// value is checked before any is applied, so a bad one changes nothing. Positions are
    /// 1-based; grouping and sorting at a base use `pos`, else the outlined base, else the view
    /// centre.
    pub fn attach_display(&mut self, p: &Value) -> Result<(), String> {
        use crate::session::{parse_color, parse_group, parse_sort};
        let text = |k: &str| p.get(k).and_then(Value::as_str).map(str::trim);
        let flag = |k: &str| -> Result<Option<bool>, String> {
            match p.get(k) {
                None | Some(Value::Null) => Ok(None),
                Some(Value::Bool(b)) => Ok(Some(*b)),
                Some(v) => Err(format!("`{k}` must be true or false, not {v}")),
            }
        };
        // The CLI's spellings are accepted too.
        let norm = |s: &str| match s {
            "insert-size" => "insert size".to_owned(),
            "base-mods" | "methylation" => "base modifications".to_owned(),
            "read-group" => "read group".to_owned(),
            "mapq" => "MAPQ".to_owned(),
            other => other.to_owned(),
        };
        let base = match p.get("pos").and_then(Value::as_u64) {
            Some(0) => return Err("`pos` is 1-based".into()),
            Some(pos) if pos as f64 > self.contig_len() => {
                let name = self.contigs.get(self.contig).map_or("", |c| c.name.as_str());
                return Err(format!("`pos` {pos} is past the end of {name} ({} bp)", self.contig_len() as u64));
            }
            Some(pos) => pos - 1,
            None => match self.highlight {
                Some((c, pos)) if c == self.contig => pos,
                _ => self.center.max(0.0) as u64,
            },
        };

        let mut style = self.style;
        if let Some(c) = text("color") {
            let name = norm(c);
            let parsed = parse_color(&name);
            if name != "none" && crate::session::color_text(parsed) != name {
                return Err(format!("color: `{c}` is not one of none, strand, haplotype, insert size, base modifications, tag:XX"));
            }
            style.color = parsed;
        }
        if let Some(g) = text("group") {
            let name = norm(g);
            // `group_pos` (1-based) is what the state reports and `before` hands back.
            let at = match p.get("group_pos").and_then(Value::as_u64) {
                Some(0) => return Err("`group_pos` is 1-based".into()),
                Some(gp) => gp,
                None => base + 1,
            };
            let parsed = parse_group(&name, Some(at));
            if name != "none" && crate::session::group_text(parsed) != name {
                return Err(format!("group: `{g}` is not one of none, strand, haplotype, read group, alleles, base, tag:XX"));
            }
            style.group = parsed;
        }
        // As the state reports it (and `before` hands back): {by, pos} or null for none.
        match p.get("sort") {
            Some(Value::Object(o)) => {
                let by = o.get("by").and_then(Value::as_str).ok_or("sort.by: give how to sort (base, strand, MAPQ, start, insert size)")?;
                let pos = o.get("pos").and_then(Value::as_u64).filter(|&x| x >= 1).ok_or("sort.pos: a 1-based position")?;
                style.sort = Some((pos - 1, parse_sort(&norm(by)).ok_or_else(|| format!("sort.by: `{by}` is not one of base, strand, MAPQ, start, insert size"))?));
            }
            Some(Value::Null) => style.sort = None,
            _ => {}
        }
        if let Some(s) = text("sort") {
            let name = norm(s);
            style.sort = if name == "none" {
                None
            } else {
                Some((base, parse_sort(&name).ok_or_else(|| format!("sort: `{s}` is not one of none, base, strand, MAPQ, start, insert size"))?))
            };
        }
        if let Some(b) = flag("link_splits")? {
            style.link_splits = b;
        }
        if let Some(b) = flag("soft_clips")? {
            style.soft_clips = b;
        }
        if let Some(b) = flag("squished")? {
            style.squished = b;
        }
        if let Some(b) = flag("junctions")? {
            style.junctions = b;
        }
        let all_transcripts = flag("all_transcripts")?;
        let read_window = match p.get("read_window") {
            None | Some(Value::Null) => None,
            Some(v) => Some(
                v.as_u64()
                    .filter(|bp| (super::MIN_READ_WINDOW..=super::MAX_READ_WINDOW).contains(bp))
                    .ok_or_else(|| format!("read_window is {}-{} bases", super::MIN_READ_WINDOW, super::MAX_READ_WINDOW))?,
            ),
        };
        let mut shared_expr: Option<String> = None;
        if let Some(f) = p.get("filter").filter(|v| !v.is_null()) {
            let f = f.as_object().ok_or("`filter` is an object: min_mapq, hide_duplicates, hide_secondary, hide_supplementary, hide_qc_fail, expression")?;
            for (k, v) in f {
                match (k.as_str(), v) {
                    ("min_mapq", v) => style.filter.min_mapq = v.as_u64().filter(|&q| q <= 255).ok_or("filter.min_mapq is 0-255")? as u8,
                    ("hide_duplicates", Value::Bool(b)) => style.filter.hide_duplicates = *b,
                    ("hide_secondary", Value::Bool(b)) => style.filter.hide_secondary = *b,
                    ("hide_supplementary", Value::Bool(b)) => style.filter.hide_supplementary = *b,
                    ("hide_qc_fail", Value::Bool(b)) => style.filter.hide_qc_fail = *b,
                    ("expression", Value::String(e)) => {
                        if !e.trim().is_empty() {
                            crate::expr::intern(e).map_err(|err| format!("filter.expression: {err}"))?;
                        }
                        shared_expr = Some(e.clone());
                    }
                    (k, v) => return Err(format!("filter.{k}: unknown setting or wrong type ({v})")),
                }
            }
        }
        // One file's own filter expression.
        let file = match (text("file"), p.get("expression").and_then(Value::as_str)) {
            (Some(f), Some(e)) => {
                let paths: Vec<std::path::PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
                let path = unique_file(&paths, f).map_err(|e| format!("file: {e}"))?;
                if !e.trim().is_empty() {
                    crate::expr::intern(e).map_err(|err| format!("expression: {err}"))?;
                }
                Some((path, e.to_owned()))
            }
            (None, None) => None,
            _ => return Err("give `file` and `expression` together (an empty expression clears that file's filter)".into()),
        };

        // Variant site colouring: `af` or an INFO key, for one VCF (`variant_file`) or for every
        // one whose header declares the key.
        let vcf_paths: Vec<std::path::PathBuf> = self.variants.iter().map(|v| v.path.clone()).collect();
        // Per-file colours as the state reports them (`before` hands these back).
        let mut variant_colors = Vec::new();
        if let Some(m) = p.get("variant_colors").filter(|v| !v.is_null()) {
            let m = m.as_object().ok_or("`variant_colors` maps a VCF (path or file name) to `af` or an INFO key")?;
            for (f, c) in m {
                let c = c.as_str().ok_or_else(|| format!("variant_colors[{f}]: give `af` or an INFO key"))?;
                let path = unique_file(&vcf_paths, f).map_err(|e| format!("variant_colors: {e}"))?;
                let i = self.variants.iter().position(|v| v.path == path).expect("listed above");
                variant_colors.push((i, crate::variants::SiteColor::parse(c)));
            }
        }
        let variant_color = match text("variant_color") {
            Some(c) => {
                let want = crate::variants::SiteColor::parse(c);
                let named = match text("variant_file") {
                    Some(f) => Some(unique_file(&vcf_paths, f).map_err(|e| format!("variant_file: {e}"))?),
                    None => None,
                };
                let mut hits = Vec::new();
                for (i, v) in self.variants.iter().enumerate() {
                    if named.as_ref().is_some_and(|f| &v.path != f) {
                        continue;
                    }
                    let known = match &want {
                        crate::variants::SiteColor::AlleleFrequency => true,
                        crate::variants::SiteColor::Info(k) => v.info_fields.iter().any(|f| &f.id == k),
                    };
                    if known {
                        hits.push(i);
                    } else if named.is_some() {
                        return Err(format!("variant_color: {} declares no INFO/{c}", v.file_name()));
                    }
                }
                if hits.is_empty() {
                    return Err(match named {
                        Some(f) => format!("variant_file: {} declares no INFO/{c}", crate::files::display_name(&f)),
                        None => format!("variant_color: no open VCF declares INFO/{c}"),
                    });
                }
                Some((want, hits))
            }
            None => None,
        };

        // Everything checked: apply.
        self.style = style;
        for (i, c) in variant_colors {
            self.variants[i].color = c;
        }
        if let Some((c, hits)) = variant_color {
            for i in hits {
                self.variants[i].color = c.clone();
            }
        }
        if let Some(e) = shared_expr {
            self.apply_filter_text(&e, None);
        }
        if let Some((path, e)) = file {
            self.apply_filter_text(&e, Some(path));
        }
        if let Some(b) = all_transcripts {
            self.all_transcripts = b;
        }
        if let Some(bp) = read_window {
            super::save_read_window(bp);
        }
        Ok(())
    }

    /// Opens `files` in this pane. With `replace`, the alignment, variant, gene and signal
    /// tracks not among them close first (the reference stays unless a FASTA is given), and
    /// those already open stay as they are. A missing or unknown file changes nothing.
    /// Clears the selection.
    pub fn attach_unselect(&mut self) {
        self.selection = None;
    }

    /// The window's bookmarks: `list` (default), `add` (at `locus`, else the selection, the
    /// outlined base or the view, with an optional `note`), `remove` (by 1-based `index` or
    /// `locus`), `clear`, or `export` to a BED file at `path`. Returns the bookmarks after.
    pub fn attach_bookmarks(&mut self, p: &Value) -> Result<Value, String> {
        let text = |k: &str| p.get(k).and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty());
        let place = |this: &Self, t: &str| -> Result<(String, u64, u64), String> {
            let r = this.resolve(t).map_err(|e| e.unwrap_or_else(|| format!("`{t}`: genes are still loading; try again")))?;
            Ok((r.contig, r.start, r.end))
        };
        match text("action").unwrap_or("list") {
            "list" => {}
            "add" => {
                let (contig, start, end) = match text("locus") {
                    Some(t) => place(self, t)?,
                    None => self.here()?,
                };
                let note = text("note").unwrap_or("").to_owned();
                // Shown as its band on the ruler: the Bookmarks window stays as the person has
                // it (opening it for every bookmark an agent adds covered the panels).
                crate::bookmarks::add(&mut self.bookmarks.borrow_mut(), crate::bookmarks::Bookmark { contig, start, end, note });
            }
            "remove" => {
                let mut list = self.bookmarks.borrow_mut();
                let i = match (p.get("index").and_then(Value::as_u64), text("locus")) {
                    (Some(n), _) => (n as usize).checked_sub(1).filter(|&i| i < list.len()).ok_or_else(|| format!("no bookmark {n}: there are {}", list.len()))?,
                    (None, Some(t)) => {
                        let (contig, start, end) = place(self, t)?;
                        list.iter().position(|b| b.contig == contig && b.start == start && b.end == end).ok_or_else(|| format!("no bookmark at {t}"))?
                    }
                    (None, None) => return Err("give the bookmark's `index` or `locus`".into()),
                };
                list.remove(i);
            }
            "clear" => self.bookmarks.borrow_mut().clear(),
            "export" => {
                let path = text("path").ok_or("give a `path` for the BED file")?;
                check_export(Path::new(path), &["bed"], &self.protected_files())?;
                crate::store::write_atomic(Path::new(path), crate::bookmarks::to_bed(&self.bookmarks.borrow()).as_bytes()).map_err(|e| format!("{path}: {e}"))?;
            }
            other => return Err(format!("unknown action `{other}`: list, add, remove, clear or export")),
        }
        let list: Vec<Value> = self
            .bookmarks
            .borrow()
            .iter()
            .enumerate()
            .map(|(i, b)| json!({"index": i + 1, "locus": b.locus(), "note": b.note}))
            .collect();
        Ok(json!({"bookmarks": list}))
    }

    /// The window's overlays: `list` (default); `add` one over `region` (a base, region or
    /// gene; else the selection, the outlined base or the view) with a `label` and `color`;
    /// `remove` one by 1-based `index`; or `clear` them. Returns the overlays after.
    pub fn attach_overlays(&mut self, p: &Value) -> Result<Value, String> {
        let text = |k: &str| p.get(k).and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty());
        match text("action").unwrap_or("list") {
            "list" => {}
            "add" => {
                let (contig, start, end) = match text("region") {
                    Some(t) => {
                        let r = self.resolve(t).map_err(|e| e.unwrap_or_else(|| format!("`{t}`: genes are still loading; try again")))?;
                        (r.contig, r.start, r.end)
                    }
                    None => self.here()?,
                };
                let color = text("color").unwrap_or("yellow");
                if crate::overlays::rgb(color).is_none() {
                    let names: Vec<&str> = crate::overlays::COLORS.iter().map(|(n, _)| *n).collect();
                    return Err(format!("color: `{color}` is not one of {}", names.join(", ")));
                }
                let label = text("label").unwrap_or("").chars().take(120).collect();
                self.overlays.borrow_mut().push(crate::overlays::Overlay { contig, start, end, label, color: color.to_owned() });
            }
            "remove" => {
                let mut list = self.overlays.borrow_mut();
                let n = p.get("index").and_then(Value::as_u64).ok_or("give the overlay's `index` (1-based)")?;
                let i = (n as usize).checked_sub(1).filter(|&i| i < list.len()).ok_or_else(|| format!("no overlay {n}: there are {}", list.len()))?;
                list.remove(i);
            }
            "clear" => self.overlays.borrow_mut().clear(),
            other => return Err(format!("unknown action `{other}`: list, add, remove or clear")),
        }
        let list: Vec<Value> = self
            .overlays
            .borrow()
            .iter()
            .enumerate()
            .map(|(i, o)| json!({"index": i + 1, "region": o.locus(), "label": o.label, "color": o.color}))
            .collect();
        Ok(json!({"overlays": list}))
    }

    /// Where a bookmark or overlay given no place goes: the selection, the outlined base or the
    /// view, as Ctrl+D does; refused with nothing open (there is no contig to name).
    fn here(&self) -> Result<(String, u64, u64), String> {
        let (c, start, end) = self.place().ok_or("no files are open in this pane")?;
        Ok((self.contigs[c].name.clone(), start, end))
    }

    /// Every file open in the window that an export must not overwrite: this pane's (see
    /// `own_files`) and the other panes' (their queue lists and decisions too).
    pub(super) fn protected_files(&self) -> Vec<PathBuf> {
        let mut v = self.own_files();
        v.extend(self.other_panes_files.iter().cloned());
        v
    }

    /// Every file open here that an export must not overwrite: tracks, the reference, the
    /// queue's list, and the review and queue decision files with their backups (every open
    /// VCF's, wherever saved: a review not running now still has its decisions).
    pub(crate) fn own_files(&self) -> Vec<PathBuf> {
        let mut v: Vec<PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        v.extend(self.variants.iter().map(|x| x.path.clone()));
        v.extend(self.genes.iter().map(|x| x.path.clone()));
        v.extend(self.signals.iter().map(|x| x.path.clone()));
        v.extend(self.reference_path.clone());
        let mut decisions: Vec<PathBuf> = self.variants.iter().flat_map(|x| crate::store::sidecars(&x.path, crate::review::EXT)).collect();
        if let Some(q) = &self.queue {
            v.extend(q.list.clone());
            decisions.extend(q.save_path.clone());
        }
        decisions.extend(self.review.as_ref().map(|r| r.path.clone()));
        v.extend(decisions.iter().map(|d| crate::store::backup_of(d)));
        v.extend(decisions);
        v
    }

    /// Saves the session (files, view, display, filters, bookmarks, track layout) to `path`.
    pub fn attach_save_session(&self, p: &Value) -> Result<Value, String> {
        let path = p.get("path").and_then(Value::as_str).filter(|s| !s.trim().is_empty()).ok_or("give a `path`, e.g. review.ggv.json")?;
        let path = std::path::Path::new(path);
        if crate::files::kind_of(&path.to_string_lossy()) != Some(crate::files::Kind::Session) {
            return Err("a session file ends in .json (e.g. review.ggv.json)".into());
        }
        check_export(path, &["json"], &self.protected_files())?;
        self.session().with_workspace().save(path).map_err(|e| format!("{e:#}"))?;
        Ok(json!({"saved": path.display().to_string()}))
    }

    /// The open tracks by kind, in order; `close` one or `move` it `up`/`down` within its kind.
    /// A track is named by its path, or its file name when that is unique.
    pub fn attach_tracks(&mut self, p: &Value) -> Result<Value, String> {
        let s = self.session();
        type Kind<'a> = (&'a str, &'a [PathBuf], fn(usize) -> Hit);
        let kinds: [Kind; 4] =
            [("alignments", &s.bams, Hit::Track), ("variants", &s.vcfs, Hit::Variants), ("genes", &s.genes, Hit::Genes), ("signals", &s.signals, Hit::Signal)];
        let action = p.get("action").and_then(Value::as_str).unwrap_or("list");
        if action != "list" {
            let want = p.get("track").and_then(Value::as_str).map(str::trim).filter(|t| !t.is_empty()).ok_or("give the `track` (path or file name)")?;
            let hits: Vec<(Hit, &str)> = kinds
                .iter()
                .flat_map(|(kind, paths, hit)| {
                    paths.iter().enumerate().filter(|(_, path)| crate::files::names(path, want)).map(move |(i, _)| (hit(i), *kind))
                })
                .collect();
            let hit = match hits.as_slice() {
                [(h, _)] => *h,
                [] => return Err(format!("no open track is `{want}`")),
                _ => return Err(format!("`{want}` names {} tracks; give the full path", hits.len())),
            };
            match action {
                "close" => self.close_hit(hit),
                "move" => match p.get("direction").and_then(Value::as_str) {
                    Some("up") => self.move_track(hit, true),
                    Some("down") => self.move_track(hit, false),
                    _ => return Err("give `direction`: up or down".into()),
                },
                other => return Err(format!("unknown action `{other}`: list, close or move")),
            }
        }
        let s = self.session();
        let list = |paths: &[PathBuf]| paths.iter().map(|p| json!({"file": p.file_name().map(|n| n.to_string_lossy().into_owned()), "path": p.display().to_string()})).collect::<Vec<_>>();
        Ok(json!({
            "alignments": list(&s.bams), "variants": list(&s.vcfs), "genes": list(&s.genes), "signals": list(&s.signals),
            "reference": s.reference.as_ref().map(|r| r.display().to_string()),
        }))
    }

    /// Variant triage decisions for each open VCF (saved next to it as `.ggv-review.tsv`), and
    /// which one is being reviewed.
    pub fn attach_review(&self) -> Value {
        let files: Vec<Value> = self
            .variants
            .iter()
            .map(|v| {
                let r = crate::review::Review::open(v.id, &v.path);
                let [real, artifact, unsure] = r.counts();
                let decisions: Vec<Value> = r
                    .decisions
                    .iter()
                    .map(|(k, d)| json!({"locus": format!("{}:{}", k.contig, k.pos), "ref": k.reference, "alt": k.alt, "decision": d.label.name(), "note": d.note}))
                    .collect();
                json!({
                    "vcf": v.path.display().to_string(),
                    "file": r.path.display().to_string(),
                    "reviewing": self.review.as_ref().is_some_and(|x| x.track == v.id),
                    "counts": {"real": real, "artifact": artifact, "unsure": unsure},
                    "decisions": decisions,
                    "error": r.error, "notice": r.notice,
                })
            })
            .collect();
        json!({"review": files})
    }

    pub fn attach_open(&mut self, files: Vec<PathBuf>, replace: bool, ctx: &egui::Context) -> Result<(), String> {
        if let Some(missing) = files.iter().find(|p| !p.exists()) {
            return Err(format!("{}: no such file", missing.display()));
        }
        if let Some(unknown) = files.iter().find(|p| crate::files::kind_of(&p.to_string_lossy()).is_none()) {
            return Err(format!("{}: not a BAM/CRAM, VCF/BCF, GTF/GFF/BED, FASTA, bigWig/bedGraph/SEG or session", unknown.display()));
        }
        // A session replaces the workspace; one that can't be read is refused, leaving the
        // workspace as it was. Other files given with it are added to it.
        let is_session = |p: &PathBuf| crate::files::kind_of(&p.to_string_lossy()) == Some(crate::files::Kind::Session);
        let (sessions, files): (Vec<PathBuf>, Vec<PathBuf>) = files.into_iter().partition(is_session);
        if sessions.len() > 1 {
            return Err("open one session at a time".into());
        }
        if let Some(session) = sessions.first() {
            self.open_session(session, ctx)?;
            if !files.is_empty() {
                self.open_files(files, ctx);
            }
            return Ok(());
        }
        // Paths are compared as the files they name: an agent may give a relative path or a
        // symlink for a file already open.
        let same = |p: &Path| std::fs::canonicalize(p).unwrap_or_else(|_| p.to_owned());
        let wanted: Vec<PathBuf> = files.iter().map(|p| same(p)).collect();
        let given = |p: &PathBuf| wanted.contains(&same(p));
        // The reference already open isn't reopened (that would reload every track and move the view).
        let s = self.session();
        let reference = s.reference.as_deref().map(same);
        let files: Vec<PathBuf> = files.into_iter().filter(|p| reference.as_ref() != Some(&same(p))).collect();
        let files = if replace {
            for (i, p) in s.bams.iter().enumerate().rev() {
                if !given(p) {
                    self.close_track(i);
                }
            }
            for (i, p) in s.vcfs.iter().enumerate().rev() {
                if !given(p) {
                    self.close_variants(i);
                }
            }
            for (i, p) in s.genes.iter().enumerate().rev() {
                if !given(p) {
                    self.close_genes(i);
                }
            }
            for (i, p) in s.signals.iter().enumerate().rev() {
                if !given(p) {
                    self.close_signal(i);
                }
            }
            let open: Vec<PathBuf> = [s.bams, s.vcfs, s.genes, s.signals].concat().iter().map(|p| same(p)).collect();
            files.into_iter().filter(|p| !open.contains(&same(p))).collect()
        } else {
            files
        };
        self.open_files(files, ctx);
        Ok(())
    }
}

/// The one open file `name` means: its path (in any spelling), or a file name only one open
/// file has. Two open files called `same.bam` in different folders need the path.
pub(super) fn unique_file(paths: &[std::path::PathBuf], name: &str) -> Result<std::path::PathBuf, String> {
    if let Some(p) = paths.iter().find(|p| crate::files::same_file(p, Path::new(name))) {
        return Ok(p.clone());
    }
    let named: Vec<&std::path::PathBuf> = paths.iter().filter(|p| crate::files::display_name(p) == name).collect();
    match named.as_slice() {
        [one] => Ok((*one).clone()),
        [] => Err(format!("`{name}` is not open in this pane")),
        many => Err(format!(
            "`{name}` is ambiguous: {} open files have that name ({}); give the full path",
            many.len(),
            many.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::{check_export, unique_file};
    use std::path::PathBuf;

    #[test]
    fn exports_never_overwrite_open_files() {
        let dir = std::env::temp_dir().join(format!("ggv-export-guard-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        let bam = dir.join("reads.bam");
        let fa = dir.join("ref.fa.gz");
        for (f, text) in [("reads.bam", "BAM"), ("reads.bam.bai", "BAI"), ("ref.fa.gz", "FA"), ("ref.fa.gz.fai", "FAI"), ("ref.fa.gz.gzi", "GZI"), ("existing.json", "{\"keep\":\"original\"}")] {
            std::fs::write(dir.join(f), text).unwrap();
        }
        let open = vec![bam.clone(), fa.clone()];
        // The inputs, their indexes and another spelling of them are refused.
        for p in [bam.clone(), dir.join("reads.bam.bai"), fa.clone(), dir.join("ref.fa.gz.fai"), dir.join("ref.fa.gz.gzi"), dir.join("sub").join("..").join("reads.bam")] {
            let as_bed = p.with_added_extension("x");
            assert!(check_export(&p, &["bam", "bai", "gz", "fai", "gzi"], &open).unwrap_err().contains("open in GGV"), "{}", p.display());
            assert!(check_export(&as_bed, &["x"], &open).is_ok(), "a new file beside it is fine");
        }
        // Unsupported extensions are refused before anything is written.
        assert!(check_export(&dir.join("existing.json"), &["bed", "tsv"], &open).unwrap_err().contains(".bed or .tsv"));
        assert!(check_export(&bam, &["bed"], &open).is_err());
        assert!(check_export(&dir.join("out.BED"), &["bed"], &open).is_ok());
        assert!(check_export(&dir.join("out.tsv"), &["bed", "tsv"], &open).is_ok());
        assert_eq!(std::fs::read_to_string(&bam).unwrap(), "BAM");
        assert_eq!(std::fs::read_to_string(dir.join("existing.json")).unwrap(), "{\"keep\":\"original\"}");
        #[cfg(unix)]
        {
            let link = dir.join("alias.bam");
            std::os::unix::fs::symlink(&bam, &link).unwrap();
            assert!(check_export(&link, &["bam"], &open).is_err(), "a link to an open file");
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn file_names_must_be_unique() {
        let paths = vec![PathBuf::from("/data/same-a/same.bam"), PathBuf::from("/data/same-b/same.bam"), PathBuf::from("/data/other.bam")];
        assert!(unique_file(&paths, "same.bam").unwrap_err().contains("ambiguous"));
        assert_eq!(unique_file(&paths, "/data/same-b/same.bam").unwrap(), paths[1]);
        assert_eq!(unique_file(&paths, "other.bam").unwrap(), paths[2]);
        assert!(unique_file(&paths, "missing.bam").is_err());
    }
}
