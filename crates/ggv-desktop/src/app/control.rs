//! What attached agents can do in a pane beyond moving the view (L07): run any native action by
//! name, operate the variant review and the review queue, run and stop a sequence search, and
//! set up the visible comparison. Every command reuses the code the window's own buttons and
//! keys run, so the UI, the palette and agents can't disagree, and each returns what it changed
//! or an error that says why nothing did.

use eframe::egui;
use serde_json::{Value, json};

use super::{App, actions::Action};
use crate::{
    queue::Verdict,
    review::Label,
    session::{color_name, group_name},
    tracks::CoverageScale,
};

/// The version of the attached command set. Bumped when a command or its arguments change in a
/// way an agent would notice; `capabilities` lists what this window understands.
pub const API_VERSION: u32 = 2;

/// Commands this window understands, for agents to check before relying on one.
pub const CAPABILITIES: &[&str] = &[
    "goto", "select", "unselect", "open", "display", "bookmarks", "overlays", "save_session", "tracks", "review", "queue", "action", "panes",
    "search", "compare", "rows",
];

/// Actions that open a native file dialog, which an agent can't answer: the dedicated commands
/// (`open`, `session`, `queue`, `bookmarks`) do these with a path instead.
fn needs_a_person(a: Action) -> bool {
    use Action::*;
    matches!(a, AddBam | AddVcf | AddGenes | SetReference | SaveSession | OpenSession | ExportPng | Report | ExportBookmarks | QueueOpen | Palette | GoToLocus)
}

fn snake(debug: &str) -> String {
    let mut out = String::new();
    for (i, ch) in debug.chars().enumerate() {
        if ch.is_ascii_uppercase() {
            if i > 0 {
                out.push('_');
            }
            out.push(ch.to_ascii_lowercase());
        } else {
            out.push(ch);
        }
    }
    out
}

impl Action {
    /// A stable machine name: `next_variant`, `color:strand`, `group:read group`, `scale:log`.
    pub fn id(self) -> String {
        match self {
            Action::Color(c) => format!("color:{}", color_name(c)),
            Action::Group(g) => format!("group:{}", group_name(g)),
            Action::Scale(s) => format!(
                "scale:{}",
                match s {
                    CoverageScale::View => "view",
                    CoverageScale::Loaded => "loaded",
                    CoverageScale::Log => "log",
                }
            ),
            a => snake(&format!("{a:?}")),
        }
    }

    pub fn by_id(id: &str) -> Option<Self> {
        Action::all().into_iter().find(|a| a.id() == id)
    }
}

fn label_named(s: &str) -> Option<Label> {
    [Label::Real, Label::Artifact, Label::Unsure].into_iter().find(|l| l.name().eq_ignore_ascii_case(s))
}

fn verdict_named(s: &str) -> Option<Verdict> {
    [Verdict::Pass, Verdict::Fail, Verdict::Unsure].into_iter().find(|v| v.name().eq_ignore_ascii_case(s))
}

fn text<'a>(p: &'a Value, k: &str) -> Option<&'a str> {
    p.get(k).and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty())
}

impl App {
    /// The actions agents can run, with what they do and their keys.
    pub fn action_list() -> Value {
        let list: Vec<Value> = Action::all()
            .into_iter()
            .filter(|a| !needs_a_person(*a))
            .map(|a| json!({"id": a.id(), "name": a.name(), "keys": a.keys().map(super::actions::keys_label)}))
            .collect();
        json!({"actions": list})
    }

    /// Runs action `id` as its menu item, key or palette entry would.
    pub fn attach_action(&mut self, id: &str, ctx: &egui::Context) -> Result<Value, String> {
        let a = Action::by_id(id).ok_or_else(|| format!("no action `{id}`: `action` with `list` gives every name"))?;
        if needs_a_person(a) {
            return Err(format!("`{id}` opens a file dialog for the person; use the open, session, queue or bookmarks command with a path instead"));
        }
        self.run(a, ctx);
        // Comparing opens the comparison first; when it couldn't (fewer than two alignment
        // tracks), the error it set is the answer, not success.
        if matches!(a, Action::Compare | Action::NextDiff | Action::PrevDiff) && self.compare.is_none() {
            return Err(self.error.clone().unwrap_or_else(|| "no comparison could be opened".into()));
        }
        Ok(json!({"ran": id, "name": a.name()}))
    }

    /// The review now: which VCF, the site at the outlined base, the note being typed, counts.
    fn review_now(&self) -> Value {
        let Some(r) = &self.review else { return json!({"reviewing": false}) };
        let site = self.current_site().map(|(k, v)| {
            json!({
                "locus": format!("{}:{}", k.contig, k.pos), "ref": k.reference, "alt": k.alt, "id": v.id, "qual": v.qual,
                "decision": r.decisions.get(&k).map(|d| d.label.name()),
            })
        });
        let [real, artifact, unsure] = r.counts();
        json!({
            "reviewing": true,
            "vcf": r.source.display().to_string(),
            "file": r.path.display().to_string(),
            "site": site,
            "finding_next": site.is_none() && self.review_track().and_then(|i| self.variants.get(i)).is_some_and(|v| v.is_loading()),
            "note": r.note,
            "counts": {"real": real, "artifact": artifact, "unsure": unsure},
            "error": r.error, "notice": r.notice,
        })
    }

    /// `review` with an `action`: start (a VCF by path or file name, else the first), decide
    /// (real, artifact or unsure, with an optional note), note, next, prev or stop. Without
    /// one, every VCF's decisions.
    pub fn attach_review_command(&mut self, p: &Value) -> Result<Value, String> {
        let action = text(p, "action").unwrap_or("status");
        if action == "status" {
            let mut all = self.attach_review();
            all["current"] = self.review_now();
            return Ok(all);
        }
        match action {
            "start" => {
                if self.variants.is_empty() {
                    return Err("no VCF is open: open one first".into());
                }
                let i = match text(p, "vcf") {
                    None => 0,
                    Some(want) => {
                        let hits: Vec<usize> = self
                            .variants
                            .iter()
                            .enumerate()
                            .filter(|(_, v)| crate::files::names(&v.path, want))
                            .map(|(i, _)| i)
                            .collect();
                        match hits.as_slice() {
                            [i] => *i,
                            [] => return Err(format!("no open VCF is `{want}`")),
                            _ => return Err(format!("`{want}` names more than one open VCF: give the path")),
                        }
                    }
                };
                self.start_review(i);
            }
            "decide" => {
                let d = text(p, "decision").ok_or("give `decision`: real, artifact or unsure")?;
                let label = label_named(d).ok_or_else(|| format!("`{d}` isn't a decision: real, artifact or unsure"))?;
                if self.review.is_none() {
                    return Err("not reviewing: `review` with action `start` first".into());
                }
                let Some((key, _)) = self.current_site() else {
                    return Err("no site is outlined yet (the review is still finding it): call `review` again".into());
                };
                if let (Some(n), Some(r)) = (p.get("note").and_then(Value::as_str), &mut self.review) {
                    r.arrive(&key);
                    r.note = n.to_owned();
                }
                self.decide(label);
                if let Some(e) = self.review.as_ref().and_then(|r| r.error.clone()) {
                    return Err(format!("decided, but it couldn't be saved: {e}"));
                }
            }
            "note" => {
                let n = p.get("note").and_then(Value::as_str).ok_or("give the `note` (\"\" clears it)")?.to_owned();
                let key = self.current_site().map(|(k, _)| k);
                let r = self.review.as_mut().ok_or("not reviewing")?;
                if let Some(k) = key {
                    r.arrive(&k);
                }
                r.note = n;
            }
            "next" | "prev" => {
                self.review.as_ref().ok_or("not reviewing")?;
                self.review_track().ok_or("the reviewed VCF is no longer open")?;
                self.review_step(action == "next");
            }
            "stop" => self.review = None,
            other => return Err(format!("unknown review action `{other}`: status, start, decide, note, next, prev or stop")),
        }
        Ok(json!({"current": self.review_now()}))
    }

    fn queue_now(&self) -> Value {
        let Some(q) = &self.queue else { return json!({"open": false}) };
        let [pass, fail, unsure] = q.counts();
        let items: Vec<Value> = q
            .items
            .iter()
            .enumerate()
            .map(|(i, it)| {
                let d = q.decision(i);
                json!({"index": i + 1, "locus": it.key(), "alleles": it.alleles, "label": it.label, "decision": d.map(|d| d.verdict.name()), "note": d.map(|d| d.note.clone())})
            })
            .collect();
        json!({
            "open": true, "source": q.source, "position": q.pos + 1, "count": q.items.len(),
            "counts": {"pass": pass, "fail": fail, "unsure": unsure},
            "note": q.note, "file": q.save_path.as_ref().map(|p| p.display().to_string()),
            "error": q.error, "notice": q.notice, "items": items,
        })
    }

    /// `queue`: status, open (a list file at `path`), paste (`text`: one locus per line),
    /// decide (pass, fail or unsure, with an optional note; moves to the next undecided), note,
    /// next, prev, go (`index`, 1-based), next_undecided, stop, export (`path`, BED or TSV by
    /// its extension).
    pub fn attach_queue(&mut self, p: &Value) -> Result<Value, String> {
        let action = text(p, "action").unwrap_or("status");
        let need = |s: &Self| if s.queue.is_some() { Ok(()) } else { Err("no queue is open: `queue` with action `open` or `paste` first".to_owned()) };
        match action {
            "status" => {}
            "open" => {
                let path = std::path::PathBuf::from(text(p, "path").ok_or("give the list's `path` (BED, VCF or one locus per line)")?);
                if !path.is_file() {
                    return Err(format!("{}: no such file", path.display()));
                }
                // A list that fails (unreadable, no places) leaves any open queue as it was, so
                // the error is what says so, not whether a queue is open.
                self.error = None;
                self.open_queue_file(&path);
                if let Some(e) = self.error.clone() {
                    return Err(e);
                }
            }
            "paste" => {
                let t = p.get("text").and_then(Value::as_str).ok_or("give `text`: one locus per line (chr1:100-200, BED or VCF lines)")?;
                self.error = None;
                self.start_queue(crate::queue::Queue::for_contigs("pasted list".into(), t, None, &self.contigs));
                if let Some(e) = self.error.clone() {
                    return Err(e);
                }
            }
            "decide" => {
                need(self)?;
                let d = text(p, "decision").ok_or("give `decision`: pass, fail or unsure")?;
                let v = verdict_named(d).ok_or_else(|| format!("`{d}` isn't a decision: pass, fail or unsure"))?;
                if let (Some(n), Some(q)) = (p.get("note").and_then(Value::as_str), &mut self.queue) {
                    q.note = n.to_owned();
                }
                self.queue_decide(v);
                if let Some(e) = self.queue.as_ref().and_then(|q| q.error.clone()) {
                    return Err(format!("decided, but it couldn't be saved: {e}"));
                }
            }
            "note" => {
                need(self)?;
                let n = p.get("note").and_then(Value::as_str).ok_or("give the `note`")?.to_owned();
                if let Some(q) = &mut self.queue {
                    q.note = n;
                }
            }
            "next" | "prev" => {
                need(self)?;
                self.queue_step(action == "next");
            }
            "next_undecided" => {
                need(self)?;
                self.queue_next_undecided();
            }
            "go" => {
                need(self)?;
                let i = p.get("index").and_then(Value::as_u64).ok_or("give the item's `index` (1-based)")? as usize;
                let n = self.queue.as_ref().map_or(0, |q| q.items.len());
                if i == 0 || i > n {
                    return Err(format!("no item {i}: the queue has {n}"));
                }
                if let Some(q) = &mut self.queue {
                    q.go(i - 1);
                }
                self.show_queue_item();
            }
            "export" => {
                need(self)?;
                let path = std::path::PathBuf::from(text(p, "path").ok_or("give the `path` to write (.bed or .tsv)")?);
                super::attach::check_export(&path, &["bed", "tsv"], &self.protected_files())?;
                let q = self.queue.as_ref().expect("checked");
                let bed = path.extension().is_some_and(|e| e.eq_ignore_ascii_case("bed"));
                let body = if bed { q.to_bed() } else { q.to_tsv(true) };
                crate::store::write_atomic(&path, body.as_bytes()).map_err(|e| format!("{}: {e}", path.display()))?;
            }
            "stop" => self.queue = None,
            other => return Err(format!("unknown queue action `{other}`")),
        }
        Ok(self.queue_now())
    }

    /// The sequence search now.
    fn search_now(&self) -> Value {
        let Some(s) = &self.seq else { return json!({"open": false}) };
        let Some(x) = &s.search else { return json!({"open": true, "searching": false, "error": s.error}) };
        let hits: Vec<Value> = x
            .hits
            .iter()
            .take(500)
            .enumerate()
            .map(|(i, h)| json!({"index": i + 1, "locus": format!("{}:{}-{}", h.contig, h.start + 1, h.start + x.len), "mismatches": h.mismatches, "strand": if h.reverse { "-" } else { "+" }}))
            .collect();
        json!({
            "open": true, "pattern": x.pattern, "mismatches": x.mismatches, "done": x.done, "progress": x.progress(),
            "hits_total": x.hits.len(), "truncated": x.truncated, "hits": hits, "error": x.error.clone().or(s.error.clone()),
        })
    }

    /// `search`: start one for `sequence` (with up to `mismatches`, in `scope` genome, contig or
    /// view), status, go to hit `index`, or stop (cancels a search in progress and closes it).
    pub fn attach_search(&mut self, p: &Value, ctx: &egui::Context) -> Result<Value, String> {
        use super::tools::{Scope, SeqUi};
        match text(p, "action").unwrap_or(if p.get("sequence").is_some() { "start" } else { "status" }) {
            "status" => {}
            "start" => {
                let seq = text(p, "sequence").ok_or("give the `sequence` (A, C, G, T, IUPAC codes)")?;
                if self.reference_path.is_none() {
                    return Err("no reference is open: searching needs one".into());
                }
                let scope = match text(p, "scope").unwrap_or("genome") {
                    "genome" => Scope::Genome,
                    "contig" => Scope::Contig,
                    "view" => Scope::View,
                    other => return Err(format!("`{other}` isn't a scope: genome, contig or view")),
                };
                let mismatches = p.get("mismatches").and_then(Value::as_u64).unwrap_or(0).min(3) as u8;
                if let Some(old) = self.seq.as_ref().and_then(|s| s.search.as_ref()) {
                    old.cancel();
                }
                let mut s = SeqUi { text: seq.to_owned(), mismatches, scope, ..Default::default() };
                self.start_seq_search(&mut s, ctx);
                if let Some(e) = &s.error {
                    return Err(e.clone());
                }
                self.seq = Some(s);
            }
            "go" => {
                let i = p.get("index").and_then(Value::as_u64).ok_or("give the hit's `index` (1-based)")? as usize;
                let (hit, len) = {
                    let x = self.seq.as_ref().and_then(|s| s.search.as_ref()).ok_or("no search: start one first")?;
                    (x.hits.get(i.wrapping_sub(1)).cloned().ok_or_else(|| format!("no hit {i}: {} so far", x.hits.len()))?, x.len)
                };
                self.jump_to(&hit.contig, hit.start + len / 2);
                self.zoom_to_show(len);
                if let Some(c) = self.contigs.iter().position(|x| x.name == hit.contig) {
                    self.selection = Some((c, hit.start, hit.start + len));
                }
            }
            "stop" => {
                if let Some(x) = self.seq.as_ref().and_then(|s| s.search.as_ref()) {
                    x.cancel();
                }
                self.seq = None;
            }
            other => return Err(format!("unknown search action `{other}`: start, status, go or stop")),
        }
        Ok(self.search_now())
    }

    /// The alignment track `want` names: its path, or a file name only one track has (two
    /// `same.bam` in different folders need the path, as for `tracks`; not the first of them).
    fn track_named(&self, want: &str) -> Result<usize, String> {
        let paths: Vec<std::path::PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        let path = super::attach::unique_file(&paths, want)?;
        Ok(paths.iter().position(|p| *p == path).unwrap_or_default())
    }

    /// `compare`: open the visible comparison of alignment tracks `a` and `b` (1-based, or path
    /// or file name) with a `threshold` and `min_depth`, step to the `next` or `prev`
    /// difference, or `close` it.
    pub fn attach_compare(&mut self, p: &Value, ctx: &egui::Context) -> Result<Value, String> {
        let track = |s: &Self, k: &str| -> Result<Option<usize>, String> {
            let Some(v) = p.get(k) else { return Ok(None) };
            if let Some(n) = v.as_u64() {
                return if n >= 1 && (n as usize) <= s.tracks.len() { Ok(Some(n as usize - 1)) } else { Err(format!("no alignment track {n}: {} are open", s.tracks.len())) };
            }
            s.track_named(v.as_str().unwrap_or_default()).map(Some)
        };
        match text(p, "action").unwrap_or("open") {
            "open" => {
                if self.tracks.len() < 2 {
                    return Err("comparing needs two alignment tracks open".into());
                }
                // Everything checked before the visible comparison changes: a refused call
                // leaves it as it was.
                let (a, b) = (track(self, "a")?, track(self, "b")?);
                let (cur_a, cur_b) = self.compare.as_ref().map_or((0, 1), |c| (c.a, c.b));
                let (a, b) = (a.unwrap_or(cur_a), b.unwrap_or(cur_b));
                if a == b {
                    return Err("`a` and `b` are the same track".into());
                }
                // The window's ranges: a difference must be more than `threshold` (at 1, a
                // full 100% difference), with at least `min_depth` reads in both.
                let (ts, ds) = (crate::diff::THRESHOLDS, crate::diff::MIN_DEPTHS);
                let threshold = match p.get("threshold") {
                    None | Some(Value::Null) => None,
                    Some(v) => Some(v.as_f64().filter(|t| ts.contains(t)).ok_or_else(|| format!("`threshold` is an allele-fraction difference from {} to {}", ts.start(), ts.end()))?),
                };
                let min_depth = match p.get("min_depth") {
                    None | Some(Value::Null) => None,
                    Some(v) => Some(
                        v.as_u64().and_then(|d| u32::try_from(d).ok()).filter(|d| ds.contains(d)).ok_or_else(|| format!("`min_depth` is {} to {} reads", ds.start(), ds.end()))?,
                    ),
                };
                self.open_compare();
                self.set_compare_tracks(a, b);
                let c = self.compare.as_mut().expect("opened");
                if let Some(t) = threshold {
                    c.settings.threshold = t;
                }
                if let Some(d) = min_depth {
                    c.settings.min_depth = d;
                }
            }
            "next" | "prev" => {
                self.sync_compare();
                if self.compare.is_none() {
                    return Err("no comparison open: `compare` with action `open` first".into());
                }
                let forward = text(p, "action") == Some("next");
                self.step_diff(forward, ctx);
            }
            "close" => self.compare = None,
            other => return Err(format!("unknown compare action `{other}`: open, next, prev or close")),
        }
        Ok(match &self.compare {
            None => json!({"open": false}),
            Some(c) => json!({
                "open": true,
                "a": self.tracks.get(c.a).map(|t| t.loader.bam_path.display().to_string()),
                "b": self.tracks.get(c.b).map(|t| t.loader.bam_path.display().to_string()),
                "threshold": c.settings.threshold, "min_depth": c.settings.min_depth,
                "scanning": c.scan.is_some(), "status": c.status,
            }),
        })
    }
}

impl App {
    /// Each alignment track's read rows: how many, which are on screen, and its groups.
    pub(super) fn rows_now(&self) -> Value {
        let tracks: Vec<Value> = self
            .tracks
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let style = self.track_style(&t.loader.bam_path);
                let Some((rl, cur)) = t.row_layout(style).zip(t.current.as_ref()) else {
                    return json!({"track": i + 1, "file": t.loader.bam_path.display().to_string(), "rows": null, "note": "no reads loaded here (zoom in)"});
                };
                let shown_h = (t.view_h - crate::tracks::READS_TOP).max(0.0);
                let first = (t.scroll_y / rl.pitch).floor() as usize;
                let last = (((t.scroll_y + shown_h) / rl.pitch).ceil() as usize).min(rl.rows);
                let groups: Vec<Value> = cur.packed.groups.iter().map(|g| json!({"label": g.label, "first_row": g.first_row + 1})).collect();
                json!({
                    "track": i + 1, "file": t.loader.bam_path.display().to_string(), "rows": rl.rows,
                    "first_shown": (rl.rows > 0).then_some(first + 1), "last_shown": (rl.rows > 0).then_some(last),
                    "groups": groups,
                })
            })
            .collect();
        json!({"tracks": tracks})
    }

    /// `rows`: scrolls an alignment track's read rows (the person drags vertically, or uses
    /// the arrow and page keys): to `row` (1-based, shown at the top), `by` some rows, `to` the
    /// `top` or `bottom`, to the row holding read `read` (which is also outlined), or to the
    /// first row of group `group`. `track` is 1-based, a path or a file name (default: every
    /// track for `to`, else the first). Without any of these, where each track is.
    pub fn attach_rows(&mut self, p: &Value) -> Result<Value, String> {
        let which: Vec<usize> = match p.get("track") {
            None => (0..self.tracks.len()).collect(),
            Some(v) => {
                let i = if let Some(n) = v.as_u64() {
                    (n >= 1 && (n as usize) <= self.tracks.len()).then(|| n as usize - 1).ok_or_else(|| format!("no alignment track {v}: {} open", self.tracks.len()))?
                } else {
                    self.track_named(v.as_str().unwrap_or_default())?
                };
                vec![i]
            }
        };
        if self.tracks.is_empty() {
            return Err("no alignment track is open".into());
        }
        let moves = ["row", "by", "to", "read", "group"].iter().filter(|k| p.get(**k).is_some()).count();
        if moves > 1 {
            return Err("give one of `row`, `by`, `to`, `read` or `group`".into());
        }
        // `row`/`by`/`group` act on one track: the named one, or the first.
        let one = if p.get("track").is_some() { which[0] } else { 0 };
        let pitch_of = |s: &Self, i: usize| s.tracks[i].row_layout(s.track_style(&s.tracks[i].loader.bam_path)).map(|rl| (rl.pitch, rl.rows));
        if let Some(r) = p.get("row").and_then(Value::as_u64) {
            let (pitch, rows) = pitch_of(self, one).ok_or("that track has no reads loaded here (zoom in)")?;
            if r == 0 || r as usize > rows.max(1) {
                return Err(format!("no row {r}: the track has {rows}"));
            }
            self.tracks[one].scroll_y = (r - 1) as f32 * pitch;
        } else if let Some(b) = p.get("by").and_then(Value::as_i64) {
            let (pitch, _) = pitch_of(self, one).ok_or("that track has no reads loaded here (zoom in)")?;
            let t = &mut self.tracks[one];
            t.scroll_y = (t.scroll_y + b as f32 * pitch).max(0.0);
        } else if let Some(to) = text(p, "to") {
            for &i in &which {
                self.tracks[i].scroll_y = match to {
                    "top" => 0.0,
                    "bottom" => f32::MAX / 4.0, // clamped to the last row when drawn
                    other => return Err(format!("`to` is top or bottom, not `{other}`")),
                };
            }
        } else if let Some(name) = text(p, "read") {
            // The first track (of those named) with the read, scrolled so its row is near the top.
            let mut found = None;
            for &i in &which {
                let t = &self.tracks[i];
                let Some(cur) = &t.current else { continue };
                let Some((reads, _)) = &cur.view.reads else { continue };
                // The rows as drawn (grouped and sorted), not the load's own packing.
                if let Some(row) = cur.packed.rows.iter().position(|row| row.iter().any(|&k| reads.get(k).is_some_and(|r| r.name == name))) {
                    found = Some((i, row));
                    break;
                }
            }
            let (i, row) = found.ok_or_else(|| format!("no read named `{name}` is loaded in view (or it is filtered out)"))?;
            let (pitch, _) = pitch_of(self, i).ok_or("no reads loaded")?;
            self.tracks[i].scroll_y = (row.saturating_sub(2)) as f32 * pitch;
            self.selected_read = Some(name.to_owned());
        } else if let Some(label) = text(p, "group") {
            let first = self.tracks[one].current.as_ref().and_then(|c| c.packed.groups.iter().find(|g| g.label == label).map(|g| g.first_row));
            let first = first.ok_or_else(|| format!("no group `{label}` in that track (`rows` lists them)"))?;
            let (pitch, _) = pitch_of(self, one).ok_or("no reads loaded")?;
            self.tracks[one].scroll_y = first as f32 * pitch;
        }
        // Clamp now, so the reply says where the rows really are.
        for i in 0..self.tracks.len() {
            let h = self.tracks[i].view_h;
            if h > 0.0 {
                let max = self.tracks[i].max_scroll(self.track_style(&self.tracks[i].loader.bam_path), h);
                self.tracks[i].scroll_y = self.tracks[i].scroll_y.clamp(0.0, max);
            }
        }
        Ok(self.rows_now())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::Hit;

    #[test]
    fn every_action_has_a_unique_id_that_finds_it() {
        let all = Action::all();
        let ids: Vec<String> = all.iter().map(|a| a.id()).collect();
        for (a, id) in all.iter().zip(&ids) {
            assert_eq!(Action::by_id(id), Some(*a), "{id}");
            assert_eq!(ids.iter().filter(|x| *x == id).count(), 1, "{id} is used twice");
        }
        assert_eq!(Action::NextVariant.id(), "next_variant");
        assert_eq!(Action::Scale(CoverageScale::Log).id(), "scale:log");
        assert!(ids.contains(&"color:strand".to_owned()));
    }

    #[test]
    fn review_follows_its_vcf_when_tracks_move_or_close() {
        let data = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let dir = std::env::temp_dir().join(format!("ggv-review-move-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut vcfs = Vec::new();
        for name in ["HG002.revio.deepvariant.phased.vcf.gz", "HG002.revio.pbsv.vcf.gz"] {
            for ext in ["", ".tbi"] {
                std::fs::copy(data.join(format!("{name}{ext}")), dir.join(format!("{name}{ext}"))).unwrap();
            }
            vcfs.push(dir.join(name));
        }
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { vcfs: vcfs.clone(), ..Default::default() }).unwrap();
        app.start_review(0);
        let sidecar = app.review.as_ref().unwrap().path.clone();
        assert!(sidecar.to_string_lossy().contains("deepvariant"));

        // Moved down: the review stays on the first VCF, now second, and its sidecar.
        app.move_track(super::super::Hit::Variants(0), false);
        assert_eq!(app.variants[1].path, vcfs[0]);
        assert_eq!(app.review_track(), Some(1));
        assert_eq!(app.review_now()["vcf"], vcfs[0].display().to_string());
        assert_eq!(app.review.as_ref().unwrap().path, sidecar);
        let listed = app.attach_review();
        assert_eq!(listed["review"][1]["reviewing"], true);
        assert_eq!(listed["review"][0]["reviewing"], false);

        // Closing the other VCF keeps it; closing the reviewed one ends the review.
        app.close_variants(0);
        assert_eq!(app.review_track(), Some(0));
        app.close_variants(0);
        assert!(app.review.is_none());
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A pasted BED row on a contig with `:` in its name (a `samtools faidx` sub-reference) goes
    /// to that row in a window with that contig, also when the genome arrives after the queue,
    /// and exports as the row.
    #[test]
    fn queue_rows_on_contigs_with_colons_go_to_the_row() {
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        let text = "chr20:10000000-20000000\t500\t501\tsiteA\n";
        let sub = vec![ggv_model::Contig { name: "chr20:10000000-20000000".into(), length: 10_000_001 }];
        app.contigs = sub.clone();
        app.pending_goto = None;
        app.attach_queue(&json!({"action": "paste", "text": text})).unwrap();
        assert_eq!(app.pending_goto.as_deref(), Some("chr20:10000000-20000000:501"));
        assert_eq!(app.queue_now()["items"][0]["locus"], "chr20:10000000-20000000:501");
        app.queue = None;
        // The genome arrives after the queue was opened.
        app.contigs.clear();
        app.attach_queue(&json!({"action": "paste", "text": text})).unwrap();
        assert_eq!(app.queue_now()["items"][0]["label"], "500 501 siteA", "no genome: the primary reading");
        app.contigs = sub;
        app.refresh_decisions();
        assert_eq!(app.queue_now()["items"][0]["label"], "siteA");
        app.queue_step(true);
        assert_eq!(app.pending_goto.as_deref(), Some("chr20:10000000-20000000:501"));
        assert_eq!(app.queue.as_mut().map(|q| { q.save_path = None; q.decide(crate::queue::Verdict::Pass); q.to_bed() }).unwrap(), "chr20:10000000-20000000\t500\t501\tpass\n");
    }

    #[test]
    fn exports_refuse_open_files_and_other_formats() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-export-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        let list = dir.join("sites.bed");
        std::fs::write(&list, "chr1\t0\t10\ta\n").unwrap();
        app.attach_queue(&json!({"action": "open", "path": list.display().to_string()})).unwrap();
        app.attach_queue(&json!({"action": "decide", "decision": "pass"})).unwrap();
        let json_file = dir.join("existing.json");
        std::fs::write(&json_file, "{\"keep\":\"original\"}").unwrap();
        let err = app.attach_queue(&json!({"action": "export", "path": json_file.display().to_string()})).unwrap_err();
        assert!(err.contains(".bed or .tsv"), "{err}");
        assert_eq!(std::fs::read_to_string(&json_file).unwrap(), "{\"keep\":\"original\"}");
        // The queue's own list is open too.
        assert!(app.attach_queue(&json!({"action": "export", "path": list.display().to_string()})).is_err());
        assert_eq!(std::fs::read_to_string(&list).unwrap(), "chr1\t0\t10\ta\n");
        app.attach_queue(&json!({"action": "export", "path": dir.join("out.bed").display().to_string()})).unwrap();
        assert_eq!(std::fs::read_to_string(dir.join("out.bed")).unwrap(), "chr1\t0\t10\tpass\n");
        // Bookmarks: only .bed, and never over an open file.
        let err = app.attach_bookmarks(&json!({"action": "export", "path": list.display().to_string()})).unwrap_err();
        assert!(err.contains("open in GGV"), "{err}");
        assert!(app.attach_bookmarks(&json!({"action": "export", "path": json_file.display().to_string()})).is_err());
        let _ = std::fs::remove_dir_all(dir);
    }

    /// One frame of the Missing indexes window, with `events`; the AccessKit tree it drew.
    fn index_frame(app: &mut App, ctx: &egui::Context, events: Vec<egui::Event>) -> Option<egui::accesskit::TreeUpdate> {
        let raw = egui::RawInput { events, screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1200.0, 800.0))), ..Default::default() };
        let mut out = ctx.run_pass(raw, |ctx| app.index_window(ctx));
        out.textures_delta.clear();
        out.platform_output.accesskit_update
    }

    /// Clicks the button labelled `label` in the Missing indexes window, as a person would.
    fn click(app: &mut App, ctx: &egui::Context, label: &str) {
        // A new window is laid out unseen on its first frame.
        index_frame(app, ctx, vec![]);
        let tree = index_frame(app, ctx, vec![]).expect("accesskit enabled");
        let b = tree
            .nodes
            .iter()
            .find(|(_, n)| n.role() == egui::accesskit::Role::Button && n.label() == Some(label))
            .and_then(|(_, n)| n.bounds())
            .unwrap_or_else(|| panic!("no `{label}` button in {:?}", tree.nodes.iter().map(|(_, n)| (n.role(), n.label())).collect::<Vec<_>>()));
        let at = egui::pos2(((b.x0 + b.x1) / 2.0) as f32, ((b.y0 + b.y1) / 2.0) as f32);
        let button = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
        index_frame(app, ctx, vec![egui::Event::PointerMoved(at), button(true)]);
        index_frame(app, ctx, vec![button(false)]);
    }

    /// Runs frames until no index is being built (or 60 s pass).
    fn until_built(app: &mut App, ctx: &egui::Context) {
        let t = std::time::Instant::now();
        while app.indexing.iter().any(|j| matches!(j.state, super::super::IndexState::Building(_))) && t.elapsed().as_secs() < 60 {
            std::thread::sleep(std::time::Duration::from_millis(20));
            index_frame(app, ctx, vec![]);
        }
        index_frame(app, ctx, vec![]);
    }

    /// Issue #11: alignments opened with a reference that has no `.fai` open once the window's
    /// Build .fai button has indexed it, offering their own indexes next when they lack them.
    #[test]
    fn building_the_reference_index_opens_the_alignments_asked_for() {
        let data = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        if !data.join("HG002.revio.slices.bam").exists() {
            return;
        }
        // (BAMs, whether they keep their .bai)
        for (case, (names, indexed)) in [(vec!["reads.bam"], true), (vec!["reads.bam"], false), (vec!["a.bam", "b.bam"], true)].into_iter().enumerate() {
            let dir = std::env::temp_dir().join(format!("ggv-fai-click-{}-{case}", std::process::id()));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            let fa = dir.join("ref.fa.gz");
            std::fs::copy(data.join("GRCh38.chr20.slices.fa.gz"), &fa).unwrap();
            // Case 1 keeps samtools' .gzi; the others build it too.
            if case == 1 {
                std::fs::copy(data.join("GRCh38.chr20.slices.fa.gz.gzi"), dir.join("ref.fa.gz.gzi")).unwrap();
            }
            let mut paths = vec![];
            for n in &names {
                std::fs::copy(data.join("HG002.revio.slices.bam"), dir.join(n)).unwrap();
                if indexed {
                    std::fs::copy(data.join("HG002.revio.slices.bam.bai"), dir.join(format!("{n}.bai"))).unwrap();
                }
                paths.push(dir.join(n));
            }
            paths.push(fa.clone());

            let ctx = egui::Context::default();
            ctx.enable_accesskit();
            let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
            app.open_paths(paths, &ctx);
            assert!(app.tracks.is_empty(), "case {case}: opened before the reference was indexed");
            click(&mut app, &ctx, "Build .fai");
            until_built(&mut app, &ctx);
            assert!(dir.join("ref.fa.gz.fai").exists() && dir.join("ref.fa.gz.gzi").exists(), "case {case}");
            assert_eq!(app.reference_path.as_ref(), Some(&fa), "case {case}");
            if !indexed {
                // The alignment's own index is offered next; building it opens the file.
                assert!(app.tracks.is_empty(), "case {case}");
                click(&mut app, &ctx, "Build .bai");
                until_built(&mut app, &ctx);
            }
            let open: Vec<_> = app.tracks.iter().map(|t| t.loader.bam_path.file_name().unwrap().to_string_lossy().into_owned()).collect();
            assert_eq!(open, names, "case {case}");
            assert!(app.indexing.is_empty() && app.deferred_bams.is_empty(), "case {case}");
            drop(app);
            let _ = std::fs::remove_dir_all(dir);
        }

        // A reference without an index opened while alignments are open on another: they stay
        // open on the one they have; dismissing the offer keeps both, building it switches.
        let dir = std::env::temp_dir().join(format!("ggv-fai-switch-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("new")).unwrap();
        let fa = dir.join("ref.fa.gz");
        for ext in ["", ".fai", ".gzi"] {
            std::fs::copy(data.join(format!("GRCh38.chr20.slices.fa.gz{ext}")), dir.join(format!("ref.fa.gz{ext}"))).unwrap();
        }
        let noidx = dir.join("new").join("noidx.fa.gz");
        std::fs::copy(data.join("GRCh38.chr20.slices.fa.gz"), &noidx).unwrap();
        let bams = [dir.join("a.bam"), dir.join("b.bam")];
        for b in &bams {
            std::fs::copy(data.join("HG002.revio.slices.bam"), b).unwrap();
            std::fs::copy(data.join("HG002.revio.slices.bam.bai"), b.with_added_extension("bai")).unwrap();
        }
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        app.open_paths(vec![bams[0].clone(), fa.clone()], &ctx);
        assert_eq!((app.tracks.len(), app.reference_path.as_ref()), (1, Some(&fa)));
        app.open_paths(vec![noidx.clone(), bams[1].clone()], &ctx);
        assert_eq!((app.tracks.len(), app.reference_path.as_ref()), (1, Some(&fa)), "the open alignments stay until it is indexed");
        assert_eq!(app.deferred_bams, vec![bams[1].clone()]);
        click(&mut app, &ctx, "Dismiss");
        assert_eq!(app.reference_path.as_ref(), Some(&fa), "dismissed: the reference they had");
        let open: Vec<_> = app.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        assert_eq!(open, bams.to_vec(), "nothing lost: the one asked for with it opens on the old reference");
        assert!(app.indexing.is_empty() && app.deferred_bams.is_empty() && app.pending_reference.is_none());
        app.open_paths(vec![noidx.clone()], &ctx);
        click(&mut app, &ctx, "Build .fai");
        until_built(&mut app, &ctx);
        assert_eq!((app.tracks.len(), app.reference_path.as_ref()), (2, Some(&noidx)), "built: switched, the alignments reopened on it");
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A BAM without its BGZF end-of-file marker opens, with a warning that doesn't count as an
    /// error (so it doesn't fail an automated run). Each truncated file's warning shows while it
    /// is open, and goes with it.
    #[test]
    fn a_truncated_bam_opens_with_a_warning() {
        let data = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        if !data.join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-no-eof-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bytes = std::fs::read(data.join("HG002.revio.slices.bam")).unwrap();
        let cut: Vec<_> = ["cut.bam", "cut2.bam"].iter().map(|n| dir.join(n)).collect();
        for bam in &cut {
            std::fs::write(bam, &bytes[..bytes.len() - 28]).unwrap();
            std::fs::copy(data.join("HG002.revio.slices.bam.bai"), bam.with_added_extension("bai")).unwrap();
        }
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        app.open_paths(vec![cut[0].clone()], &ctx);
        assert_eq!(app.tracks.len(), 1);
        assert!(app.error.is_none(), "{:?}", app.error);
        assert!(app.warning().is_some_and(|w| w.contains("end-of-file")), "{:?}", app.warning());
        app.open_paths(vec![cut[1].clone()], &ctx);
        let w = app.warning().unwrap();
        assert!(w.contains("cut.bam") && w.contains("cut2.bam"), "both files' warnings: {w}");
        app.close_track(0);
        let w = app.warning().unwrap();
        assert!(!w.contains("cut.bam:") && w.contains("cut2.bam"), "the closed file's warning goes with it: {w}");
        assert_eq!(app.attach_state()["warning"], json!(w));
        app.apply_session(crate::session::Session::default(), false, &ctx);
        assert_eq!(app.warning(), None, "a workspace without them has no warning");
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    fn data_dir() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio")
    }

    /// Copies test files (with their indexes) into `dir`.
    fn copy_data(dir: &std::path::Path, names: &[&str]) -> Vec<std::path::PathBuf> {
        std::fs::create_dir_all(dir).unwrap();
        for entry in std::fs::read_dir(data_dir()).unwrap().flatten() {
            let n = entry.file_name().to_string_lossy().into_owned();
            if names.iter().any(|w| n.starts_with(w)) {
                std::fs::copy(entry.path(), dir.join(&n)).unwrap();
            }
        }
        names.iter().map(|n| dir.join(n)).collect()
    }

    /// Loads (without drawing) until `done`, or fails after 60 s.
    fn settle(app: &mut App, done: impl Fn(&App) -> bool) {
        let t = std::time::Instant::now();
        while !done(app) && t.elapsed().as_secs() < 60 {
            app.load_hidden(1.0);
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(done(app), "timed out loading");
    }

    /// The review's site once it has one: (1-based position, REF length, ALT).
    fn site(app: &mut App) -> (u64, usize, String) {
        settle(app, |a| a.variants.iter().all(|v| !v.is_seeking()) && a.current_site().is_some());
        let (k, _) = app.current_site().unwrap();
        (k.pos, k.reference.len(), k.alt)
    }

    /// Two records start at chr20:31185371 (an insertion and a deletion): the review visits and
    /// decides each, stepping forward and back.
    #[test]
    fn review_steps_through_every_record_at_one_position() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-review-same-pos-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir, &["HG002.revio.slices.bam", "HG002.revio.pbsv.vcf.gz"]);
        let ctx = egui::Context::default();
        let args = super::super::Args { bams: vec![files[0].clone()], vcfs: vec![files[1].clone()], ..Default::default() };
        let mut app = App::new(&ctx, args).unwrap();
        app.goto("chr20:31185371");
        settle(&mut app, |a| a.variants.iter().all(|v| !v.is_loading() && v.current.is_some()));
        app.start_review(0);
        assert_eq!(site(&mut app).0, 31_185_371);
        assert_eq!(site(&mut app).1, 1, "the insertion first, in file order");
        app.attach_review_command(&json!({"action": "decide", "decision": "real"})).unwrap();
        assert_eq!(site(&mut app), (31_185_371, 51, "A".to_owned()), "then the deletion at the same base");
        app.attach_review_command(&json!({"action": "decide", "decision": "artifact"})).unwrap();
        assert_eq!(site(&mut app).0, 31_185_843, "then the next base with a record");
        // Back: the last record at the previous base, then the one before it there.
        app.attach_review_command(&json!({"action": "prev"})).unwrap();
        assert_eq!(site(&mut app), (31_185_371, 51, "A".to_owned()));
        app.attach_review_command(&json!({"action": "prev"})).unwrap();
        assert_eq!(site(&mut app).1, 1);
        assert_eq!(app.review_now()["site"]["decision"], "real");
        app.attach_review_command(&json!({"action": "prev"})).unwrap();
        assert_eq!(site(&mut app).0, 31_178_448);
        assert_eq!(app.review.as_ref().unwrap().counts(), [1, 1, 0], "both records at one base decided");
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A review carries on when its VCF is reopened (a session, a reload), and ends cleanly when
    /// a session doesn't have it.
    #[test]
    fn review_survives_reopening_its_vcf() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-review-reopen-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir.join("data"), &["HG002.revio.slices.bam", "HG002.revio.pbsv.vcf.gz", "HG002.revio.deepvariant.phased.vcf.gz"]);
        let ctx = egui::Context::default();
        let args = super::super::Args { bams: vec![files[0].clone()], vcfs: vec![files[1].clone(), files[2].clone()], ..Default::default() };
        let mut app = App::new(&ctx, args).unwrap();
        app.goto("chr20:31185371");
        settle(&mut app, |a| a.variants.iter().all(|v| !v.is_loading() && v.current.is_some()));
        app.attach_review_command(&json!({"action": "start", "vcf": "HG002.revio.pbsv.vcf.gz"})).unwrap();
        site(&mut app);
        app.attach_review_command(&json!({"action": "decide", "decision": "real"})).unwrap();
        // The VCF changed on disk (auto-reload) or Reload all.
        app.reload_path(&files[1], &ctx);
        assert_eq!(app.review_track(), Some(0), "bound to the reopened track");
        assert_eq!(app.attach_review()["review"][0]["reviewing"], true);
        assert_eq!(site(&mut app).1, 51, "still on the second record at that base");
        app.attach_review_command(&json!({"action": "decide", "decision": "unsure"})).unwrap();
        // A session with the same files.
        let session = dir.join("s1.ggv.json");
        app.attach_save_session(&json!({"path": session.display().to_string()})).unwrap();
        app.attach_open(vec![session.clone()], false, &ctx).unwrap();
        assert_eq!(app.review_track(), Some(0));
        app.attach_review_command(&json!({"action": "next"})).expect("the reviewed VCF is open");
        site(&mut app);
        app.attach_review_command(&json!({"action": "decide", "decision": "artifact"})).unwrap();
        assert_eq!(app.review.as_ref().unwrap().counts(), [1, 1, 1], "every decision recorded, in one file");
        // A workspace without it: the review ends, saying so.
        app.apply_session(crate::session::Session { bams: vec![files[0].clone()], ..Default::default() }, false, &ctx);
        assert!(app.review.is_none());
        assert!(app.notice.as_ref().is_some_and(|(n, _)| n.contains("review ended")), "{:?}", app.notice);
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// After a session is opened (relative paths joined to its folder), window commands still
    /// find files by any spelling of their path.
    #[test]
    fn files_are_found_by_any_spelling_after_a_session() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-session-spelling-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir.join("data"), &["HG002.revio.slices.bam", "HG002.revio.pbsv.vcf.gz", "HG002.revio.deepvariant.phased.vcf.gz"]);
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        let ctx = egui::Context::default();
        let args = super::super::Args { bams: vec![files[0].clone()], vcfs: vec![files[1].clone(), files[2].clone()], ..Default::default() };
        let mut app = App::new(&ctx, args).unwrap();
        let session = dir.join("s1.ggv.json");
        app.attach_save_session(&json!({"path": session.display().to_string()})).unwrap();
        app.attach_open(vec![session], false, &ctx).unwrap();
        #[cfg(windows)]
        for v in &app.variants {
            assert!(!v.path.to_string_lossy().contains('/'), "one separator: {}", v.path.display());
        }
        let slashes = |p: &std::path::Path| p.to_string_lossy().replace('\\', "/");
        let dotted = |p: &std::path::Path| dir.join("sub").join("..").join("data").join(p.file_name().unwrap()).display().to_string();
        for spell in [files[1].display().to_string(), slashes(&files[1]), dotted(&files[1])] {
            app.attach_review_command(&json!({"action": "start", "vcf": spell})).unwrap();
            assert_eq!(app.review_track(), Some(0), "{spell}");
            app.attach_display(&json!({"variant_color": "af", "variant_file": spell})).unwrap();
            app.attach_tracks(&json!({"action": "move", "track": spell, "direction": "up"})).unwrap();
        }
        for spell in [slashes(&files[0]), dotted(&files[0])] {
            app.attach_rows(&json!({"track": spell})).unwrap();
            app.attach_display(&json!({"file": spell, "expression": "mapq > 1"})).unwrap();
        }
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Exports never overwrite the saved review decisions (or their backup) of an open VCF,
    /// reviewed now or not.
    #[test]
    fn exports_refuse_every_open_vcfs_review_file() {
        let dir = std::env::temp_dir().join(format!("ggv-protect-review-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir, &["HG002.revio.pbsv.vcf.gz"]);
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { vcfs: files.clone(), ..Default::default() }).unwrap();
        let saved = "#contig\tpos\tref\talt\tdecision\tnote\nchr20\t31185371\tA\tAG\treal\t\n";
        let sidecar = crate::store::sidecar(&files[0], crate::review::EXT);
        let bak = crate::store::backup_of(&sidecar);
        std::fs::write(&sidecar, saved).unwrap();
        std::fs::write(&bak, saved).unwrap();
        app.attach_queue(&json!({"action": "paste", "text": "chr20:31180455\n"})).unwrap();
        // Not running a review now: before, only the running review's file was refused, and
        // writing over it moved the decisions into the backup, the next write losing them.
        assert!(app.review.is_none());
        let err = app.attach_queue(&json!({"action": "export", "path": sidecar.display().to_string()})).unwrap_err();
        assert!(err.contains("open in GGV"), "{err}");
        assert!(app.protected_files().contains(&bak));
        for p in [&sidecar, &bak] {
            assert_eq!(std::fs::read_to_string(p).unwrap(), saved);
        }
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// With no files open there is no place: adding a bookmark or overlay there is refused (it
    /// was saved at an empty contig, `:1`), as Ctrl+D does nothing.
    #[test]
    fn bookmarks_and_overlays_need_a_place() {
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        let err = app.attach_bookmarks(&json!({"action": "add", "note": "x"})).unwrap_err();
        assert!(err.contains("no files are open"), "{err}");
        assert!(app.attach_overlays(&json!({"action": "add", "label": "y"})).unwrap_err().contains("no files are open"));
        assert!(app.bookmarks.borrow().is_empty() && app.overlays.borrow().is_empty());
    }

    /// Reloading a VCF keeps how its sites are coloured.
    #[test]
    fn reload_keeps_the_site_colouring() {
        let dir = std::env::temp_dir().join(format!("ggv-reload-color-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir, &["HG002.revio.pbsv.vcf.gz"]);
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { vcfs: files.clone(), ..Default::default() }).unwrap();
        app.attach_display(&json!({"variant_color": "SVLEN", "variant_file": "HG002.revio.pbsv.vcf.gz"})).unwrap();
        let id = app.variants[0].id;
        app.reload_path(&files[0], &ctx);
        assert_ne!(app.variants[0].id, id, "reopened");
        assert_eq!(app.variants[0].color, crate::variants::SiteColor::Info("SVLEN".into()));
        assert_eq!(app.attach_state()["display"]["variant_colors"]["HG002.revio.pbsv.vcf.gz"], "SVLEN");
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Compare: 100% finds a full difference, the status says what a step found, settings are
    /// refused outside the window's ranges, and actions that can't compare say so.
    #[test]
    fn compare_at_100_percent_and_its_status() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-compare-full-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let files = copy_data(&dir, &["HG002.revio.slices.bam", "GRCh38.chr20.slices.fa.gz"]);
        let b = dir.join("b.bam");
        std::fs::copy(&files[0], &b).unwrap();
        std::fs::copy(files[0].with_added_extension("bai"), b.with_added_extension("bai")).unwrap();
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        for a in ["compare", "next_diff", "prev_diff"] {
            let err = app.attach_action(a, &ctx).unwrap_err();
            assert!(err.contains("two alignment tracks"), "{a}: {err}");
        }
        app.open_paths(vec![files[0].clone(), b.clone(), files[1].clone()], &ctx);
        // One haplotype in each: at a heterozygous site one is all alt, the other all reference.
        app.attach_display(&json!({"file": files[0].display().to_string(), "expression": "hp == 1"})).unwrap();
        app.attach_display(&json!({"file": b.display().to_string(), "expression": "hp == 2"})).unwrap();
        for bad in [json!({"threshold": 0}), json!({"threshold": 1.5}), json!({"min_depth": 0}), json!({"min_depth": 1001})] {
            assert!(app.attach_compare(&bad, &ctx).is_err(), "{bad}");
        }
        assert!(app.compare.is_none(), "a refused call opens nothing");
        app.attach_compare(&json!({"action": "open", "threshold": 1.0, "min_depth": 3}), &ctx).unwrap();
        app.goto("chr20:31180000");
        settle(&mut app, |a| {
            a.tracks.iter().all(|t| !t.is_loading() && t.current.as_ref().is_some_and(|c| c.view.region.start < 31_180_000 && c.view.region.end > 31_180_500))
        });
        app.compare.as_mut().unwrap().status = Some("No more differences in that direction".into());
        let r = app.attach_compare(&json!({"action": "next"}), &ctx).unwrap();
        assert_eq!(app.highlight.map(|h| h.1 + 1), Some(31_180_455), "a full difference, found at 100%");
        assert!(r["status"].as_str().is_some_and(|s| s.contains("31,180,455") && s.contains("100%")), "{}", r["status"]);
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The selection and outlined base name their contig by index: when the genome changes (a
    /// BAM opened in an annotation's window, the BAM closed again) they follow it by name, not
    /// becoming another contig (chr1) or one past the end (a panic).
    #[test]
    fn marks_follow_their_contig_when_the_genome_changes() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { genes: vec![data_dir().join("ncbiRefSeq.slices.gtf")], ..Default::default() }).unwrap();
        let t = std::time::Instant::now();
        while app.contigs.is_empty() && t.elapsed().as_secs() < 60 {
            app.poll_genes();
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert_eq!(app.contigs.len(), 1, "the annotation's genome: chr20 only");
        app.attach_goto("chr20:31180455").unwrap();
        app.attach_select("chr20:31180440-31180460").unwrap();
        let marks = |app: &App| {
            let s = app.attach_state();
            (s["selection"].clone(), s["highlight"].clone())
        };
        let want = (json!("chr20:31180440-31180460"), json!("chr20:31180455"));
        assert_eq!(marks(&app), want);
        app.add_tracks(vec![data_dir().join("HG002.revio.slices.bam")], &ctx);
        assert!(app.contigs.len() > 1 && app.contigs[app.contig].name == "chr20");
        assert_eq!(marks(&app), want, "a BAM's genome: chr20 is another index");
        app.remember_place();
        app.close_track(0);
        assert_eq!(app.contigs.len(), 1);
        assert_eq!(marks(&app), want, "back to the annotation's");
        assert!(!app.back.is_empty() && app.back.iter().all(|b| b.0 == 0));
        let r = app.attach_bookmarks(&json!({"action": "add", "note": "x"}));
        assert!(r.is_ok(), "{r:?}");
        assert!(app.bookmarks.borrow()[0].locus().starts_with("chr20:"));
    }

    /// Find read by name: a read loaded but off screen (reads load well beyond the view) is
    /// gone to, not reported on screen where it isn't; one on screen stays put.
    #[test]
    fn find_read_goes_to_a_loaded_read_off_screen() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { bams: vec![data_dir().join("HG002.revio.slices.bam")], ..Default::default() }).unwrap();
        app.goto("chr20:31180400-31180500");
        let name = "m84011_220902_175841_s1/237439347/ccs";
        settle(&mut app, |a| a.tracks[0].current.as_ref().and_then(|c| c.view.reads.as_ref()).is_some_and(|(rs, _)| rs.iter().any(|r| r.name == name)));
        let (vs, ve) = app.visible();
        // As Find does, and --find-read at startup.
        app.find = Some(super::super::tools::FindUi { text: name.into(), ..Default::default() });
        app.run_find_now(&ctx);
        let f = app.find.as_ref().unwrap();
        assert_eq!(f.loaded_hit.as_deref(), Some(name));
        assert!(f.search.is_none(), "the loaded copy says where");
        let (s, e) = app.visible();
        assert!((s, e) != (vs, ve) && s < 31_171_052.0 && e > 31_161_625.0, "now on screen: {s}-{e}");
        // On screen: nothing moves.
        settle(&mut app, |a| !a.tracks[0].is_loading() && a.tracks[0].current.as_ref().is_some_and(|c| (c.view.region.start as f64) < s));
        let seen = app.tracks[0]
            .current
            .as_ref()
            .and_then(|c| c.view.reads.as_ref())
            .and_then(|(rs, _)| rs.iter().find(|r| (r.start as f64) < e && (r.end as f64) > s && r.name != name))
            .map(|r| r.name.clone())
            .unwrap();
        app.find.as_mut().unwrap().text = seen.clone();
        app.run_find_now(&ctx);
        assert_eq!((app.find.as_ref().unwrap().loaded_hit.as_deref(), app.visible()), (Some(seen.as_str()), (s, e)));
    }

    /// `compare` and `rows` refuse a file name two open alignment tracks share (as `tracks`
    /// does), rather than taking the first; the full path picks one.
    #[test]
    fn compare_and_rows_refuse_an_ambiguous_file_name() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-ambiguous-bam-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let left = copy_data(&dir.join("left"), &["HG002.revio.slices.bam"]).remove(0);
        let right = copy_data(&dir.join("right"), &["HG002.revio.slices.bam"]).remove(0);
        let third = dir.join("third.bam");
        std::fs::copy(&left, &third).unwrap();
        std::fs::copy(left.with_added_extension("bai"), third.with_added_extension("bai")).unwrap();
        let ctx = egui::Context::default();
        let bams = vec![left.clone(), right.clone(), third];
        let mut app = App::new(&ctx, super::super::Args { bams, locus: Some("chr20:31180000-31181000".into()), ..Default::default() }).unwrap();
        app.attach_compare(&json!({"a": right.display().to_string(), "b": 3}), &ctx).unwrap();
        assert_eq!(app.compare.as_ref().map(|c| (c.a, c.b)), Some((1, 2)));
        let err = app.attach_compare(&json!({"a": "HG002.revio.slices.bam", "b": 3}), &ctx).unwrap_err();
        assert!(err.contains("ambiguous") && err.contains("left") && err.contains("right"), "{err}");
        assert_eq!(app.compare.as_ref().map(|c| (c.a, c.b)), Some((1, 2)), "unchanged");
        let err = app.attach_rows(&json!({"track": "HG002.revio.slices.bam", "to": "top"})).unwrap_err();
        assert!(err.contains("ambiguous"), "{err}");
        assert!(app.attach_rows(&json!({"track": "third.bam", "to": "top"})).is_ok());
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Signal tracks reload like the others: Reload all files (and the watcher, which lists
    /// them) reads the file again.
    #[test]
    fn signal_tracks_reload() {
        let dir = std::env::temp_dir().join(format!("ggv-signal-reload-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let sig = dir.join("s.bedgraph");
        std::fs::write(&sig, "chr20\t1000\t2000\t5\n").unwrap();
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { signals: vec![sig.clone()], ..Default::default() }).unwrap();
        let value = |app: &mut App| {
            let t = std::time::Instant::now();
            loop {
                let s = &mut app.signals[0];
                s.update("chr20", 1000.0, 2000.0, 100.0);
                if let Some(b) = s.at(0, 1500.0) {
                    return b.mean;
                }
                assert!(t.elapsed().as_secs() < 60, "timed out loading");
                std::thread::sleep(std::time::Duration::from_millis(10));
            }
        };
        assert_eq!(value(&mut app), 5.0);
        std::fs::write(&sig, "chr20\t1000\t2000\t9\n").unwrap();
        app.reload_all(&ctx);
        assert_eq!(value(&mut app), 9.0);
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Loads (without drawing, taking in genes and signal files as frames do) until `done`.
    fn settle_all(app: &mut App, done: impl Fn(&App) -> bool) {
        let t = std::time::Instant::now();
        while !done(app) && t.elapsed().as_secs() < 60 {
            app.poll_genes();
            app.poll_signals();
            app.load_hidden(1.0);
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(done(app), "timed out loading");
    }

    /// Closing the last alignment track keeps a signal track's genome (as a window opened with
    /// only signal tracks has): the view stays and can still be moved.
    #[test]
    fn closing_the_last_bam_keeps_the_signal_tracks_genome() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-close-last-bam-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let sig = dir.join("s.bedgraph");
        std::fs::write(&sig, "chr20\t31180000\t31190000\t5\n").unwrap();
        let ctx = egui::Context::default();
        let args = super::super::Args { bams: vec![data_dir().join("HG002.revio.slices.bam")], signals: vec![sig], ..Default::default() };
        let mut app = App::new(&ctx, args).unwrap();
        settle_all(&mut app, |a| !a.signals[0].contigs().is_empty());
        app.attach_goto("chr20:31181000-31182000").unwrap();
        let label = app.locus_label();
        app.close_track(0);
        assert!(!app.contigs.is_empty(), "the signal file's genome");
        assert_eq!(app.locus_label(), label);
        app.attach_goto("chr20:31185000-31186000").unwrap();
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// GNAS c.393C>T at chr20:58903752, with the BAM, the DeepVariant calls and RefSeq open.
    fn gnas_window(ctx: &egui::Context, reference: Option<std::path::PathBuf>) -> App {
        let args = super::super::Args {
            bams: vec![data_dir().join("HG002.revio.slices.bam")],
            vcfs: vec![data_dir().join("HG002.revio.deepvariant.phased.vcf.gz")],
            genes: vec![data_dir().join("ncbiRefSeq.slices.gtf")],
            reference,
            ..Default::default()
        };
        let mut app = App::new(ctx, args).unwrap();
        settle_all(&mut app, |a| !a.contigs.is_empty());
        app.attach_goto("chr20:58903752").unwrap();
        settle_all(&mut app, |a| a.variants.iter().all(|v| !v.is_loading() && v.current.is_some()) && a.genes.iter().all(|g| g.built.is_some()));
        app
    }

    /// Consequences cached before a reference was open are worked out again once it is: the
    /// codon is known then.
    #[test]
    fn consequences_use_a_reference_opened_later() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_window(&ctx, None);
        let lines = app.consequence_lines(58_903_751).join("\n");
        assert!(lines.contains("GNAS") && lines.contains("codon unknown"), "{lines}");
        app.set_reference(data_dir().join("GRCh38.chr20.slices.fa.gz"), &ctx);
        let lines = app.consequence_lines(58_903_751).join("\n");
        assert!(lines.contains("p.Ile131=") && !lines.contains("codon unknown"), "{lines}");
    }

    /// A variant's details window works out what the variant does at its own place, wherever
    /// the view goes while it is open (another contig, or out of the genes there).
    #[test]
    fn variant_details_consequences_stay_with_the_variant() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_window(&ctx, Some(data_dir().join("GRCh38.chr20.slices.fa.gz")));
        let var = app.variants[0].current.as_ref().unwrap().variants.iter().find(|v| v.start == 58_903_751).cloned().unwrap();
        app.open_variant_popup(0, var);
        let want = app.popup_consequences();
        assert!(want.iter().any(|l| l.contains("GNAS") && l.contains("p.Ile131=")), "{want:?}");
        for place in ["chr20:50000000", "chr1:58903000-58904500"] {
            app.attach_goto(place).unwrap();
            settle_all(&mut app, |a| a.genes.iter().all(|g| g.built.as_ref().is_some_and(|b| b.contig == a.contigs[a.contig].name)));
            // Worked out again (the cache would hide a view-following answer).
            app.variant_popup.as_mut().unwrap().consequences = None;
            assert_eq!(app.popup_consequences(), want, "at {place}");
        }
    }

    /// The evidence panel names transcripts without UCSC's copy suffix, as consequences do.
    #[test]
    fn evidence_genes_drop_ucsc_copy_suffixes() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-copy-suffix-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let gtf = dir.join("copy.gtf");
        let attrs = "gene_id \"PTEN\"; transcript_id \"NM_000314.8_2\"; gene_name \"PTEN\";";
        std::fs::write(&gtf, format!("chr20\tt\ttranscript\t31180001\t31181000\t.\t+\t.\t{attrs}\nchr20\tt\texon\t31180001\t31181000\t.\t+\t.\t{attrs}\n")).unwrap();
        let ctx = egui::Context::default();
        let args = super::super::Args { bams: vec![data_dir().join("HG002.revio.slices.bam")], genes: vec![gtf], ..Default::default() };
        let mut app = App::new(&ctx, args).unwrap();
        settle_all(&mut app, |a| !a.contigs.is_empty());
        app.attach_goto("chr20:31180455").unwrap();
        settle_all(&mut app, |a| a.genes.iter().all(|g| g.built.is_some()));
        assert_eq!(app.genes_at(31_180_454), vec!["PTEN NM_000314.8  exon 1/1".to_owned()]);
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// `gnas_window` (with these VCFs) with its reads loaded.
    fn gnas_with_reads(ctx: &egui::Context, vcfs: &[&str]) -> App {
        let args = super::super::Args {
            bams: vec![data_dir().join("HG002.revio.slices.bam")],
            vcfs: vcfs.iter().map(|v| data_dir().join(v)).collect(),
            genes: vec![data_dir().join("ncbiRefSeq.slices.gtf")],
            ..Default::default()
        };
        let mut app = App::new(ctx, args).unwrap();
        settle_all(&mut app, |a| !a.contigs.is_empty());
        app.attach_goto("chr20:58903752").unwrap();
        settle_all(&mut app, |a| {
            a.variants.iter().all(|v| !v.is_loading() && v.current.is_some())
                && a.genes.iter().all(|g| g.built.is_some())
                && a.tracks[0].current.as_ref().is_some_and(|c| c.view.reads.is_some() && !c.packed.rows.is_empty())
        });
        app
    }

    /// A window `h` points tall (the tracks start `TRACKS_TOP` below its top).
    fn window(h: f32) -> egui::Rect {
        egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1000.0, super::super::TRACKS_TOP + h))
    }

    /// Room the tracks leave goes to the alignment track, whose reads don't all fit, even when
    /// every track has a height the person set (it isn't the last track): no empty band below
    /// the tracks, no stretched track with nothing to show, and no edge to find below the last.
    #[test]
    fn the_last_track_fills_the_window() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_with_reads(&ctx, &["HG002.revio.deepvariant.phased.vcf.gz"]);
        for key in app.track_keys() {
            app.arrangement.heights.insert(key, 100.0);
        }
        let reads = app.content_height(Hit::Track(0)).unwrap().0;
        assert!(reads > 200.0, "enough reads to overflow: {reads}");
        let kinds = |l: &super::super::Layout| l.order.iter().map(|(k, r)| (k.0, r.height())).collect::<Vec<_>>();
        // Order: variants, alignments, genes (the last).
        let rect = window(200.0 + reads - 50.0);
        let layout = app.arrange(rect);
        use super::super::arrange::Kind;
        assert_eq!(kinds(&layout), vec![(Kind::Variants, 100.0), (Kind::Alignments, reads - 50.0), (Kind::Genes, 100.0)], "the reads take the room");
        assert_eq!(layout.order[2].1.bottom(), rect.bottom(), "no band below the last track");
        assert!(layout.stack.is_none());
        assert_eq!(layout.resize_edges(), 2, "the last track's bottom isn't a handle");
        assert!(layout.edge_at(egui::pos2(500.0, rect.bottom() - 1.0)).is_none());
        assert!(layout.edge_at(egui::pos2(500.0, layout.order[1].1.bottom())).is_some());

        // More room than the reads need: they get what they need; the gene track, with nothing
        // more to show, isn't stretched; what is left stays below, and isn't a handle either.
        let rect = window(200.0 + reads + 300.0);
        let layout = app.arrange(rect);
        assert_eq!(layout.order[1].1.height(), reads);
        let genes = app.content_height(Hit::Genes(0)).unwrap().0.max(100.0);
        assert_eq!(layout.order[2].1.height(), genes);
        assert!(layout.order[2].1.bottom() < rect.bottom());
        assert!(layout.edge_at(egui::pos2(500.0, layout.order[2].1.bottom())).is_none());

        // Zoomed out past the read window: the alignment track shows coverage and a note. A
        // height the person set stays (so the tracks below don't jump), and the genes, now with
        // more rows, use what is left; with none set it keeps no more room than it needs.
        app.attach_goto("chr20:10000000-60000000").unwrap();
        settle_all(&mut app, |a| a.genes.iter().all(|g| g.built.as_ref().is_some_and(|b| b.rows.len() > 1)));
        let (coverage, capped) = app.content_height(Hit::Track(0)).unwrap();
        assert!(capped);
        for k in app.track_keys() {
            let h = match k.0 {
                Kind::Alignments => 300.0,
                Kind::Genes => 40.0,
                _ => 100.0,
            };
            app.arrangement.heights.insert(k, h);
        }
        let gene_rows = app.content_height(Hit::Genes(0)).unwrap().0;
        assert!(gene_rows > 40.0);
        let layout = app.arrange(window(800.0));
        assert_eq!(layout.order[1].1.height(), 300.0, "the height set");
        assert_eq!(layout.order[2].1.height(), gene_rows.min(400.0), "the gene rows take the room: {:?}", kinds(&layout));
        let bam = layout.order[1].0.clone();
        app.arrangement.heights.remove(&bam);
        assert_eq!(app.arrange(window(800.0)).order[1].1.height(), coverage, "none set");
        for k in app.track_keys() {
            app.arrangement.heights.insert(k, 100.0);
        }

        // Too little room: heights give way, nothing is stretched.
        assert!(app.arrange(window(150.0)).order.iter().all(|(_, r)| r.height() <= 100.0));
    }

    /// With room below the last track, dragging the edge under the alignment track moves the
    /// genes down whole, back up again, and they stay put when zooming out hides the reads.
    #[test]
    fn dragging_the_edge_moves_the_genes_and_zooming_out_leaves_them() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_with_reads(&ctx, &[]);
        let reads = app.content_height(Hit::Track(0)).unwrap().0;
        let rect = window(reads + 600.0);
        let frame = |app: &mut App, events: Vec<egui::Event>| {
            let raw = egui::RawInput { events, screen_rect: Some(rect), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| {
                let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
                let layout = app.arrange(rect);
                app.arrange_input(ui, &response, &layout);
            });
            out.textures_delta.clear();
        };
        let drag = |app: &mut App, dy: f32| {
            let from = egui::pos2(500.0, app.arrange(rect).order[0].1.bottom());
            let to = from + egui::vec2(0.0, dy);
            let button = |pos, pressed| egui::Event::PointerButton { pos, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
            frame(app, vec![egui::Event::PointerMoved(from)]);
            frame(app, vec![button(from, true)]);
            frame(app, vec![egui::Event::PointerMoved(from + egui::vec2(0.0, dy / 2.0))]);
            frame(app, vec![egui::Event::PointerMoved(to)]);
            frame(app, vec![button(to, false)]);
            frame(app, vec![]);
        };
        let start = app.arrange(rect);
        let genes = start.order[1].1;
        assert!(genes.bottom() < rect.bottom() - 200.0, "room below the genes: {genes:?}");
        let moved = |app: &App| {
            let g = app.arrange(rect).order[1].1;
            assert!((g.height() - genes.height()).abs() < 0.5, "the same height: {g:?}, was {genes:?}");
            g.top() - genes.top()
        };

        drag(&mut app, 200.0);
        assert!((moved(&app) - 200.0).abs() < 0.5, "pushed down: {}", moved(&app));
        assert!(!app.arrangement.heights.contains_key(&start.order[1].0), "the genes keep their automatic height");
        drag(&mut app, -100.0);
        assert!((moved(&app) - 100.0).abs() < 0.5, "back up: {}", moved(&app));

        app.attach_goto("chr20:10000000-60000000").unwrap();
        assert!(app.content_height(Hit::Track(0)).unwrap().1, "reads hidden");
        assert!((moved(&app) - 100.0).abs() < 0.5, "where they were: {}", moved(&app));
    }

    /// In a small window, the alignment track (its reads scroll) gives way first; the small
    /// tracks keep their header and a whole row; only when even that doesn't fit does the stack
    /// scroll, and scrolled to the end, room left below the last tracks goes to the reads.
    #[test]
    fn small_windows_squeeze_the_reads_first() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_with_reads(&ctx, &["HG002.revio.deepvariant.phased.vcf.gz", "HG002.revio.pbsv.vcf.gz"]);
        let pbsv = app.variants[1].path.clone();
        use super::super::arrange::Kind;
        let (bam, gtf) = (app.tracks[0].loader.bam_path.clone(), app.genes[0].path.clone());
        let dv = app.variants[0].path.clone();
        // The person's session: heights set on a big screen.
        app.arrangement.order = vec![(Kind::Variants, dv.clone()), (Kind::Alignments, bam.clone()), (Kind::Variants, pbsv.clone()), (Kind::Genes, gtf.clone())];
        for (k, h) in app.arrangement.order.clone().into_iter().zip([63.0, 690.0, 67.0, 92.0]) {
            app.arrangement.heights.insert(k, h);
        }
        let heights = |l: &super::super::Layout| l.order.iter().map(|(_, r)| r.height()).collect::<Vec<_>>();
        let layout = app.arrange(window(330.0));
        assert_eq!(heights(&layout), vec![63.0, 108.0, 67.0, 92.0], "only the reads give way");
        assert!(layout.stack.is_none());

        // Smaller: the reads keep their coverage and two rows, the others give way together,
        // each keeping its header and a row.
        let hits = [Hit::Variants(0), Hit::Track(0), Hit::Variants(1), Hit::Genes(0)];
        let readable: Vec<f32> = hits.iter().map(|h| app.readable_height(*h)).collect();
        let room = readable.iter().sum::<f32>() + 10.0;
        let layout = app.arrange(window(room));
        assert!(layout.stack.is_none());
        assert_eq!(layout.order.last().unwrap().1.bottom(), window(room).bottom(), "what rounding leaves goes to the reads");
        for ((h, min), hit) in heights(&layout).into_iter().zip(&readable).zip(hits) {
            assert!(h >= *min, "{hit:?}: {h} < {min}");
        }
        assert!(heights(&layout)[1] - readable[1] < 4.0, "the reads gave way first: {:?}", heights(&layout));

        // Too small for that: the stack scrolls. The reads moved last; scrolled to the end, the
        // room the last two leave goes to the reads, not to a band below them.
        app.arrangement.order = vec![(Kind::Variants, dv), (Kind::Variants, pbsv), (Kind::Genes, gtf), (Kind::Alignments, bam)];
        app.stack_first = 99;
        let rect = window(readable[3] + readable[1] + 10.0);
        let layout = app.arrange(rect);
        let bar = layout.stack.expect("a track-stack scrollbar");
        assert_eq!(bar.first, bar.last_first);
        assert_eq!(layout.order.len(), 2);
        assert_eq!(layout.order[1].1.bottom(), rect.bottom(), "{:?}", heights(&layout));
    }

    /// Over reads (or gene or sample rows) that don't all fit, the wheel scrolls them, and the
    /// header and a scrollbar say where they are; over the coverage, the ruler, or a track
    /// whose rows fit, it zooms.
    #[test]
    fn the_wheel_scrolls_rows_that_dont_fit() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_with_reads(&ctx, &["HG002.revio.deepvariant.phased.vcf.gz"]);
        let layout = app.arrange(window(400.0));
        let section = layout.sections[0];
        let reads = egui::pos2(500.0, section.top() + crate::tracks::READS_TOP + 20.0);
        let (hit, rs) = app.rows_under(&layout, reads).expect("reads that scroll");
        assert_eq!(hit, Hit::Track(0));
        assert_eq!(rs.shown().0, 1);
        assert!(rs.shown().1 < rs.rows);
        assert!(rs.bar().right() <= section.right() && rs.thumb().top() >= rs.bar().top());
        assert!(app.rows_under(&layout, egui::pos2(500.0, section.top() + crate::tracks::COVERAGE_TOP + 20.0)).is_none(), "coverage zooms");
        assert!(app.rows_under(&layout, egui::pos2(500.0, 10.0)).is_none(), "the ruler zooms");
        app.scroll_rows(hit, 120.0);
        let rs = app.row_scroll(hit, section).unwrap();
        assert_eq!(rs.shown().0, 11);
        assert!(rs.thumb().top() > rs.bar().top());
        // A window tall enough for every read: nothing scrolls, the wheel zooms.
        let tall = app.arrange(window(2000.0));
        assert!(app.rows_under(&tall, egui::pos2(500.0, tall.sections[0].top() + crate::tracks::READS_TOP + 20.0)).is_none());
    }

    /// Structural-variant tags (and review ticks and clicks) come only from variants on the
    /// contig in view, in a view narrow enough to show variants, as the sites themselves.
    #[test]
    fn variant_tags_stay_with_their_contig_and_zoom() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_window(&ctx, None);
        assert!(app.shown_variants(0).is_some());
        app.attach_goto("chr20:10000000-60000000").unwrap();
        assert!(app.shown_variants(0).is_none(), "wider than the variant window");
        let other = app.contigs.iter().position(|c| c.name != "chr20").expect("another contig");
        app.contig = other;
        app.center = 20_000.0;
        app.ppb = 1.0;
        assert!(app.shown_variants(0).is_none(), "on another contig");
    }

    /// A gene track with nothing in view says what is in view, as the other tracks do.
    #[test]
    fn the_gene_track_says_what_is_in_view() {
        if !data_dir().join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut app = gnas_window(&ctx, None);
        assert!(app.genes_in_view(0).is_some_and(|n| n > 0));
        assert!(app.gene_header(0).contains(" in view"), "{}", app.gene_header(0));
        app.attach_goto("chr20:31165000-31184000").unwrap();
        settle_all(&mut app, |a| a.genes[0].built.as_ref().is_some_and(|b| b.start <= 31_165_000 && b.end >= 31_184_000));
        assert_eq!(app.genes_in_view(0), Some(0));
        assert!(app.gene_header(0).contains("·  0 in view  ·"), "{}", app.gene_header(0));
    }

    /// One frame of `draw` in a `w` × 900 window, as the window draws it.
    fn ui_frame(ctx: &egui::Context, w: f32, events: Vec<egui::Event>, draw: impl FnMut(&mut egui::Ui)) -> egui::FullOutput {
        let raw = egui::RawInput { events, screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(w, 900.0))), ..Default::default() };
        let mut out = ctx.run_ui(raw, draw);
        out.textures_delta.clear();
        out
    }

    /// The AccessKit nodes of a frame: role, label or value, and bounds.
    fn nodes(out: &egui::FullOutput) -> Vec<(egui::accesskit::Role, String, egui::Rect)> {
        let tree = out.platform_output.accesskit_update.as_ref().expect("accesskit enabled");
        tree.nodes
            .iter()
            .filter_map(|(_, n)| {
                let b = n.bounds()?;
                let text = n.label().or(n.value()).unwrap_or_default().to_owned();
                Some((n.role(), text, egui::Rect::from_min_max(egui::pos2(b.x0 as f32, b.y0 as f32), egui::pos2(b.x1 as f32, b.y1 as f32))))
            })
            .collect()
    }

    fn press(at: egui::Pos2) -> Vec<Vec<egui::Event>> {
        let button = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
        vec![vec![egui::Event::PointerMoved(at)], vec![button(true)], vec![button(false)]]
    }

    /// A status-bar error names the file once (ggv-io's error already starts with its path), its
    /// × comes first and stays on screen however long the text, the text is cut to the window,
    /// and clicking it copies the whole message.
    #[test]
    fn a_long_status_error_keeps_its_dismiss_and_copies_whole() {
        let path = std::path::Path::new("C:/Users/someone/AppData/Local/Temp/a/very/deep/folder/with/data/calls_plain.vcf");
        let io = format!("{}: no index found (looked for .tbi, .csi); a plain-text VCF can't be indexed: compress it with bgzip, then index that with tabix -p vcf", path.display());
        let e = super::super::path_error(path, &io);
        assert_eq!(e, io, "the path once");
        assert_eq!(e.matches("calls_plain.vcf").count(), 1);
        assert_eq!(super::super::path_error(path, "failed to fill whole buffer"), format!("{}: failed to fill whole buffer", path.display()));
        // On the bar the path is cut to the file's name, leaving room for the reason.
        assert!(super::super::status_line(&e).starts_with("calls_plain.vcf: no index found"), "{}", super::super::status_line(&e));
        // A path in this system's own spelling (Windows backslashes are no separator on Linux).
        let a = std::path::Path::new("data").join("a.bam");
        assert_eq!(super::super::status_line(&format!("{}: truncated\nb.bam: too", a.display())), "a.bam: truncated  ·  b.bam: too");
        assert_eq!(super::super::status_line("No more variants in that direction"), "No more variants in that direction");

        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        app.error = Some(e.clone());
        for dark in [false, true] {
            ctx.set_theme(if dark { egui::Theme::Dark } else { egui::Theme::Light });
            let mut out = None;
            for _ in 0..2 {
                out = Some(ui_frame(&ctx, 800.0, vec![], |ui| {
                    egui::Panel::bottom("status").show(ui, |ui| app.status_bar(ui, Some(0.001)));
                }));
            }
            let n = nodes(out.as_ref().unwrap());
            let dismiss = n.iter().find(|(r, t, _)| *r == egui::accesskit::Role::Button && t == "×").map(|x| x.2).expect("a × button");
            let label = n.iter().find(|(_, t, _)| t.contains("calls_plain")).map(|x| x.2).expect("the error's label");
            assert!(dismiss.right() <= 800.0 && dismiss.left() >= 0.0, "× on screen: {dismiss:?}");
            assert!(dismiss.left() < label.left(), "× before the text");
            assert!(label.right() <= 800.0, "the text is cut to the window: {label:?}");
            let frame = n.iter().find(|(_, t, _)| t.starts_with("frame ")).map(|x| x.2).expect("the frame time");
            assert!(frame.right() <= 800.0 && frame.left() > label.right(), "frame time still on screen, right of the error");
            let ink = super::super::error_rgb(dark);
            let bg = if dark { [27, 27, 27] } else { [248, 248, 248] };
            assert!(crate::variants::contrast([ink.r(), ink.g(), ink.b()], bg) >= 4.5, "readable red (dark: {dark})");
        }
        // A click on the text copies it whole; the × dismisses it.
        let out = ui_frame(&ctx, 800.0, vec![], |ui| {
            egui::Panel::bottom("status").show(ui, |ui| app.status_bar(ui, None));
        });
        let n = nodes(&out);
        let label = n.iter().find(|(_, t, _)| t.contains("calls_plain")).map(|x| x.2).unwrap();
        let mut copied = Vec::new();
        for events in press(label.center()) {
            let out = ui_frame(&ctx, 800.0, events, |ui| {
                egui::Panel::bottom("status").show(ui, |ui| app.status_bar(ui, None));
            });
            copied.extend(out.platform_output.commands.into_iter().filter_map(|c| match c {
                egui::OutputCommand::CopyText(t) => Some(t),
                _ => None,
            }));
        }
        assert_eq!(copied, vec![e.clone()]);
        assert!(app.error.is_some(), "copying keeps it");
        let dismiss = nodes(&ui_frame(&ctx, 800.0, vec![], |ui| {
            egui::Panel::bottom("status").show(ui, |ui| app.status_bar(ui, None));
        }))
        .into_iter()
        .find(|(r, t, _)| *r == egui::accesskit::Role::Button && t == "×")
        .unwrap()
        .2;
        for events in press(dismiss.center()) {
            ui_frame(&ctx, 800.0, events, |ui| {
                egui::Panel::bottom("status").show(ui, |ui| app.status_bar(ui, None));
            });
        }
        assert!(app.error.is_none(), "dismissed");
    }

    /// Files opened on their own that aren't there are listed in Files not found, as a
    /// session's are, rather than offered an index build they can't have; locating them opens
    /// them.
    #[test]
    fn a_file_that_isnt_there_is_reported_not_found_not_offered_an_index() {
        let data = data_dir();
        if !data.join("HG002.revio.slices.bam").exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-not-found-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let (vcf, bam) = (dir.join("missing.vcf.gz"), dir.join("moved.bam"));
        assert!(!super::super::missing_index(&vcf), "nothing to index");
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args { bams: vec![bam.clone()], vcfs: vec![vcf.clone()], ..Default::default() }).unwrap();
        assert!(app.indexing.is_empty(), "no index offered");
        assert!(app.tracks.is_empty() && app.variants.is_empty());
        let m = app.missing.as_ref().expect("listed in Files not found");
        assert_eq!(m.paths, vec![bam.clone(), vcf.clone()]);
        assert!(m.opened);
        assert!(!app.recent.contains(&vcf) && !app.recent.contains(&bam), "not remembered as recent");
        // Dropped again: listed once.
        app.open_paths(vec![vcf.clone()], &ctx);
        assert_eq!(app.missing.as_ref().unwrap().paths.len(), 2);
        assert!(app.indexing.is_empty());

        // Found in a folder: both open, beside what's open.
        for (from, to) in [("HG002.revio.slices.bam", &bam), ("HG002.revio.deepvariant.phased.vcf.gz", &vcf)] {
            std::fs::copy(data.join(from), to).unwrap();
        }
        std::fs::copy(data.join("HG002.revio.slices.bam.bai"), bam.with_added_extension("bai")).unwrap();
        std::fs::copy(data.join("HG002.revio.deepvariant.phased.vcf.gz.tbi"), vcf.with_added_extension("tbi")).unwrap();
        app.look_in(&dir, &ctx);
        assert!(app.missing.is_none(), "{:?}", app.missing.as_ref().map(|m| m.paths.clone()));
        assert_eq!((app.tracks.len(), app.variants.len()), (1, 1));
        assert!(app.error.is_none(), "{:?}", app.error);
        drop(app);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// After a window resize with no input, the locus box shows the region now on screen: the
    /// frame that finds the new width asks for another, which the toolbar draws with it.
    #[test]
    fn the_locus_box_follows_a_resize() {
        let ctx = egui::Context::default();
        let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
        app.contigs = vec![ggv_model::Contig { name: "chr20".into(), length: 64_444_167 }];
        app.contig = 0;
        app.pending_goto = None;
        app.center = 31_170_000.0;
        app.ppb = 0.5;
        let frame = |app: &mut App, w: f32| {
            let c = ctx.clone();
            ui_frame(&ctx, w, vec![], |ui| {
                egui::Panel::top("toolbar").show(ui, |ui| app.toolbar(ui, &c));
                egui::CentralPanel::default().show(ui, |ui| app.body(ui));
            })
        };
        frame(&mut app, 1400.0);
        frame(&mut app, 1400.0);
        assert_eq!(Some(app.locus_text.clone()), app.locus_now());
        let out = frame(&mut app, 900.0);
        assert_ne!(Some(app.locus_text.clone()), app.locus_now(), "this frame's toolbar came before the new width");
        assert_eq!(out.viewport_output[&egui::ViewportId::ROOT].repaint_delay, std::time::Duration::ZERO, "so it asks for another frame");
        frame(&mut app, 900.0);
        assert_eq!(Some(app.locus_text.clone()), app.locus_now(), "which shows the region on screen");
    }

    /// The toolbar keeps the whole locus visible and stays on one row at 800, 1024 and 1366
    /// points wide, folding the chips and then Split, Evidence and Commands as it narrows.
    #[test]
    fn the_toolbar_keeps_the_whole_locus_at_narrow_widths() {
        use super::super::{ChipFit, fit_toolbar};
        // (row, fixed, locus need, chips each/folded/short, tail full/More)
        for row in [400.0, 784.0, 1008.0, 1350.0, 2400.0] {
            let fit = fit_toolbar(row, 420.0, 40.0, 210.0, [400.0, 140.0, 80.0], [220.0, 50.0]);
            assert!(fit.locus_w >= 210.0 && fit.locus_w <= 520.0, "{row}: {fit:?}");
        }
        assert_eq!(fit_toolbar(2400.0, 420.0, 40.0, 210.0, [400.0, 140.0, 80.0], [220.0, 50.0]).chips, ChipFit::Each);
        let narrow = fit_toolbar(820.0, 420.0, 40.0, 210.0, [400.0, 140.0, 80.0], [220.0, 50.0]);
        assert_eq!((narrow.chips, narrow.more, narrow.go), (ChipFit::Short, true, true));
        let narrower = fit_toolbar(760.0, 420.0, 40.0, 210.0, [400.0, 140.0, 80.0], [220.0, 50.0]);
        assert_eq!((narrower.chips, narrower.more, narrower.go), (ChipFit::Short, true, false), "Go goes last");
        // A locus longer than the box's usual maximum still fits whole.
        assert_eq!(fit_toolbar(2400.0, 420.0, 40.0, 600.0, [0.0; 3], [220.0, 50.0]).locus_w, 600.0);

        let data = data_dir();
        let vcf = data.join("HG002.revio.deepvariant.phased.vcf.gz");
        if !vcf.exists() {
            return;
        }
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let mut app = App::new(&ctx, super::super::Args { vcfs: vec![vcf], ..Default::default() }).unwrap();
        app.contigs = vec![ggv_model::Contig { name: "chr20".into(), length: 64_444_167 }];
        app.contig = 0;
        app.pending_goto = None;
        app.center = 31_174_500.0;
        app.ppb = 0.3;
        app.width = 1000.0;
        // Three settings changed, as in the session the bug was seen with.
        app.style.squished = true;
        app.style.soft_clips = true;
        app.style.junctions = true;
        for w in [800.0, 1024.0, 1366.0] {
            let c = ctx.clone();
            let mut height = 0.0;
            let mut out = None;
            for _ in 0..3 {
                out = Some(ui_frame(&ctx, w, vec![], |ui| {
                    height = egui::Panel::top("toolbar").show(ui, |ui| app.toolbar(ui, &c)).response.rect.height();
                }));
            }
            let n = nodes(out.as_ref().unwrap());
            let locus = app.locus_now().unwrap();
            let edit = n.iter().find(|(r, _, _)| *r == egui::accesskit::Role::TextInput).map(|x| x.2).expect("the locus box");
            let mono = egui::TextStyle::Monospace.resolve(&ctx.global_style());
            let text_w = ctx.fonts_mut(|f| f.layout_no_wrap(locus.clone(), mono, egui::Color32::WHITE).size().x);
            assert!(edit.width() >= text_w + 8.0, "{w}: the box ({}) holds {locus} ({text_w})", edit.width());
            assert!(edit.right() <= w, "{w}: on screen");
            let row = ctx.global_style().spacing.interact_size.y;
            assert!(height < row * 2.0 + 4.0, "{w}: one row ({height})");
            let labels: Vec<&str> = n.iter().map(|x| x.1.as_str()).collect();
            assert!(labels.contains(&"View") && labels.contains(&"Open"), "{w}: the menus stay");
            let commands = labels.contains(&"Commands") || labels.contains(&"More");
            assert!(commands, "{w}: Commands, or More holding it");
        }
    }

    /// The status bar counts every track on screen by kind, not just the alignments.
    #[test]
    fn the_status_bar_counts_every_track() {
        use super::super::tracks_text;
        assert_eq!(tracks_text(1, 0, 0, 0), "1 alignment track");
        assert_eq!(tracks_text(0, 2, 0, 0), "2 variant tracks");
        assert_eq!(tracks_text(1, 2, 1, 0), "4 tracks (1 alignment, 2 variant, 1 gene)");
        assert_eq!(tracks_text(0, 0, 1, 1), "2 tracks (1 gene, 1 signal)");
        let data = data_dir();
        if !data.join("HG002.revio.slices.bam").exists() {
            return;
        }
        let ctx = egui::Context::default();
        let args = super::super::Args {
            bams: vec![data.join("HG002.revio.slices.bam")],
            vcfs: vec![data.join("HG002.revio.deepvariant.phased.vcf.gz"), data.join("HG002.revio.pbsv.vcf.gz")],
            genes: vec![data.join("ncbiRefSeq.slices.gtf")],
            ..Default::default()
        };
        let app = App::new(&ctx, args).unwrap();
        let text = app.counts_text();
        assert!(text.starts_with("4 tracks (1 alignment, 2 variant, 1 gene) · "), "{text}");
    }
}

