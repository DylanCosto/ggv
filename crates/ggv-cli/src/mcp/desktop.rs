//! `ggv mcp --attach`: tools answer about what the person is looking at in an GGV desktop window.
//!
//! Before each call the window's state is read (a few hundred microseconds over loopback), and
//! whatever the call leaves out is filled in from it: the files, the base (the outlined base,
//! else the view centre: what the app's Evidence panel shows), the region (the selection, else
//! the view) and each file's read filter. Every answer says what was taken from the window, so
//! the agent, and the person reading along, know which state produced the numbers. Arguments the
//! call does give always win.
//!
//! With several windows open, the one in use is the focused window, else the one used most
//! recently, unless `--attach PID` or a call's `window` names one.

use std::{
    collections::HashMap,
    net::{Shutdown, TcpStream},
    path::PathBuf,
    sync::{Arc, Mutex},
};

use ggv_io::attach::{Client, Endpoint, endpoints, endpoints_in};
use serde_json::{Value, json};

use super::tools::Defaults;

pub struct Desktop {
    /// `--attach PID`: always this window.
    pinned: Option<u32>,
    clients: HashMap<u32, Client>,
    /// Where endpoint files are (tests); else the usual folder.
    dir: Option<PathBuf>,
}

/// Lets another thread stop a wait on a window: cancelling (or the client going away) shuts the
/// connection the wait is reading, so it ends at once.
#[derive(Clone, Default)]
pub struct Cancel(Arc<Mutex<Cancelling>>);

#[derive(Default)]
struct Cancelling {
    cancelled: bool,
    /// The client has gone: a wait still going is stopped, but a finished one is kept.
    closing: bool,
    /// A wait was cut short by either.
    interrupted: bool,
    waiting_on: Option<TcpStream>,
}

impl Cancel {
    pub fn cancel(&self) {
        self.stop(true);
    }

    /// Stops a wait still going, without cancelling one that has finished.
    pub fn close(&self) {
        self.stop(false);
    }

    fn stop(&self, cancel: bool) {
        let mut c = self.0.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        c.cancelled |= cancel;
        c.closing = true;
        if let Some(s) = c.waiting_on.take() {
            let _ = s.shutdown(Shutdown::Both);
            c.interrupted = true;
        }
    }

    /// Whether `other` is this one (a clone of it).
    pub fn same(&self, other: &Cancel) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.lock().unwrap_or_else(std::sync::PoisonError::into_inner).cancelled
    }

    /// Whether a wait was cut short (cancelled, or the client went away).
    pub fn interrupted(&self) -> bool {
        self.0.lock().unwrap_or_else(std::sync::PoisonError::into_inner).interrupted
    }

    /// Notes the connection a wait is about to read; false if already stopped.
    fn watch(&self, stream: Option<TcpStream>) -> bool {
        let mut c = self.0.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        if c.closing {
            c.interrupted = true;
            return false;
        }
        c.waiting_on = stream;
        true
    }
}

/// A window's state, as read for one call.
pub struct Window {
    pub pid: u32,
    pub state: Value,
}

impl Window {
    /// The pane a call is about: `pane` (1-based) or the active one.
    pub fn pane(&self, pane: Option<u64>) -> Result<&Value, String> {
        let panes = self.state["panes"].as_array().ok_or("the window sent no panes")?;
        let n = pane.unwrap_or_else(|| self.state["active_pane"].as_u64().unwrap_or(1));
        panes.get((n as usize).wrapping_sub(1)).ok_or_else(|| format!("no pane {n}: the window has {}", panes.len()))
    }
}

impl Desktop {
    pub fn new(pinned: Option<u32>) -> Self {
        Self { pinned, clients: HashMap::new(), dir: None }
    }

    /// Reads endpoint files from `dir` instead of the usual folder.
    #[cfg(test)]
    pub fn in_dir(pinned: Option<u32>, dir: PathBuf) -> Self {
        Self { pinned, clients: HashMap::new(), dir: Some(dir) }
    }

    /// A new `Desktop` for the same windows, with its own connections (for another thread).
    pub fn fresh(&self) -> Self {
        Self { pinned: self.pinned, clients: HashMap::new(), dir: self.dir.clone() }
    }

    fn endpoints(&self) -> Vec<Endpoint> {
        match &self.dir {
            Some(d) => endpoints_in(d),
            None => endpoints(),
        }
    }

    fn request(&mut self, pid: u32, method: &str, params: Value) -> Result<Value, String> {
        // A kept connection may have gone stale (the window restarted or closed): once more on a
        // fresh one, which removes the endpoint file if the window is gone. Only then: an error
        // the window sent back (or a reply that timed out) means it has the request, and sending
        // it again would carry out a command twice. A timed-out connection is still dropped, so
        // its late reply can't answer the next request.
        if self.clients.contains_key(&pid) {
            match self.clients.get_mut(&pid).expect("present").request(method, params.clone()) {
                Ok(v) => return Ok(v),
                Err(e) if e.stale => {
                    self.clients.remove(&pid);
                }
                Err(e) => {
                    if e.lost {
                        self.clients.remove(&pid);
                    }
                    return Err(e.message);
                }
            }
        }
        {
            let ep = self.endpoints().into_iter().find(|e| e.pid == pid).ok_or_else(|| format!("no GGV window with process id {pid} is open"))?;
            let c = Client::connect(&ep).map_err(|e| format!("window {pid}: {e}"))?;
            self.clients.insert(pid, c);
        }
        let r = self.clients.get_mut(&pid).expect("just connected").request(method, params);
        if r.as_ref().is_err_and(|e| e.lost) {
            // Reconnect next time (the window may have restarted or gone, or a late reply may
            // still be on its way).
            self.clients.remove(&pid);
        }
        r.map_err(|e| e.message)
    }

    /// Every open window's state (unreachable ones are skipped and their files removed).
    pub fn windows(&mut self) -> Vec<Window> {
        let pids: Vec<u32> = self.endpoints().iter().map(|e| e.pid).collect();
        pids.into_iter().filter_map(|pid| self.request(pid, "state", Value::Null).ok().map(|state| Window { pid, state })).collect()
    }

    /// The window a call is about: `window` if given, else `--attach PID`, else the focused or
    /// most recently used one.
    pub fn window(&mut self, window: Option<u32>) -> Result<Window, String> {
        if let Some(pid) = window.or(self.pinned) {
            return self.request(pid, "state", Value::Null).map(|state| Window { pid, state });
        }
        let mut all = self.windows();
        if all.is_empty() {
            return Err("no GGV desktop window is open: start ggv-desktop (and open files), or pass files to the call".into());
        }
        all.sort_by_key(|w| (!w.state["focused"].as_bool().unwrap_or(false), w.state["idle_ms"].as_u64().unwrap_or(u64::MAX)));
        Ok(all.swap_remove(0))
    }

    pub fn command(&mut self, pid: u32, method: &str, params: Value) -> Result<Value, String> {
        self.request(pid, method, params)
    }

    /// Waits up to `wait` for the person's next events (on a fresh request, so a stale
    /// connection is replaced first). `cancel` stops the wait early.
    pub fn events(&mut self, pid: u32, params: Value, wait: std::time::Duration, cancel: &Cancel) -> Result<Value, String> {
        self.request(pid, "state", Value::Null)?;
        let c = self.clients.get_mut(&pid).ok_or("no connection to the window")?;
        if !cancel.watch(c.interrupter().ok()) {
            return Err("cancelled".into());
        }
        let r = c.request_waiting("events", params, wait + std::time::Duration::from_secs(5));
        cancel.watch(None);
        if r.as_ref().is_err_and(|e| e.lost) {
            self.clients.remove(&pid);
        }
        r.map_err(|e| if cancel.interrupted() { "cancelled".into() } else { e.message })
    }

    /// The `events` tool: the person's next events in `window` (else the usual one).
    pub fn events_tool(&mut self, window: Option<u32>, args: &Value, cancel: &Cancel) -> Result<Value, String> {
        let wait_s = args.get("wait_s").and_then(Value::as_f64).unwrap_or(30.0).clamp(0.0, 300.0);
        let wait = std::time::Duration::from_secs_f64(wait_s);
        let mut params = json!({"wait_ms": wait.as_millis() as u64});
        for k in ["since", "kinds"] {
            if let Some(v) = args.get(k) {
                params[k] = v.clone();
            }
        }
        let w = self.window(window)?;
        let mut out = self.events(w.pid, params, wait, cancel)?;
        out["window"] = json!(w.pid);
        Ok(out)
    }
}

fn strings(v: &Value) -> Vec<PathBuf> {
    v.as_array().map(|a| a.iter().filter_map(Value::as_str).map(PathBuf::from).collect()).unwrap_or_default()
}

/// The pane's files as tool defaults.
pub fn defaults_of(pane: &Value) -> Defaults {
    let f = &pane["files"];
    Defaults { bams: strings(&f["bams"]), reference: f["reference"].as_str().map(PathBuf::from), vcfs: strings(&f["vcfs"]), genes: strings(&f["genes"]), signals: strings(&f["signals"]) }
}

/// Tools that take a base, a region, read filters, and both a base and a region.
const BASE_TOOLS: &[&str] = &["evidence", "consequence"];
const REGION_TOOLS: &[&str] = &["summarize", "coverage", "variants", "genes", "phase", "compare", "signal", "junctions", "svs", "repeats"];
const FILTER_TOOLS: &[&str] = &["evidence", "reads", "summarize", "coverage", "phase", "compare", "render", "junctions", "svs", "repeats"];
const FILTER_KEYS: &[&str] = &["filter", "min_mapq", "include_duplicates", "exclude_flags"];

fn has(args: &Value, key: &str) -> bool {
    args.get(key).is_some_and(|v| !v.is_null() && v.as_str() != Some(""))
}

/// Fills in what `args` leaves out from the pane, returning the new arguments and a record of
/// what was taken from it.
pub fn fill(tool: &str, args: &Value, pane: &Value) -> (Value, serde_json::Map<String, Value>, Vec<String>) {
    let mut a = if args.is_object() { args.clone() } else { json!({}) };
    let mut applied = serde_json::Map::new();
    let mut notes = Vec::new();
    // The outlined base while it is on screen (one panned away from is no longer what the
    // person looks at), else the view centre.
    let v = &pane["view"];
    let on_screen = |h: &str| -> bool {
        let Some((contig, pos)) = h.rsplit_once(':') else { return false };
        let pos: u64 = pos.replace(',', "").parse().unwrap_or(0);
        Some(contig) == v["contig"].as_str() && v["start"].as_u64().is_some_and(|s| s <= pos) && v["end"].as_u64().is_some_and(|e| pos <= e)
    };
    let base = pane["highlight"]
        .as_str()
        .filter(|h| on_screen(h))
        .map(|h| (h.to_owned(), "outlined base"))
        .or_else(|| Some((format!("{}:{}", v["contig"].as_str()?, v["centre"].as_u64()?), "view centre")));
    // A short selection made last: every base in it (up to 50, as `evidence` takes).
    let selected_bases: Option<Vec<String>> = (pane["latest"] == json!("selection"))
        .then(|| pane["selection"].as_str())
        .flatten()
        .and_then(|s| {
            let (contig, span) = s.rsplit_once(':')?;
            let (a, b) = span.split_once('-')?;
            let (a, b): (u64, u64) = (a.replace(',', "").parse().ok()?, b.replace(',', "").parse().ok()?);
            (b >= a && b - a < 50).then(|| (a..=b).map(|p| format!("{contig}:{p}")).collect())
        });
    let region = pane["selection"].as_str().map(|s| (s.to_owned(), "selection")).or_else(|| pane["view"]["region"].as_str().map(|r| (r.to_owned(), "view")));

    let set = |a: &mut Value, applied: &mut serde_json::Map<String, Value>, key: &str, (value, from): (String, &str)| {
        a[key] = json!(value);
        applied.insert(key.into(), json!({"value": value, "from": from}));
    };
    let gave_place = has(&a, "locus") || has(&a, "loci") || has(&a, "region");
    if BASE_TOOLS.contains(&tool) && !has(&a, "locus") && !has(&a, "loci") {
        match (&selected_bases, base.clone()) {
            (Some(loci), _) if tool == "evidence" => {
                a["loci"] = json!(loci);
                applied.insert("loci".into(), json!({"value": pane["selection"], "from": "selection (each base)"}));
            }
            (_, Some(b)) => set(&mut a, &mut applied, "locus", b),
            _ => {}
        }
    } else if REGION_TOOLS.contains(&tool) && !has(&a, "region") && !(tool == "genes" && has(&a, "name")) {
        if let Some(r) = region.clone() {
            set(&mut a, &mut applied, "region", r);
        }
    } else if tool == "reads" && !gave_place {
        // Whichever the person pointed at last: the selection or the outlined base.
        let selection = pane["selection"].as_str();
        let base_last = pane["latest"] != json!("selection") && !pane["highlight"].is_null();
        match (selection, base.clone()) {
            (Some(s), _) if !base_last => set(&mut a, &mut applied, "region", (s.to_owned(), "selection")),
            (_, Some(b)) => set(&mut a, &mut applied, "locus", b),
            _ => {}
        }
    } else if tool == "render" && !gave_place {
        // A picture of what the person sees.
        if let Some(r) = pane["view"]["region"].as_str() {
            set(&mut a, &mut applied, "region", (r.to_owned(), "view"));
        }
    }
    if tool == "render" && !has(&a, "color_by") {
        let color = match pane["display"]["color"].as_str() {
            Some("strand") => Some("strand"),
            Some("haplotype") => Some("haplotype"),
            Some("insert size") => Some("insert-size"),
            Some("base modifications") => Some("base-mods"),
            Some(t) if t.starts_with("tag:") => Some(t),
            _ => None,
        };
        if let Some(c) = color {
            a["color_by"] = json!(c);
            applied.insert("color_by".into(), json!({"value": c, "from": "display"}));
        }
    }
    if tool == "render" && !has(&a, "link") && pane["display"]["link_splits"] == json!(true) {
        a["link"] = json!(true);
        applied.insert("link".into(), json!({"value": true, "from": "display"}));
    }
    // Grouping too, except by base (a picture of the view has no single base to group at).
    if tool == "render" && !has(&a, "group_by")
        && let Some(g) = pane["display"]["group"].as_str().filter(|g| ["strand", "haplotype", "read group", "alleles"].contains(g) || g.starts_with("tag:"))
    {
        a["group_by"] = json!(g);
        applied.insert("group_by".into(), json!({"value": g, "from": "display"}));
    }

    // Files: the defaults are the pane's for the call; say which were used.
    let d = defaults_of(pane);
    let names_files = ["bams", "bam", "a", "b", "files"].iter().any(|k| has(&a, k));
    if !names_files && !d.bams.is_empty() && tool != "variants" && tool != "genes" && tool != "consequence" {
        applied.insert("bams".into(), json!(d.bams));
    }
    if !has(&a, "reference") && d.reference.is_some() && tool != "variants" && tool != "genes" {
        applied.insert("reference".into(), json!(d.reference));
    }
    if !has(&a, "vcfs") && !d.vcfs.is_empty() && ["evidence", "variants", "svs"].contains(&tool) {
        applied.insert("vcfs".into(), json!(d.vcfs));
    }
    if !has(&a, "genes") && !d.genes.is_empty() && ["evidence", "genes", "consequence"].contains(&tool) {
        applied.insert("genes".into(), json!(d.genes));
    }
    if !has(&a, "signals") && !d.signals.is_empty() && tool == "signal" {
        applied.insert("signals".into(), json!(d.signals));
    }

    // Read filters, as the app applies them to the files the call is about.
    if FILTER_TOOLS.contains(&tool) && !FILTER_KEYS.iter().any(|k| has(&a, k)) {
        let mut used: Vec<PathBuf> = ["bams", "bam", "a", "b"].iter().flat_map(|k| match &a[*k] {
            Value::String(s) => vec![PathBuf::from(s)],
            v => strings(v),
        }).collect();
        if used.is_empty() {
            used = match tool {
                "reads" | "phase" => d.bams.iter().take(1).cloned().collect(),
                "compare" => d.bams.iter().take(2).cloned().collect(),
                _ => d.bams.clone(),
            };
        }
        let filters = pane["filters"].as_array().cloned().unwrap_or_default();
        let of = |p: &PathBuf| filters.iter().find(|f| f["path"].as_str().is_some_and(|x| std::path::Path::new(x) == p.as_path())).cloned();
        let found: Vec<Value> = used.iter().filter_map(of).collect();
        let differ = |f: &Value, g: &Value| f["exclude_flags"] != g["exclude_flags"] || f["min_mapq"] != g["min_mapq"] || f["expression"] != g["expression"];
        if tool == "compare" && used.len() == 2 && found.len() == 2 && differ(&found[0], &found[1]) {
            // The window compares each track under its own filter; so does the tool.
            let own = |f: &Value| json!({"exclude_flags": f["exclude_flags"], "min_mapq": f["min_mapq"], "filter": f["expression"].as_str().unwrap_or("")});
            a["a_filters"] = own(&found[0]);
            a["b_filters"] = own(&found[1]);
            let shown = |f: &Value| json!({"exclude_flags": f["exclude_flags"], "min_mapq": f["min_mapq"], "expression": f["expression"]});
            applied.insert("filters".into(), json!({"a": shown(&found[0]), "b": shown(&found[1])}));
            notes.push("the two files have different read filters in the app; each was counted with its own, as the window's comparison does".into());
        } else if let Some(first) = found.first() {
            if found.iter().any(|f| differ(f, first)) {
                notes.push(format!(
                    "the files have different read filters in the app; all were counted with {}'s; pass one file at a time for each file's own",
                    first["path"].as_str().unwrap_or("the first file")
                ));
            }
            if found.len() < used.len() {
                notes.push("some files aren't open in the app; they were counted with its filter for the others".into());
            }
            a["exclude_flags"] = first["exclude_flags"].clone();
            a["min_mapq"] = first["min_mapq"].clone();
            if let Some(e) = first["expression"].as_str() {
                a["filter"] = json!(e);
            }
            applied.insert("filters".into(), json!({"exclude_flags": first["exclude_flags"], "min_mapq": first["min_mapq"], "expression": first["expression"]}));
        }
    }
    if pane["ready"] == json!(false) {
        notes.push("the window was still loading; the app's own numbers may not be on screen yet".into());
    }
    (a, applied, notes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pane() -> Value {
        json!({
            "view": {"region": "chr20:1,000-2,000", "contig": "chr20", "start": 1000, "end": 2000, "centre": 1500},
            "highlight": "chr20:1234",
            "selection": null,
            "files": {"bams": ["/d/a.bam", "/d/b.bam"], "reference": "/d/ref.fa", "vcfs": ["/d/v.vcf.gz"], "genes": []},
            "filters": [
                {"path": "/d/a.bam", "exclude_flags": "0x704", "min_mapq": 20, "expression": "!supplementary"},
                {"path": "/d/b.bam", "exclude_flags": "0x704", "min_mapq": 20, "expression": "!supplementary"}
            ],
            "ready": true,
        })
    }

    /// The window's comparison counts each track under its own filter; so does the tool, where
    /// it used to count both under the first file's and find nothing.
    #[test]
    fn compare_takes_each_tracks_own_filter() {
        let mut p = pane();
        p["filters"][0]["expression"] = json!("hp == 1");
        p["filters"][1]["expression"] = json!("hp == 2");
        let (a, applied, notes) = fill("compare", &json!({}), &p);
        assert_eq!((a["a_filters"]["filter"].as_str(), a["b_filters"]["filter"].as_str()), (Some("hp == 1"), Some("hp == 2")));
        assert_eq!(a["a_filters"]["min_mapq"], 20);
        assert!(a.get("filter").is_none(), "no one filter for both: {a}");
        assert_eq!(applied["filters"]["b"]["expression"], "hp == 2");
        assert!(notes.iter().any(|n| n.contains("each was counted with its own")), "{notes:?}");
        // The same filter on both: one filter, as for any other tool.
        let (a, _, notes) = fill("compare", &json!({}), &pane());
        assert!(a.get("a_filters").is_none() && a["filter"] == "!supplementary" && notes.is_empty(), "{a} {notes:?}");
        // Another tool over the two files still says which filter it used for both.
        let (_, _, notes) = fill("summarize", &json!({}), &p);
        assert!(notes.iter().any(|n| n.contains("different read filters")), "{notes:?}");
    }

    #[test]
    fn fills_base_region_files_and_filters() {
        let (a, applied, notes) = fill("evidence", &json!({}), &pane());
        assert_eq!(a["locus"], "chr20:1234");
        assert_eq!(applied["locus"]["from"], "outlined base");
        assert_eq!(a["min_mapq"], 20);
        assert_eq!(a["exclude_flags"], "0x704");
        assert_eq!(a["filter"], "!supplementary");
        assert_eq!(applied["bams"], json!(["/d/a.bam", "/d/b.bam"]));
        assert!(notes.is_empty(), "{notes:?}");

        let (a, applied, _) = fill("summarize", &json!({}), &pane());
        assert_eq!(a["region"], "chr20:1,000-2,000");
        assert_eq!(applied["region"]["from"], "view");

        let mut p = pane();
        p["selection"] = json!("chr20:1200-1300");
        p["highlight"] = Value::Null;
        let (a, _, _) = fill("reads", &json!({}), &p);
        assert_eq!(a["region"], "chr20:1200-1300");
        // Both set: whichever came last.
        p["highlight"] = json!("chr20:1250");
        p["latest"] = json!("highlight");
        assert_eq!(fill("reads", &json!({}), &p).0["locus"], "chr20:1250");
        p["latest"] = json!("selection");
        assert_eq!(fill("reads", &json!({}), &p).0["region"], "chr20:1200-1300");
        // A picture is of the view, selection or not.
        assert_eq!(fill("render", &json!({}), &p).0["region"], "chr20:1,000-2,000");
        p["highlight"] = Value::Null;
        let (a, applied, _) = fill("evidence", &json!({}), &p);
        assert_eq!(a["locus"], "chr20:1500");
        assert_eq!(applied["locus"]["from"], "view centre");
    }

    #[test]
    fn base_follows_what_is_on_screen() {
        let mut p = pane();
        // Outlined base panned off screen: the view centre instead.
        p["highlight"] = json!("chr20:5000");
        let (a, applied, _) = fill("evidence", &json!({}), &p);
        assert_eq!(a["locus"], "chr20:1500");
        assert_eq!(applied["locus"]["from"], "view centre");
        // A short selection made last: each of its bases.
        p["selection"] = json!("chr20:1200-1203");
        p["latest"] = json!("selection");
        let (a, _, _) = fill("evidence", &json!({}), &p);
        assert_eq!(a["loci"], json!(["chr20:1200", "chr20:1201", "chr20:1202", "chr20:1203"]));
        // A long one: the base rules again.
        p["selection"] = json!("chr20:1200-1400");
        assert_eq!(fill("evidence", &json!({}), &p).0["locus"], "chr20:1500");
    }

    #[test]
    fn pictures_use_the_window_colouring() {
        let mut p = pane();
        p["display"] = json!({"color": "base modifications"});
        let (a, applied, _) = fill("render", &json!({}), &p);
        assert_eq!(a["color_by"], "base-mods");
        assert_eq!(applied["color_by"]["from"], "display");
        assert_eq!(a["region"], "chr20:1,000-2,000");
        // Asked for explicitly: that wins.
        assert_eq!(fill("render", &json!({"color_by": "strand"}), &p).0["color_by"], "strand");
        p["display"] = json!({"color": "none", "group": "haplotype"});
        let (a, _, _) = fill("render", &json!({}), &p);
        assert!(a.get("color_by").is_none());
        assert_eq!(a["group_by"], "haplotype");
    }

    #[test]
    fn given_arguments_win() {
        let (a, applied, _) = fill("evidence", &json!({"locus": "chr20:99", "min_mapq": 0}), &pane());
        assert_eq!(a["locus"], "chr20:99");
        assert_eq!(a["min_mapq"], 0);
        assert!(!applied.contains_key("locus") && !applied.contains_key("filters"));
        let (_, applied, _) = fill("evidence", &json!({"bams": ["/x.bam"]}), &pane());
        assert!(!applied.contains_key("bams"));
    }

    #[test]
    fn svs_and_repeats_inherit_filters_and_region() {
        let mut p = pane();
        p["filters"][0]["expression"] = json!("mapq > 60");
        p["filters"][0]["min_mapq"] = json!(5);
        p["filters"][0]["exclude_flags"] = json!("0x904");
        for tool in ["svs", "repeats"] {
            // One file: its own expression, MAPQ and flags.
            let (a, applied, _) = fill(tool, &json!({"bams": ["/d/a.bam"]}), &p);
            assert_eq!(a["filter"], "mapq > 60", "{tool}");
            assert_eq!(a["min_mapq"], 5, "{tool}");
            assert_eq!(a["exclude_flags"], "0x904", "{tool}");
            assert_eq!(applied["filters"]["expression"], "mapq > 60", "{tool}");
            assert_eq!(fill(tool, &json!({"bams": ["/d/b.bam"]}), &p).0["filter"], "!supplementary", "{tool}");
            // The region from the view, or the selection.
            let (a, applied, notes) = fill(tool, &json!({}), &p);
            assert_eq!(a["region"], "chr20:1,000-2,000", "{tool}");
            assert_eq!(applied["region"]["from"], "view", "{tool}");
            assert!(notes.iter().any(|n| n.contains("different read filters")), "{tool}: {notes:?}");
            // Given filters win.
            let (a, applied, _) = fill(tool, &json!({"bams": ["/d/a.bam"], "min_mapq": 0}), &p);
            assert_eq!(a["min_mapq"], 0, "{tool}");
            assert!(a.get("filter").is_none() && !applied.contains_key("filters"), "{tool}");
        }
        let mut s = p.clone();
        s["selection"] = json!("chr20:1200-1300");
        assert_eq!(fill("repeats", &json!({}), &s).0["region"], "chr20:1200-1300");
    }

    #[test]
    fn different_filters_are_reported() {
        let mut p = pane();
        p["filters"][1]["min_mapq"] = json!(0);
        let (_, _, notes) = fill("coverage", &json!({}), &p);
        assert!(notes.iter().any(|n| n.contains("different read filters")), "{notes:?}");
        // One file asked for: its own filter.
        let (a, _, notes) = fill("coverage", &json!({"bams": ["/d/b.bam"]}), &p);
        assert_eq!(a["min_mapq"], 0);
        assert!(notes.is_empty(), "{notes:?}");
    }
}
