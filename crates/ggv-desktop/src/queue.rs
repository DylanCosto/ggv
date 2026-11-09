//! Review queues: step through any list of places (BED, VCF, a list of loci, or pasted text)
//! and mark each pass, fail or unsure with a note. Decisions are saved on every change, next to
//! the list (`<list>.ggv-queue.tsv`) or, for pasted lists, in GGV's config folder, so a review
//! can be stopped and resumed. Results export as TSV or BED.

use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

/// The places a queue steps through.
pub use ggv_io::places::{Place as Item, parse, read_file};
use ggv_io::places::{Format, parse_as};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    Pass,
    Fail,
    Unsure,
}

impl Verdict {
    pub const ALL: [Self; 3] = [Self::Pass, Self::Fail, Self::Unsure];

    pub fn name(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Fail => "fail",
            Self::Unsure => "unsure",
        }
    }

    fn parse(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|v| v.name() == s)
    }

    pub fn rgb(self) -> [u8; 3] {
        match self {
            Self::Pass => [30, 160, 60],
            Self::Fail => [210, 50, 50],
            Self::Unsure => [230, 170, 20],
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decision {
    pub verdict: Verdict,
    pub note: String,
}

pub struct Queue {
    /// What the list came from, for the window title.
    pub source: String,
    /// The places, each in its reading for the window's contigs (see `settle`).
    pub items: Vec<Item>,
    /// The places as read from the list, with both readings of an ambiguous row.
    listed: Vec<Item>,
    /// The contig names `items` were settled with.
    settled_for: Vec<String>,
    /// Some place has two readings, so the contigs matter.
    ambiguous: bool,
    /// Each item's identity (see `ids`), in step with `items`.
    ids: Vec<String>,
    /// The item shown.
    pub pos: usize,
    /// By item identity.
    pub decisions: BTreeMap<String, Decision>,
    /// Items decided here not saved yet (the save failed): kept over what is on disk.
    unsaved: BTreeSet<String>,
    /// The saves of `save_path` in this process the decisions include (see `refresh`).
    seen: u64,
    pub save_path: Option<PathBuf>,
    /// The note being typed for the current item.
    pub note: String,
    pub error: Option<String>,
    /// Where the decisions went when they couldn't be saved where they were.
    pub notice: Option<String>,
    /// The list file, when the queue came from one.
    pub list: Option<PathBuf>,
    /// Lines that weren't places (headers, junk).
    pub skipped: usize,
}

const HEADER: &str = "#locus\talleles\tentry\tlabel\tdecision\tnote";
/// Saves from before entries had identities: one decision per locus.
const OLD_HEADER: &str = "#locus\tlabel\tdecision\tnote";
const EXT: &str = "ggv-queue.tsv";

/// Each item's identity: its place, its alleles (VCF) and which copy of that pair it is, so
/// two variants at one site, or a place listed twice, keep their own decisions.
fn ids(items: &[Item]) -> Vec<String> {
    let mut seen: BTreeMap<(String, &str), usize> = BTreeMap::new();
    items
        .iter()
        .map(|it| {
            let alleles = it.alleles.as_deref().unwrap_or("");
            let n = seen.entry((it.key(), alleles)).or_default();
            *n += 1;
            format!("{}\t{}\t{n}", it.key(), clean(alleles))
        })
        .collect()
}

/// Where decisions for a list are kept.
pub fn save_path_for(list: Option<&Path>, text: &str) -> Option<PathBuf> {
    match list {
        Some(p) => Some(crate::store::sidecar(p, EXT)),
        None => {
            // FNV-1a of the pasted text: the same paste resumes the same review.
            let mut h: u64 = 0xcbf2_9ce4_8422_2325;
            for b in text.bytes() {
                h = (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
            }
            crate::session::last_session_path().map(|p| p.with_file_name("queues").join(format!("pasted-{h:016x}.tsv")))
        }
    }
}

/// Each place in its reading for a file with `contigs` ([`Item::settle`]); with no contigs,
/// as listed.
fn settled(listed: &[Item], contigs: &[ggv_model::Contig]) -> Vec<Item> {
    listed.iter().map(|p| if contigs.is_empty() { p.clone() } else { p.settle(contigs).clone() }).collect()
}

impl Queue {
    /// A queue over `text`, resuming any decisions saved for it.
    #[cfg(test)]
    pub fn new(source: String, text: &str, list: Option<&Path>) -> Self {
        Self::for_contigs(source, text, list, &[])
    }

    /// A queue over `text` shown in a window with `contigs` (each place read as it is meant
    /// there, see [`Item::settle`]), resuming any decisions saved for it.
    pub fn for_contigs(source: String, text: &str, list: Option<&Path>, contigs: &[ggv_model::Contig]) -> Self {
        let (listed, skipped) = parse_as(text, list.map_or(Format::Auto, Format::of_path));
        let items = settled(&listed, contigs);
        let ids = ids(&items);
        let save_path = save_path_for(list, text);
        let (mut notice, mut error) = (None, None);
        let mut decisions = BTreeMap::new();
        if let Some(p) = save_path.as_ref().filter(|p| p.exists()) {
            match load(p, &items, &ids) {
                Ok((d, said, _)) => {
                    decisions = d;
                    notice = said;
                }
                Err(e) => error = Some(e),
            }
        }
        let list = list.map(Path::to_path_buf);
        let seen = save_path.as_deref().map_or(0, crate::store::saves);
        let mut q = Self {
            source,
            items,
            settled_for: contigs.iter().map(|c| c.name.clone()).collect(),
            ambiguous: listed.iter().any(|p| p.other.is_some()),
            listed,
            ids,
            pos: 0,
            decisions,
            unsaved: BTreeSet::new(),
            seen,
            save_path,
            note: String::new(),
            error,
            notice,
            list,
            skipped,
        };
        // Resume at the first undecided item, with its saved note in the editor: a completed
        // queue reopens on item 1, and re-marking it must not replace its note with "".
        q.pos = q.next_undecided(0).unwrap_or(0);
        q.load_note();
        q
    }

    /// Reads each place as it is meant in a window with `contigs` (a BED row on a contig with
    /// `:` in its name, or a locus with a label: see [`Item::settle`]), when they changed. A
    /// decision stays with its entry.
    pub fn settle(&mut self, contigs: &[ggv_model::Contig]) {
        if !self.ambiguous || self.settled_for.iter().map(String::as_str).eq(contigs.iter().map(|c| c.name.as_str())) {
            return;
        }
        self.settled_for = contigs.iter().map(|c| c.name.clone()).collect();
        self.items = settled(&self.listed, contigs);
        let new = ids(&self.items);
        if new == self.ids {
            return;
        }
        let mut decisions = BTreeMap::new();
        let mut unsaved = BTreeSet::new();
        for (old, id) in self.ids.iter().zip(&new) {
            if let Some(d) = self.decisions.remove(old) {
                decisions.insert(id.clone(), d);
            }
            if self.unsaved.contains(old) {
                unsaved.insert(id.clone());
            }
        }
        self.decisions = decisions;
        self.unsaved = unsaved;
        self.ids = new;
        // Decisions saved under this reading (by an earlier session, or another pane).
        if let Some(path) = self.save_path.clone() {
            self.take_in(&path, false);
        }
        self.load_note();
    }

    pub fn current(&self) -> Option<&Item> {
        self.items.get(self.pos)
    }

    pub fn decision(&self, i: usize) -> Option<&Decision> {
        self.ids.get(i).and_then(|id| self.decisions.get(id))
    }

    pub fn counts(&self) -> [usize; 3] {
        let mut c = [0; 3];
        for i in 0..self.items.len() {
            if let Some(d) = self.decision(i) {
                c[d.verdict as usize] += 1;
            }
        }
        c
    }

    /// The first undecided item at or after `from`, wrapping around.
    pub fn next_undecided(&self, from: usize) -> Option<usize> {
        let n = self.items.len();
        (0..n).map(|k| (from + k) % n).find(|&i| self.decision(i).is_none())
    }

    pub fn step(&mut self, forward: bool) {
        if self.items.is_empty() {
            return;
        }
        self.pos = if forward { (self.pos + 1).min(self.items.len() - 1) } else { self.pos.saturating_sub(1) };
        self.load_note();
    }

    pub fn go(&mut self, i: usize) {
        if i < self.items.len() {
            self.pos = i;
            self.load_note();
        }
    }

    fn load_note(&mut self) {
        self.note = self.decision(self.pos).map(|d| d.note.clone()).unwrap_or_default();
    }

    /// Records a verdict for the current item and moves to the next undecided one.
    pub fn decide(&mut self, verdict: Verdict) {
        let Some(key) = self.ids.get(self.pos).cloned() else { return };
        let note = std::mem::take(&mut self.note);
        self.decisions.insert(key.clone(), Decision { verdict, note });
        self.unsaved.insert(key);
        // Another pane with this list may have saved since: the save takes its decisions in,
        // and this one, the latest, goes over them.
        self.relocate();
        self.save();
        if let Some(next) = self.next_undecided(self.pos + 1) {
            self.pos = next;
        }
        self.load_note();
    }

    /// Takes in decisions another pane with the same list saved since this one last read or
    /// saved them (also when its saves had to move elsewhere).
    pub fn refresh(&mut self) {
        if self.relocate() || self.save_path.as_deref().is_some_and(|p| crate::store::saves(p) != self.seen) {
            self.merge_saved();
        }
    }

    /// Follows the list's decisions to where they were last saved (by another pane, after a
    /// save had to go elsewhere); whether they moved.
    fn relocate(&mut self) -> bool {
        let Some(list) = &self.list else { return false };
        match crate::store::located(list, EXT) {
            Some(p) if self.save_path.as_ref() != Some(&p) => {
                self.save_path = Some(p);
                true
            }
            _ => false,
        }
    }

    /// The saved decisions taken in (see `take_in`).
    fn merge_saved(&mut self) {
        let Some(path) = self.save_path.clone() else { return };
        self.seen = crate::store::saves(&path);
        self.take_in(&path, false);
    }

    /// Takes in the decisions saved at `path` ([`crate::store::merge_decisions`]; a file
    /// that's gone or can't be read keeps what this queue has). What it has wins over a
    /// damaged file's older backup, and with `moving` (`path` isn't where they were read).
    /// Whether the file is damaged.
    fn take_in(&mut self, path: &Path, moving: bool) -> bool {
        if !path.exists() {
            return false;
        }
        match load(path, &self.items, &self.ids) {
            Ok((saved, _, recovered)) => {
                self.decisions = crate::store::merge_decisions(&self.decisions, &self.unsaved, saved, recovered || moving);
                recovered
            }
            Err(_) => false,
        }
    }

    fn save(&mut self) {
        let Some(mut path) = self.save_path.take() else { return };
        let from = path.clone();
        let list = self.list.clone();
        let mut contents = |p: &Path| {
            let damaged = self.take_in(p, p != from);
            (self.to_tsv(false).into_bytes(), damaged)
        };
        let result = match list {
            Some(list) => crate::store::save_sidecar(&mut path, &list, EXT, contents),
            None => crate::store::save_merged(&path, &mut contents).map(|()| None).map_err(|e| format!("{}: {e}", path.display())),
        };
        match result {
            Ok(moved) => {
                self.error = None;
                self.unsaved.clear();
                self.seen = crate::store::saved(&path);
                if moved.is_some() {
                    self.notice = moved;
                }
            }
            Err(e) => self.error = Some(e),
        }
        self.save_path = Some(path);
    }

    /// TSV of items with their decisions (`all`: undecided ones too).
    pub fn to_tsv(&self, all: bool) -> String {
        let mut out = String::from(HEADER);
        out.push('\n');
        for (i, (it, id)) in self.items.iter().zip(&self.ids).enumerate() {
            let d = self.decision(i);
            if d.is_none() && !all {
                continue;
            }
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                id,
                clean(&it.label),
                d.map_or("", |d| d.verdict.name()),
                d.map_or(String::new(), |d| clean(&d.note))
            ));
        }
        out
    }

    /// BED of decided items (name: `verdict` or `verdict:note`). Whole-contig items need the
    /// contig length, so they are left out.
    pub fn to_bed(&self) -> String {
        let mut out = String::new();
        for (i, it) in self.items.iter().enumerate() {
            let (Some(d), Some(end)) = (self.decision(i), it.end) else { continue };
            let name = if d.note.is_empty() { d.verdict.name().to_owned() } else { format!("{}:{}", d.verdict.name(), clean(&d.note)) };
            // An insertion site (a zero-length row, shown as the base before it) goes back out
            // as it was listed, so the export re-imports and joins on the same coordinates.
            let start = if it.insertion { end } else { it.start };
            out.push_str(&format!("{}\t{}\t{}\t{}\n", it.contig, start, end, name.replace(' ', "_")));
        }
        out
    }
}

fn clean(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

/// Decisions saved for a queue: by item identity, what to tell the person about them, and
/// whether they are the backup of a damaged file.
type Loaded = (BTreeMap<String, Decision>, Option<String>, bool);

/// The decisions saved at `path` for `items`, and what to tell the person about them (the file
/// was damaged, or an older save's decisions were shared), and whether the file was damaged.
fn load(path: &Path, items: &[Item], ids: &[String]) -> Result<Loaded, String> {
    let (saved, recovered) = crate::store::read_recovering(path, parse_decisions).map_err(|e| format!("{} can't be read, nor its backup: {e}", path.display()))?;
    let mut notice = recovered.then(|| format!("{} was damaged; restored the previous save", path.display()));
    let decisions = match saved {
        Saved::ById(d) => d,
        Saved::ByLocus(d) => {
            let (d, shared) = from_loci(items, ids, &d);
            if shared > 0 {
                notice = Some(format!("{shared} decision(s) in an older save were for a place listed more than once; each entry there shows it, so check them"));
            }
            d
        }
    };
    Ok((decisions, notice, recovered))
}

/// Saved decisions: by item identity, or (older saves) by locus.
enum Saved {
    ById(BTreeMap<String, Decision>),
    ByLocus(BTreeMap<String, Decision>),
}

/// A saved queue: empty, or a header then one decision per line (blank lines and `#` comments
/// are fine). Anything else is damage, so the backup is tried.
fn parse_decisions(text: &str) -> anyhow::Result<Saved> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(first) = lines.next() else { return Ok(Saved::ById(BTreeMap::new())) };
    // Columns before the decision: locus, alleles, entry, label (older saves: locus, label).
    let before = match first.trim_end() {
        HEADER => 4,
        OLD_HEADER => 2,
        _ => anyhow::bail!("not a saved review queue"),
    };
    let mut out = BTreeMap::new();
    for l in lines.filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.trim_end_matches('\r').split('\t').collect();
        let ok = (f.len() == before + 1 || f.len() == before + 2) && !f[0].is_empty() && (before == 2 || f[2].parse::<usize>().is_ok_and(|n| n > 0));
        let verdict = f.get(before).and_then(|v| Verdict::parse(v)).filter(|_| ok).ok_or_else(|| anyhow::anyhow!("damaged line: {l}"))?;
        let key = if before == 4 { f[..3].join("\t") } else { f[0].to_owned() };
        out.insert(key, Decision { verdict, note: f.get(before + 1).unwrap_or(&"").to_string() });
    }
    Ok(if before == 4 { Saved::ById(out) } else { Saved::ByLocus(out) })
}

/// Decisions saved by locus, by item identity, and how many were for a locus listed more than
/// once (every entry there shared it then, so each keeps it, to be checked).
fn from_loci(items: &[Item], ids: &[String], saved: &BTreeMap<String, Decision>) -> (BTreeMap<String, Decision>, usize) {
    let mut out = BTreeMap::new();
    let mut shared = std::collections::BTreeSet::new();
    for (it, id) in items.iter().zip(ids) {
        let key = it.key();
        if let Some(d) = saved.get(&key) {
            out.insert(id.clone(), d.clone());
            if items.iter().filter(|x| x.key() == key).count() > 1 {
                shared.insert(key);
            }
        }
    }
    (out, shared.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decisions_resume_and_export() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("sites.bed");
        let text = "chr1\t0\t10\ta\nchr1\t20\t30\tb\nchr1\t40\t50\tc\n";
        let mut q = Queue::new("sites.bed".into(), text, Some(&list));
        assert_eq!(q.pos, 0);
        q.note = "clean".into();
        q.decide(Verdict::Pass);
        assert_eq!(q.pos, 1, "moves to the next undecided item");
        q.decide(Verdict::Fail);
        assert_eq!(q.counts(), [1, 1, 0]);
        assert_eq!(q.to_bed(), "chr1\t0\t10\tpass:clean\nchr1\t20\t30\tfail\n");
        assert!(q.to_tsv(true).ends_with("chr1:41-50\t\t1\tc\t\t\n"));

        let again = Queue::new("sites.bed".into(), text, Some(&list));
        assert_eq!(again.pos, 2, "resumes at the first undecided item");
        assert_eq!(again.decision(0).unwrap().note, "clean");
        assert_eq!(again.save_path.as_deref(), Some(dir.join("sites.bed.ggv-queue.tsv").as_path()));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn reconfirming_a_completed_queue_keeps_its_note() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-done-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("done.bed");
        let text = "chr1\t0\t10\ta\nchr1\t20\t30\tb\n";
        let mut q = Queue::new("done.bed".into(), text, Some(&list));
        q.note = "KEEP_THIS_NOTE".into();
        q.decide(Verdict::Pass);
        q.decide(Verdict::Fail);

        // Reopened complete: it starts on item 1 with that item's note ready to edit.
        let mut again = Queue::new("done.bed".into(), text, Some(&list));
        assert_eq!(again.pos, 0);
        assert_eq!(again.note, "KEEP_THIS_NOTE");
        again.decide(Verdict::Pass);
        let reread = Queue::new("done.bed".into(), text, Some(&list));
        assert_eq!(reread.decision(0).unwrap().note, "KEEP_THIS_NOTE");
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn entries_at_one_locus_keep_their_own_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-ids-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("queue.vcf");
        // Two variants at one site, and a place listed twice.
        let text = "chrQA\t110\tfirst\tA\tC\t.\t.\t.\nchrQA\t110\tsecond\tA\tG\t.\t.\t.\nchrQA:200\tx\nchrQA:200\ty\n";
        let mut q = Queue::new("queue.vcf".into(), text, Some(&list));
        assert_eq!(q.items.len(), 4);
        q.note = "first only".into();
        q.decide(Verdict::Pass);
        assert_eq!((q.decision(0).map(|d| d.verdict), q.decision(1)), (Some(Verdict::Pass), None));
        assert_eq!(q.counts(), [1, 0, 0]);
        q.go(1);
        q.note = "second only".into();
        q.decide(Verdict::Fail);
        q.go(3);
        q.decide(Verdict::Unsure);
        assert_eq!(q.counts(), [1, 1, 1]);
        assert!(q.decision(2).is_none(), "the first copy of a repeated place is still undecided");

        let again = Queue::new("queue.vcf".into(), text, Some(&list));
        assert_eq!(again.decision(0).unwrap(), &Decision { verdict: Verdict::Pass, note: "first only".into() });
        assert_eq!(again.decision(1).unwrap(), &Decision { verdict: Verdict::Fail, note: "second only".into() });
        assert!(again.decision(2).is_none());
        assert_eq!(again.decision(3).unwrap().verdict, Verdict::Unsure);
        assert_eq!(again.pos, 2);
        let tsv = again.to_tsv(true);
        assert!(tsv.starts_with(HEADER));
        assert!(tsv.contains("chrQA:110\tA>C\t1\tfirst\tpass\tfirst only\n"), "{tsv}");
        assert!(tsv.contains("chrQA:110\tA>G\t1\tsecond\tfail\tsecond only\n"), "{tsv}");
        assert!(tsv.contains("chrQA:200\t\t2\ty\tunsure\t\n"), "{tsv}");
        assert_eq!(again.to_bed(), "chrQA\t109\t110\tpass:first_only\nchrQA\t109\t110\tfail:second_only\nchrQA\t199\t200\tunsure\n");
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The same list open in two panes: each pane's save keeps the other's decisions, and
    /// each sees the other's once refreshed.
    #[test]
    fn two_queues_on_one_list_keep_each_others_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-two-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("sites.bed");
        let text = "chr1	0	10	a
chr1	20	30	b
chr1	40	50	c
";
        let mut one = Queue::new("sites.bed".into(), text, Some(&list));
        let mut two = Queue::new("sites.bed".into(), text, Some(&list));
        one.decide(Verdict::Pass);
        two.go(1);
        two.decide(Verdict::Fail);
        assert_eq!(two.counts(), [1, 1, 0], "the other pane's decision is kept");
        one.refresh();
        assert_eq!(one.counts(), [1, 1, 0], "and seen there");
        assert_eq!(one.pos, 1);
        one.go(2);
        one.decide(Verdict::Unsure);
        // The latest decision for an item wins.
        two.go(0);
        two.decide(Verdict::Fail);
        let again = Queue::new("sites.bed".into(), text, Some(&list));
        assert_eq!(again.counts(), [0, 2, 1]);
        assert_eq!(again.decision(0).unwrap().verdict, Verdict::Fail);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A list in a folder that can't be written (a file), open in two panes before any
    /// decision: both follow the decisions to GGV's own folder, keeping each other's.
    #[test]
    fn two_queues_on_a_list_in_an_unwritable_folder_keep_each_others_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-readonly-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        crate::store::TEST_CONFIG.set(Some(dir.join("config")));
        std::fs::write(dir.join("share"), "not a folder").unwrap();
        let list = dir.join("share").join("sites.bed");
        let text = "chr1\t0\t10\ta\nchr1\t20\t30\tb\nchr1\t40\t50\tc\n";
        let mut one = Queue::new("sites.bed".into(), text, Some(&list));
        let mut two = Queue::new("sites.bed".into(), text, Some(&list));
        one.decide(Verdict::Pass);
        one.decide(Verdict::Pass);
        assert!(one.save_path.as_ref().is_some_and(|p| p.starts_with(dir.join("config"))), "{:?}", one.save_path);
        two.refresh();
        assert_eq!((two.counts(), &two.save_path), ([2, 0, 0], &one.save_path));
        two.go(2);
        two.decide(Verdict::Fail);
        one.refresh();
        assert_eq!(one.counts(), [2, 1, 0]);
        assert_eq!(Queue::new("sites.bed".into(), text, Some(&list)).counts(), [2, 1, 0]);
        crate::store::TEST_CONFIG.set(None);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The decisions file damaged while reviewing: the next decision keeps the one before it,
    /// and the good backup stays.
    #[test]
    fn a_damaged_save_keeps_the_latest_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-damaged-live-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("sites.bed");
        let text = "chr1\t0\t10\ta\nchr1\t20\t30\tb\nchr1\t40\t50\tc\n";
        let mut q = Queue::new("sites.bed".into(), text, Some(&list));
        q.decide(Verdict::Pass);
        q.decide(Verdict::Pass);
        // Changing a decision: the backup has the old one.
        q.go(1);
        q.note = "important".into();
        q.decide(Verdict::Fail);
        let saved = q.save_path.clone().unwrap();
        std::fs::write(&saved, "garbage\n").unwrap();
        q.go(2);
        q.decide(Verdict::Unsure);
        assert_eq!(q.counts(), [1, 1, 1]);
        let again = Queue::new("sites.bed".into(), text, Some(&list));
        assert_eq!((again.counts(), again.notice.clone()), ([1, 1, 1], None));
        assert_eq!(again.decision(1).unwrap().note, "important");
        let backup = std::fs::read_to_string(crate::store::backup_of(&saved)).unwrap();
        assert!(!backup.contains("garbage") && backup.contains("\tpass\t"), "{backup}");
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A BED row on a contig with `:` in its name (a `samtools faidx` sub-reference) in a .txt
    /// list reads as a locus with a label until the window's contigs say otherwise; decisions
    /// stay with their entries, and BED export writes the rows.
    #[test]
    fn rows_on_contigs_with_colons_settle_with_the_genome() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-settle-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("sub.txt");
        let text = "chr20:10000000-20000000\t500\t501\tsiteA\nchr20:10000000-20000000\t1000\t1000\tins\n";
        let sub = [ggv_model::Contig { name: "chr20:10000000-20000000".into(), length: 10_000_001 }];
        let mut q = Queue::new("sub.txt".into(), text, Some(&list));
        assert_eq!(q.items[0].contig, "chr20", "no genome yet: the primary reading");
        q.decide(Verdict::Pass);
        q.settle(&sub);
        assert_eq!((q.items[0].key(), q.items[0].label.as_str()), ("chr20:10000000-20000000:501".to_owned(), "siteA"));
        assert_eq!(q.decision(0).map(|d| d.verdict), Some(Verdict::Pass), "the decision stays with its entry");
        q.go(1);
        q.decide(Verdict::Fail);
        assert_eq!(q.to_bed(), "chr20:10000000-20000000\t500\t501\tpass\nchr20:10000000-20000000\t1000\t1000\tfail\n");
        // Reopened in a window with that genome: read as rows, with their decisions.
        let again = Queue::for_contigs("sub.txt".into(), text, Some(&list), &sub);
        assert_eq!(again.counts(), [1, 1, 0]);
        assert_eq!(again.items[1].key(), "chr20:10000000-20000000:1000");
        // Another genome (chr20 itself): the locus reading again, decisions kept.
        q.settle(&[ggv_model::Contig { name: "chr20".into(), length: 64_444_167 }]);
        assert_eq!((q.items[0].contig.as_str(), q.counts()), ("chr20", [1, 1, 0]));
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A zero-length BED row (an insertion site) exports as it was listed.
    #[test]
    fn insertion_sites_export_as_zero_length_bed() {
        let mut q = Queue::new("pasted list".into(), "chr20	105	105	ins-site
chr20	200	210	region
", None);
        q.save_path = None;
        q.decide(Verdict::Pass);
        q.decide(Verdict::Fail);
        assert_eq!(q.to_bed(), "chr20	105	105	pass
chr20	200	210	fail
");
        let (back, _) = parse(&q.to_bed());
        assert_eq!((back[0].start, back[0].end, back[0].insertion), (q.items[0].start, q.items[0].end, true), "re-imports as the same site");
    }

    #[test]
    fn old_locus_keyed_saves_still_load() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-old-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("sites.bed");
        let text = "chr1\t0\t10\ta\nchr1\t20\t30\tb\nchr1\t20\t30\tc\n";
        std::fs::write(dir.join("sites.bed.ggv-queue.tsv"), "#locus\tlabel\tdecision\tnote\nchr1:1-10\ta\tpass\tok\nchr1:21-30\tb\tfail\t\n").unwrap();
        let q = Queue::new("sites.bed".into(), text, Some(&list));
        assert_eq!(q.error, None);
        assert_eq!(q.decision(0).unwrap(), &Decision { verdict: Verdict::Pass, note: "ok".into() });
        // The old decision covered both entries at chr1:21-30: each keeps it, and it's said.
        assert_eq!(q.counts(), [1, 2, 0]);
        assert!(q.notice.as_deref().is_some_and(|n| n.contains("1 decision(s)")), "{:?}", q.notice);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn damaged_saves_fall_back_to_the_backup() {
        let dir = std::env::temp_dir().join(format!("ggv-queue-damaged-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("qa.txt");
        let text = "chrQA:110\nchrQA:150\n";
        let mut q = Queue::new("qa.txt".into(), text, Some(&list));
        q.note = "first".into();
        q.decide(Verdict::Pass);
        q.note = "second".into();
        q.decide(Verdict::Fail);
        let saved = q.save_path.clone().unwrap();
        let bak = crate::store::backup_of(&saved);
        assert!(std::fs::read_to_string(&bak).unwrap().contains("pass\tfirst"));
        let good = std::fs::read_to_string(&saved).unwrap();
        // Garbage, a truncated header, truncated or bad records: the backup is used and said.
        let header = good.lines().next().unwrap().to_owned();
        for damaged in [
            "not a valid saved queue".to_owned(),
            header[..10].to_owned(),
            format!("{header}\nchrQA:110\t\t1\t\tpa"),
            format!("{header}\nchrQA:110\t\t1"),
            format!("{header}\nchrQA:110\t\tone\t\tpass\t"),
            "#locus\tlabel\tdecision\tnote\nchrQA:110\t\tmaybe\t\n".to_owned(),
        ] {
            std::fs::write(&saved, &damaged).unwrap();
            let q = Queue::new("qa.txt".into(), text, Some(&list));
            assert_eq!(q.decision(0).unwrap(), &Decision { verdict: Verdict::Pass, note: "first".into() }, "{damaged:?}");
            assert!(q.notice.as_deref().is_some_and(|n| n.contains("damaged")), "{damaged:?}");
            assert_eq!(q.error, None);
        }
        // Blank lines, comments and CRLF are fine; so is an empty file.
        std::fs::write(&saved, format!("\n{}\n# a comment\n\n", good.trim_end().replace('\n', "\r\n"))).unwrap();
        let q = Queue::new("qa.txt".into(), text, Some(&list));
        assert_eq!((q.counts(), q.notice.clone()), ([1, 1, 0], None));
        std::fs::write(&saved, "").unwrap();
        let q = Queue::new("qa.txt".into(), text, Some(&list));
        assert_eq!((q.counts(), q.notice.clone(), q.error.clone()), ([0, 0, 0], None, None));
        // Neither readable: an error, not a fresh review.
        std::fs::write(&saved, "junk").unwrap();
        std::fs::write(&bak, "junk").unwrap();
        let q = Queue::new("qa.txt".into(), text, Some(&list));
        assert!(q.error.as_deref().is_some_and(|e| e.contains("can't be read")), "{:?}", q.error);
        let _ = std::fs::remove_dir_all(dir);
    }
}
