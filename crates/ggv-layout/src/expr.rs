//! Read filter expressions, e.g. `mapq >= 20 && hp == 1 && !dup`.
//!
//! An expression is parsed once into a small tree and evaluated per read with no allocation.
//! Compiled expressions are interned in a process-wide table and referred to by [`ExprId`], so
//! the display style that carries one stays `Copy` and cheap to compare.
//!
//! Values are numbers, strings, booleans, or missing (`hp` on a read with no `HP` tag). A
//! missing value equals nothing: `hp == 1` is false and `hp != 1` is true for untagged reads.

use std::sync::{Arc, OnceLock, RwLock};

use crate::AlignedRead;

/// Fields a filter can test, with what they mean (shown in the filter window).
pub const FIELDS: &[(&str, &str)] = &[
    ("mapq", "mapping quality"),
    ("flag", "SAM flags (use & to test bits: flag & 0x400)"),
    ("start, end", "1-based first and last aligned base"),
    ("len", "aligned length on the reference"),
    ("tlen, insert", "template length, signed / absolute"),
    ("strand", "\"+\" or \"-\"  (also: reverse, forward)"),
    ("dup, secondary, supplementary, qcfail", "flag tests"),
    ("paired, proper, read1, read2", "pairing flag tests"),
    ("mate_unmapped, mate_reverse, mate_elsewhere", "mate tests (elsewhere: other contig)"),
    ("hp, rg  (or tag(HP), tag(RG))", "haplotype and read group tags"),
    ("name", "read name; ~ tests containment: name ~ \"ccs\""),
    ("mismatches, deletions, insertions", "event counts (mismatches need a reference)"),
    ("clips, clip_left, clip_right", "soft-clipped bases"),
    ("mods, split", "base-modification calls; has an SA tag"),
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Field {
    Mapq,
    Flag,
    Start,
    End,
    Len,
    Tlen,
    Insert,
    Strand,
    Reverse,
    Forward,
    Dup,
    Secondary,
    Supplementary,
    QcFail,
    Paired,
    Proper,
    Read1,
    Read2,
    MateUnmapped,
    MateReverse,
    MateElsewhere,
    Hp,
    Rg,
    Name,
    Mismatches,
    Deletions,
    Insertions,
    Clips,
    ClipLeft,
    ClipRight,
    Mods,
    Split,
}

fn field(name: &str) -> Option<Field> {
    use Field::*;
    Some(match name {
        "mapq" | "mq" => Mapq,
        "flag" | "flags" => Flag,
        "start" | "pos" => Start,
        "end" => End,
        "len" | "length" => Len,
        "tlen" => Tlen,
        "insert" | "isize" => Insert,
        "strand" => Strand,
        "reverse" => Reverse,
        "forward" => Forward,
        "dup" | "duplicate" => Dup,
        "secondary" => Secondary,
        "supplementary" | "supp" => Supplementary,
        "qcfail" => QcFail,
        "paired" => Paired,
        "proper" => Proper,
        "read1" | "first" => Read1,
        "read2" | "second" => Read2,
        "mate_unmapped" => MateUnmapped,
        "mate_reverse" => MateReverse,
        "mate_elsewhere" => MateElsewhere,
        "hp" => Hp,
        "rg" => Rg,
        "name" => Name,
        "mismatches" | "mm" => Mismatches,
        "deletions" | "dels" => Deletions,
        "insertions" | "ins" => Insertions,
        "clips" => Clips,
        "clip_left" => ClipLeft,
        "clip_right" => ClipRight,
        "mods" => Mods,
        "split" => Split,
        _ => return None,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    And,
    Or,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Contains,
    BitAnd,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Expr(Node);

#[derive(Clone, Debug, PartialEq)]
enum Node {
    Num(f64),
    Str(String),
    Bool(bool),
    Field(Field),
    Not(Box<Node>),
    Bin(Op, Box<Node>, Box<Node>),
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Val<'a> {
    Num(f64),
    Str(&'a str),
    Bool(bool),
    Missing,
}

impl Val<'_> {
    fn truthy(self) -> bool {
        match self {
            Val::Num(n) => n != 0.0,
            Val::Str(s) => !s.is_empty(),
            Val::Bool(b) => b,
            Val::Missing => false,
        }
    }

    fn num(self) -> Option<f64> {
        match self {
            Val::Num(n) => Some(n),
            Val::Bool(b) => Some(f64::from(u8::from(b))),
            _ => None,
        }
    }
}

impl Expr {
    pub fn parse(text: &str) -> Result<Self, String> {
        let tokens = lex(text)?;
        if tokens.is_empty() {
            return Err("empty expression".into());
        }
        let mut p = Parser { tokens, i: 0 };
        let node = p.or()?;
        if p.i < p.tokens.len() {
            return Err(format!("unexpected {} at {}", p.tokens[p.i].0.describe(), p.tokens[p.i].1 + 1));
        }
        Ok(Self(node))
    }

    /// Whether `r` passes.
    pub fn matches(&self, r: &AlignedRead) -> bool {
        eval(&self.0, r).truthy()
    }

    /// Whether the expression tests `mismatches`, which only a reference can give.
    pub fn uses_mismatches(&self) -> bool {
        fn walk(n: &Node) -> bool {
            match n {
                Node::Field(f) => *f == Field::Mismatches,
                Node::Not(x) => walk(x),
                Node::Bin(_, a, b) => walk(a) || walk(b),
                Node::Num(_) | Node::Str(_) | Node::Bool(_) => false,
            }
        }
        walk(&self.0)
    }
}

fn read_field(f: Field, r: &AlignedRead) -> Val<'_> {
    let flag = |bit: u16| Val::Bool(r.flags & bit != 0);
    let n = |x: usize| Val::Num(x as f64);
    match f {
        Field::Mapq => Val::Num(f64::from(r.mapq)),
        Field::Flag => Val::Num(f64::from(r.flags)),
        Field::Start => Val::Num((r.start + 1) as f64),
        Field::End => Val::Num(r.end as f64),
        Field::Len => Val::Num((r.end - r.start) as f64),
        Field::Tlen => Val::Num(f64::from(r.template_length)),
        Field::Insert => Val::Num(f64::from(r.template_length.unsigned_abs())),
        Field::Strand => Val::Str(if r.reverse { "-" } else { "+" }),
        Field::Reverse => Val::Bool(r.reverse),
        Field::Forward => Val::Bool(!r.reverse),
        Field::Dup => flag(0x400),
        Field::Secondary => flag(0x100),
        Field::Supplementary => flag(0x800),
        Field::QcFail => flag(0x200),
        Field::Paired => flag(0x1),
        Field::Proper => flag(0x2),
        Field::Read1 => flag(0x40),
        Field::Read2 => flag(0x80),
        Field::MateUnmapped => flag(0x8),
        Field::MateReverse => flag(0x20),
        Field::MateElsewhere => Val::Bool(r.mate.is_some_and(|m| !m.same_contig)),
        Field::Hp => r.haplotype.map_or(Val::Missing, |h| Val::Num(f64::from(h))),
        Field::Rg => r.read_group.as_deref().map_or(Val::Missing, Val::Str),
        Field::Name => Val::Str(&r.name),
        Field::Mismatches => n(r.mismatches.len()),
        Field::Deletions => n(r.deletions.len()),
        Field::Insertions => n(r.insertions.len()),
        Field::Clips => Val::Num(f64::from(r.soft_clip_left + r.soft_clip_right)),
        Field::ClipLeft => Val::Num(f64::from(r.soft_clip_left)),
        Field::ClipRight => Val::Num(f64::from(r.soft_clip_right)),
        Field::Mods => n(r.mods.len()),
        Field::Split => Val::Bool(!r.supplementary.is_empty()),
    }
}

fn eval<'a>(node: &'a Node, r: &'a AlignedRead) -> Val<'a> {
    match node {
        Node::Num(n) => Val::Num(*n),
        Node::Str(s) => Val::Str(s),
        Node::Bool(b) => Val::Bool(*b),
        Node::Field(f) => read_field(*f, r),
        Node::Not(x) => Val::Bool(!eval(x, r).truthy()),
        Node::Bin(Op::And, a, b) => Val::Bool(eval(a, r).truthy() && eval(b, r).truthy()),
        Node::Bin(Op::Or, a, b) => Val::Bool(eval(a, r).truthy() || eval(b, r).truthy()),
        Node::Bin(op, a, b) => {
            let (x, y) = (eval(a, r), eval(b, r));
            match op {
                Op::Eq => Val::Bool(equal(x, y)),
                Op::Ne => Val::Bool(!equal(x, y)),
                Op::Contains => Val::Bool(matches!((x, y), (Val::Str(s), Val::Str(t)) if s.contains(t))),
                Op::BitAnd => match (x.num(), y.num()) {
                    (Some(p), Some(q)) => Val::Num(((p as u64) & (q as u64)) as f64),
                    _ => Val::Missing,
                },
                _ => {
                    let (Some(p), Some(q)) = (x.num(), y.num()) else { return Val::Bool(false) };
                    Val::Bool(match op {
                        Op::Lt => p < q,
                        Op::Le => p <= q,
                        Op::Gt => p > q,
                        _ => p >= q,
                    })
                }
            }
        }
    }
}

fn equal(x: Val, y: Val) -> bool {
    match (x, y) {
        (Val::Missing, _) | (_, Val::Missing) => false,
        (Val::Str(a), Val::Str(b)) => a == b,
        (Val::Str(_), _) | (_, Val::Str(_)) => false,
        _ => x.num() == y.num(),
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(f64),
    Str(String),
    Ident(String),
    Op(Op),
    Not,
    Minus,
    Open,
    Close,
}

impl Tok {
    fn describe(&self) -> String {
        match self {
            Tok::Num(n) => format!("number {n}"),
            Tok::Str(s) => format!("string \"{s}\""),
            Tok::Ident(s) => format!("'{s}'"),
            Tok::Op(_) | Tok::Not | Tok::Minus => "operator".into(),
            Tok::Open => "'('".into(),
            Tok::Close => "')'".into(),
        }
    }
}

fn lex(text: &str) -> Result<Vec<(Tok, usize)>, String> {
    let b = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        let start = i;
        let two = |s: &str| text[i..].starts_with(s);
        let tok = match c {
            b' ' | b'\t' | b'\n' | b'\r' => {
                i += 1;
                continue;
            }
            b'(' => {
                i += 1;
                Tok::Open
            }
            b'-' => {
                i += 1;
                Tok::Minus
            }
            b')' => {
                i += 1;
                Tok::Close
            }
            _ if two("&&") => {
                i += 2;
                Tok::Op(Op::And)
            }
            _ if two("||") => {
                i += 2;
                Tok::Op(Op::Or)
            }
            _ if two("==") => {
                i += 2;
                Tok::Op(Op::Eq)
            }
            _ if two("!=") => {
                i += 2;
                Tok::Op(Op::Ne)
            }
            _ if two("<=") => {
                i += 2;
                Tok::Op(Op::Le)
            }
            _ if two(">=") => {
                i += 2;
                Tok::Op(Op::Ge)
            }
            b'<' | b'>' | b'=' | b'~' | b'&' | b'!' => {
                i += 1;
                match c {
                    b'<' => Tok::Op(Op::Lt),
                    b'>' => Tok::Op(Op::Gt),
                    b'=' => Tok::Op(Op::Eq),
                    b'~' => Tok::Op(Op::Contains),
                    b'&' => Tok::Op(Op::BitAnd),
                    _ => Tok::Not,
                }
            }
            b'"' | b'\'' => {
                let end = text[i + 1..].find(c as char).ok_or_else(|| format!("unclosed string at {}", i + 1))?;
                let s = text[i + 1..i + 1 + end].to_owned();
                i += end + 2;
                Tok::Str(s)
            }
            b'0'..=b'9' | b'.' => {
                if two("0x") || two("0X") {
                    let mut j = i + 2;
                    while j < b.len() && b[j].is_ascii_hexdigit() {
                        j += 1;
                    }
                    let n = u64::from_str_radix(&text[i + 2..j], 16).map_err(|_| format!("bad hex number at {}", i + 1))?;
                    i = j;
                    Tok::Num(n as f64)
                } else {
                    let mut j = i;
                    while j < b.len() && (b[j].is_ascii_digit() || b[j] == b'.' || b[j] == b'_' || b[j] == b',') {
                        j += 1;
                    }
                    let digits: String = text[i..j].chars().filter(|c| *c != '_' && *c != ',').collect();
                    let n = digits.parse().map_err(|_| format!("bad number at {}", i + 1))?;
                    i = j;
                    Tok::Num(n)
                }
            }
            c if c.is_ascii_alphabetic() || c == b'_' => {
                let mut j = i;
                while j < b.len() && (b[j].is_ascii_alphanumeric() || b[j] == b'_') {
                    j += 1;
                }
                let word = text[i..j].to_ascii_lowercase();
                i = j;
                match word.as_str() {
                    "and" => Tok::Op(Op::And),
                    "or" => Tok::Op(Op::Or),
                    "not" => Tok::Not,
                    _ => Tok::Ident(word),
                }
            }
            _ => return Err(format!("unexpected '{}' at {}", text[i..].chars().next().unwrap_or('?'), i + 1)),
        };
        out.push((tok, start));
    }
    Ok(out)
}

struct Parser {
    tokens: Vec<(Tok, usize)>,
    i: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.tokens.get(self.i).map(|t| &t.0)
    }

    fn next(&mut self) -> Result<Tok, String> {
        let t = self.tokens.get(self.i).map(|t| t.0.clone()).ok_or("expression ends too soon")?;
        self.i += 1;
        Ok(t)
    }

    fn or(&mut self) -> Result<Node, String> {
        let mut left = self.and()?;
        while self.peek() == Some(&Tok::Op(Op::Or)) {
            self.i += 1;
            left = Node::Bin(Op::Or, Box::new(left), Box::new(self.and()?));
        }
        Ok(left)
    }

    fn and(&mut self) -> Result<Node, String> {
        let mut left = self.unary()?;
        while self.peek() == Some(&Tok::Op(Op::And)) {
            self.i += 1;
            left = Node::Bin(Op::And, Box::new(left), Box::new(self.unary()?));
        }
        Ok(left)
    }

    fn unary(&mut self) -> Result<Node, String> {
        if self.peek() == Some(&Tok::Not) {
            self.i += 1;
            return Ok(Node::Not(Box::new(self.unary()?)));
        }
        self.cmp()
    }

    fn cmp(&mut self) -> Result<Node, String> {
        let left = self.bit()?;
        if let Some(Tok::Op(op @ (Op::Eq | Op::Ne | Op::Lt | Op::Le | Op::Gt | Op::Ge | Op::Contains))) = self.peek().cloned() {
            self.i += 1;
            return Ok(Node::Bin(op, Box::new(left), Box::new(self.bit()?)));
        }
        Ok(left)
    }

    fn bit(&mut self) -> Result<Node, String> {
        let mut left = self.primary()?;
        while self.peek() == Some(&Tok::Op(Op::BitAnd)) {
            self.i += 1;
            left = Node::Bin(Op::BitAnd, Box::new(left), Box::new(self.primary()?));
        }
        Ok(left)
    }

    fn primary(&mut self) -> Result<Node, String> {
        match self.next()? {
            Tok::Num(n) => Ok(Node::Num(n)),
            // Negative numbers only (`tlen < -100`); there is no arithmetic.
            Tok::Minus => match self.next()? {
                Tok::Num(n) => Ok(Node::Num(-n)),
                t => Err(format!("'-' goes before a number, not {}", t.describe())),
            },
            Tok::Str(s) => Ok(Node::Str(s)),
            Tok::Open => {
                let inner = self.or()?;
                match self.next()? {
                    Tok::Close => Ok(inner),
                    t => Err(format!("expected ')' but found {}", t.describe())),
                }
            }
            Tok::Ident(w) if w == "true" => Ok(Node::Bool(true)),
            Tok::Ident(w) if w == "false" => Ok(Node::Bool(false)),
            Tok::Ident(w) if w == "tag" && self.peek() == Some(&Tok::Open) => {
                self.i += 1;
                let tag = match self.next()? {
                    Tok::Ident(t) => t,
                    Tok::Str(t) => t.to_ascii_lowercase(),
                    t => return Err(format!("expected a tag name, found {}", t.describe())),
                };
                if self.next()? != Tok::Close {
                    return Err("expected ')' after the tag name".into());
                }
                match tag.as_str() {
                    "hp" => Ok(Node::Field(Field::Hp)),
                    "rg" => Ok(Node::Field(Field::Rg)),
                    _ => Err(format!("tag({}) isn't kept per read; HP and RG are", tag.to_ascii_uppercase())),
                }
            }
            Tok::Ident(w) => field(&w).map(Node::Field).ok_or_else(|| format!("unknown field '{w}'")),
            Tok::Op(_) | Tok::Not | Tok::Close => Err("expected a value or field".into()),
        }
    }
}

/// A pileup read test from a parsed expression: the read is decoded (against the reference
/// when the caller has it) and the expression evaluated, so queries, coverage and evidence
/// keep exactly the reads the view shows.
pub fn predicate(expr: Arc<Expr>, text: &str) -> ggv_query::ReadPredicate {
    ggv_query::ReadPredicate {
        text: text.trim().to_owned(),
        needs_reference: expr.uses_mismatches(),
        test: Arc::new(move |read, reference| expr.matches(&crate::decode(read, reference))),
    }
}

/// An interned expression: cheap to copy and compare, resolved with [`get`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ExprId(u32);

fn table() -> &'static RwLock<Vec<(String, Arc<Expr>)>> {
    static TABLE: OnceLock<RwLock<Vec<(String, Arc<Expr>)>>> = OnceLock::new();
    TABLE.get_or_init(Default::default)
}

/// Parses `text` (or reuses an earlier parse of the same text).
pub fn intern(text: &str) -> Result<ExprId, String> {
    let text = text.trim();
    if let Some(i) = table().read().unwrap().iter().position(|(t, _)| t == text) {
        return Ok(ExprId(i as u32));
    }
    let expr = Arc::new(Expr::parse(text)?);
    let mut t = table().write().unwrap();
    t.push((text.to_owned(), expr));
    Ok(ExprId((t.len() - 1) as u32))
}

pub fn get(id: ExprId) -> Arc<Expr> {
    table().read().unwrap()[id.0 as usize].1.clone()
}

pub fn text(id: ExprId) -> String {
    table().read().unwrap()[id.0 as usize].0.clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read() -> AlignedRead {
        AlignedRead {
            name: "m64011_ccs/123".into(),
            start: 99,
            end: 250,
            reverse: true,
            mapq: 42,
            flags: 0x1 | 0x2 | 0x10 | 0x40 | 0x400,
            blocks: vec![],
            deletions: vec![],
            ref_skips: vec![],
            insertions: vec![],
            mismatches: vec![],
            soft_clip_left: 5,
            soft_clip_right: 0,
            clip_mismatches: vec![],
            haplotype: Some(2),
            read_group: Some("rg1".into()),
            template_length: -300,
            mate: None,
            supplementary: vec![],
            mods: vec![],
            low_quality: vec![],
            tag: None,
            bases: None,
        }
    }

    fn t(text: &str) -> bool {
        Expr::parse(text).unwrap_or_else(|e| panic!("{text}: {e}")).matches(&read())
    }

    #[test]
    fn knows_when_a_reference_is_needed() {
        assert!(Expr::parse("mapq > 20 && !(mismatches > 0)").unwrap().uses_mismatches());
        assert!(Expr::parse("mm == 0").unwrap().uses_mismatches());
        assert!(!Expr::parse("mapq > 20 && deletions > 0").unwrap().uses_mismatches());
    }

    #[test]
    fn fields_and_operators() {
        assert!(t("mapq >= 20 && hp == 2"));
        assert!(!t("mapq >= 20 && !dup"));
        assert!(t("dup"));
        assert!(t("tag(HP) == 2 and rg == \"rg1\""));
        assert!(t("strand == '-' && reverse && !forward"));
        assert!(t("start == 100 && end == 250 && len == 151"));
        assert!(t("insert == 300 && tlen < 0"));
        assert!(t("flag & 0x400"));
        assert!(!t("flag & 0x800"));
        assert!(t("name ~ \"ccs\""));
        assert!(t("(mapq < 10 || proper) && read1 && !read2"));
        assert!(t("clips == 5 && clip_left > clip_right"));
        assert!(t("not secondary"));
        assert!(t("1,000 > mapq"));
        assert!(t("tlen < -100 && tlen == -300 && tlen > - 301"));
        assert!(!t("tlen < -300"));
    }

    #[test]
    fn missing_values_equal_nothing() {
        let mut r = read();
        r.haplotype = None;
        let e = |s: &str| Expr::parse(s).unwrap().matches(&r);
        assert!(!e("hp == 1"));
        assert!(e("hp != 1"));
        assert!(!e("hp > 0"));
        assert!(!e("hp"));
    }

    #[test]
    fn errors_say_what_and_where() {
        assert!(Expr::parse("mapq >=").unwrap_err().contains("ends too soon"));
        assert!(Expr::parse("mapqq > 1").unwrap_err().contains("unknown field 'mapqq'"));
        assert!(Expr::parse("tag(NM) > 1").unwrap_err().contains("HP and RG"));
        assert!(Expr::parse("(mapq > 1").unwrap_err().contains("ends too soon"));
        assert!(Expr::parse("mapq > 1 )").unwrap_err().contains("unexpected"));
        assert!(Expr::parse("name == \"x").unwrap_err().contains("unclosed"));
        assert!(Expr::parse("").is_err());
        assert!(Expr::parse("tlen < -mapq").unwrap_err().contains("'-' goes before a number"));
    }

    #[test]
    fn interning_reuses_ids() {
        let a = intern("mapq > 5").unwrap();
        let b = intern("  mapq > 5 ").unwrap();
        assert_eq!(a, b);
        assert_eq!(text(a), "mapq > 5");
        assert!(get(a).matches(&read()));
        assert!(intern("bogus >").is_err());
    }
}
