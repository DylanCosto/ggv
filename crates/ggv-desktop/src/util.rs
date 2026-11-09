pub fn commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// The smallest 1, 2 or 5 × 10ⁿ at least `raw`.
pub fn nice_step(raw: f64) -> u64 {
    let raw = raw.max(1.0);
    let mag = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|m| m * mag)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * mag)
        .round()
        .max(1.0) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formatting() {
        assert_eq!(commas(1_895_892), "1,895,892");
        assert_eq!(nice_step(13.0), 20);
        assert_eq!(nice_step(4100.0), 5000);
    }
}
