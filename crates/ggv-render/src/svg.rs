use std::fmt::Write;

use crate::{Anchor, Color, Prim, Scene};

pub(crate) fn render(scene: &Scene) -> String {
    let mut s = String::with_capacity(64 * scene.prims.len() + 256);
    let _ = write!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" shape-rendering="crispEdges">"#,
        w = scene.width,
        h = scene.height
    );
    let _ = write!(s, r#"<rect width="100%" height="100%" {}/>"#, fill(scene.background));

    for prim in &scene.prims {
        let _ = match prim {
            Prim::Rect { x, y, w, h, fill: c } => {
                write!(s, r#"<rect x="{}" y="{}" width="{}" height="{}" {}/>"#, n(*x), n(*y), n(*w), n(*h), fill(*c))
            }
            Prim::StrokeRect { x, y, w, h, color, width } => write!(
                s,
                r#"<rect x="{}" y="{}" width="{}" height="{}" fill="none" {} stroke-width="{}"/>"#,
                n(*x), n(*y), n(*w), n(*h), stroke(*color), n(*width)
            ),
            Prim::Line { x1, y1, x2, y2, color, width } => write!(
                s,
                r#"<line x1="{}" y1="{}" x2="{}" y2="{}" {} stroke-width="{}"/>"#,
                n(*x1), n(*y1), n(*x2), n(*y2), stroke(*color), n(*width)
            ),
            Prim::Polygon { points, fill: c } => {
                let pts: Vec<String> = points.iter().map(|(x, y)| format!("{},{}", n(*x), n(*y))).collect();
                write!(s, r#"<polygon points="{}" {} shape-rendering="auto"/>"#, pts.join(" "), fill(*c))
            }
            Prim::Text { x, y, text, scale, color, anchor } => {
                let size = 8.5 * f32::from(*scale);
                let anchor = match anchor {
                    Anchor::Start => "start",
                    Anchor::Middle => "middle",
                    Anchor::End => "end",
                };
                write!(
                    s,
                    r#"<text x="{}" y="{}" font-family="ui-monospace,Consolas,Menlo,monospace" font-size="{}" text-anchor="{anchor}" {}>{}</text>"#,
                    n(*x),
                    n(y + 7.0 * f32::from(*scale)),
                    n(size),
                    fill(*color),
                    escape(text)
                )
            }
        };
    }
    s.push_str("</svg>\n");
    s
}

fn n(v: f32) -> String {
    let r = (v * 100.0).round() / 100.0;
    if r.fract() == 0.0 { format!("{}", r as i64) } else { format!("{r}") }
}

fn fill(c: Color) -> String {
    paint("fill", c)
}

fn stroke(c: Color) -> String {
    paint("stroke", c)
}

fn paint(attr: &str, c: Color) -> String {
    if c.3 == 255 {
        format!(r#"{attr}="rgb({},{},{})""#, c.0, c.1, c.2)
    } else {
        format!(r#"{attr}="rgb({},{},{})" {attr}-opacity="{:.3}""#, c.0, c.1, c.2, f32::from(c.3) / 255.0)
    }
}

fn escape(t: &str) -> String {
    t.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}
