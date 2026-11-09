use tiny_skia::{FillRule, Paint, PathBuilder, Pixmap, Rect, Stroke, Transform};

use crate::{Anchor, Color, Prim, Scene, font, text_width};

pub(crate) fn render(scene: &Scene, pixel_ratio: f32) -> Result<Vec<u8>, String> {
    let ratio = pixel_ratio.max(0.25);
    let w = (scene.width as f32 * ratio).round() as u32;
    let h = (scene.height as f32 * ratio).round() as u32;
    let mut pixmap = Pixmap::new(w.max(1), h.max(1)).ok_or("image is too large")?;
    pixmap.fill(sk(scene.background));
    let t = Transform::from_scale(ratio, ratio);

    for prim in &scene.prims {
        match prim {
            Prim::Rect { x, y, w, h, fill } => fill_rect(&mut pixmap, *x, *y, *w, *h, *fill, t),
            Prim::StrokeRect { x, y, w, h, color, width } => {
                if let Some(r) = Rect::from_xywh(*x, *y, *w, *h) {
                    let path = PathBuilder::from_rect(r);
                    let stroke = Stroke { width: *width, ..Stroke::default() };
                    pixmap.stroke_path(&path, &paint(*color, false), &stroke, t, None);
                }
            }
            Prim::Line { x1, y1, x2, y2, color, width } => {
                let mut pb = PathBuilder::new();
                pb.move_to(*x1, *y1);
                pb.line_to(*x2, *y2);
                if let Some(path) = pb.finish() {
                    let stroke = Stroke { width: *width, ..Stroke::default() };
                    pixmap.stroke_path(&path, &paint(*color, false), &stroke, t, None);
                }
            }
            Prim::Polygon { points, fill } => {
                let mut pb = PathBuilder::new();
                for (i, (x, y)) in points.iter().enumerate() {
                    if i == 0 { pb.move_to(*x, *y) } else { pb.line_to(*x, *y) }
                }
                pb.close();
                if let Some(path) = pb.finish() {
                    pixmap.fill_path(&path, &paint(*fill, true), FillRule::Winding, t, None);
                }
            }
            Prim::Text { x, y, text, scale, color, anchor } => {
                let s = f32::from(*scale);
                let width = text_width(text, *scale);
                let x0 = match anchor {
                    Anchor::Start => *x,
                    Anchor::Middle => x - width / 2.0,
                    Anchor::End => x - width,
                };
                for (i, c) in text.chars().enumerate() {
                    let gx = x0 + i as f32 * 6.0 * s;
                    for (row, bits) in font::glyph(c).iter().enumerate() {
                        for col in 0..5 {
                            if bits & (0x10 >> col) != 0 {
                                fill_rect(
                                    &mut pixmap,
                                    gx + col as f32 * s,
                                    y + row as f32 * s,
                                    s,
                                    s,
                                    *color,
                                    t,
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    pixmap.encode_png().map_err(|e| e.to_string())
}

fn fill_rect(pixmap: &mut Pixmap, x: f32, y: f32, w: f32, h: f32, c: Color, t: Transform) {
    if let Some(r) = Rect::from_xywh(x, y, w.max(0.0), h.max(0.0)) {
        pixmap.fill_rect(r, &paint(c, false), t, None);
    }
}

fn paint(c: Color, anti_alias: bool) -> Paint<'static> {
    let mut p = Paint::default();
    p.set_color(sk(c));
    p.anti_alias = anti_alias;
    p
}

fn sk(c: Color) -> tiny_skia::Color {
    tiny_skia::Color::from_rgba8(c.0, c.1, c.2, c.3)
}
