//! Draws GGV views.
//!
//! [`build_scene`] turns a [`ggv_layout::ViewData`] into a [`Scene`]: a flat list of rectangles,
//! lines, polygons and text in pixel coordinates. Backends draw a scene and know nothing about
//! genomics: [`Scene::to_png`] (tiny-skia, CPU, no GPU or display needed) and
//! [`Scene::to_svg`]. A wgpu backend for the interactive viewer will take the same scene, so
//! what an agent snapshots is what a person sees.

mod build;
mod font;
mod png;
pub mod report;
mod svg;

pub use build::{ColorBy, RenderOptions, build_scene};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color(pub u8, pub u8, pub u8, pub u8);

impl Color {
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self(r, g, b, 255)
    }

    pub const fn with_alpha(self, a: u8) -> Self {
        Self(self.0, self.1, self.2, a)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    Start,
    Middle,
    End,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Prim {
    Rect { x: f32, y: f32, w: f32, h: f32, fill: Color },
    StrokeRect { x: f32, y: f32, w: f32, h: f32, color: Color, width: f32 },
    Line { x1: f32, y1: f32, x2: f32, y2: f32, color: Color, width: f32 },
    Polygon { points: Vec<(f32, f32)>, fill: Color },
    /// `y` is the top of the text. `scale` 1 is 7 px tall with a 6 px advance.
    Text { x: f32, y: f32, text: String, scale: u8, color: Color, anchor: Anchor },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Scene {
    pub width: u32,
    pub height: u32,
    pub background: Color,
    pub prims: Vec<Prim>,
}

/// Width in pixels of `text` at `scale` (same metrics in both backends).
pub fn text_width(text: &str, scale: u8) -> f32 {
    let n = text.chars().count() as f32;
    if n == 0.0 { 0.0 } else { (n * 6.0 - 1.0) * f32::from(scale) }
}

impl Scene {
    /// Rasterizes to PNG bytes. `pixel_ratio` 2 gives a sharper image at twice the size.
    pub fn to_png(&self, pixel_ratio: f32) -> Result<Vec<u8>, String> {
        png::render(self, pixel_ratio)
    }

    pub fn to_svg(&self) -> String {
        svg::render(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scene() -> Scene {
        Scene {
            width: 40,
            height: 20,
            background: Color::rgb(255, 255, 255),
            prims: vec![
                Prim::Rect { x: 0.0, y: 0.0, w: 10.0, h: 10.0, fill: Color::rgb(255, 0, 0) },
                Prim::Text { x: 12.0, y: 2.0, text: "A<&".into(), scale: 1, color: Color::rgb(0, 0, 0), anchor: Anchor::Start },
            ],
        }
    }

    #[test]
    fn png_has_signature_and_scales() {
        let png = scene().to_png(2.0).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        // IHDR width and height at bytes 16..24.
        assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 80);
        assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 40);
    }

    #[test]
    fn svg_escapes_text() {
        let svg = scene().to_svg();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains(r#"fill="rgb(255,0,0)""#));
        assert!(svg.contains("A&lt;&amp;"));
    }

    #[test]
    fn text_metrics() {
        assert_eq!(text_width("", 1), 0.0);
        assert_eq!(text_width("ab", 1), 11.0);
        assert_eq!(text_width("ab", 2), 22.0);
    }
}
