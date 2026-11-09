// Instanced track renderer. Every instance is a rectangle (or arrowhead) whose x extent is in
// bases relative to the loaded data's origin and whose y extent is in points. Panning and
// zooming change only the uniforms, so nothing is re-uploaded until new data arrives.

struct Uniforms {
    size: vec2<f32>,      // canvas size in points
    offset_bp: f32,       // whole bases from the data origin to the canvas's left edge
    // (offset_frac holds the rest: one f32 can't keep a fraction of a base millions of bases out)
    ppb: f32,             // points per base
    scroll_y: f32,        // vertical scroll of the read rows, in points
    scroll_top: f32,      // scrolled instances are hidden above this y (the reads area top)
    min_px: f32,          // one physical pixel, in points
    linear_out: u32,      // 1 when the framebuffer is sRGB (expects linear colour)
    cov_scale: f32,       // coverage: points per depth unit (autoscaling without rebuilding)
    cov_bottom: f32,      // coverage: the baseline, in points
    offset_frac: f32,     // the fraction of a base, in [0, 1)
    _pad: f32,
};

@group(0) @binding(0) var<uniform> u: Uniforms;

struct VsIn {
    @builtin(vertex_index) vi: u32,
    @location(0) x: vec2<f32>,
    @location(1) y: vec2<f32>,
    @location(2) color: u32,
    @location(3) flags: u32,
    @location(4) lod: vec2<f32>,
};

struct VsOut {
    @builtin(position) pos: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) uv: vec2<f32>,
    @location(2) clip_y: f32,
    @location(3) @interpolate(flat) kind: u32,
    // Rounded instances: their size in points (0 for square ones).
    @location(4) @interpolate(flat) size: vec2<f32>,
};

const FLAG_SCROLL: u32 = 1u;   // moves with the read rows' vertical scroll
// bits 1-2: kind (0 rect, 1 arrow pointing right, 2 arrow pointing left)
const FLAG_GAP: u32 = 8u;      // leave a 1 px gap on the right when wide (coverage columns)
const FLAG_MARKER: u32 = 16u;  // a 2 px marker centred on x0 (insertions)
const FLAG_MIN3: u32 = 32u;    // at least 3 px wide, centred (variants)
const FLAG_COV: u32 = 64u;     // y is in depth units above cov_bottom, scaled by cov_scale
const FLAG_ROUND: u32 = 128u;  // rounded corners with anti-aliased edges (sites, reference bases)
const RADIUS: f32 = 2.5;
const TIP: f32 = 4.0;

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    let lo = c / 12.92;
    let hi = pow((c + vec3<f32>(0.055)) / 1.055, vec3<f32>(2.4));
    return select(hi, lo, c <= vec3<f32>(0.04045));
}

@vertex
fn vs_main(v: VsIn) -> VsOut {
    var out: VsOut;
    out.uv = vec2<f32>(0.0);
    out.clip_y = 1.0;
    out.kind = 0u;
    out.color = vec4<f32>(0.0);
    out.size = vec2<f32>(0.0);

    // Level of detail: hidden outside [lod.x, lod.y) points per base.
    if (u.ppb < v.lod.x || u.ppb >= v.lod.y) {
        out.pos = vec4<f32>(-2.0, -2.0, 0.0, 1.0);
        return out;
    }

    // Two triangles: corners (0,0) (1,0) (0,1) / (0,1) (1,0) (1,1).
    let cx = f32((0x32u >> v.vi) & 1u);
    let cy = f32((0x2Cu >> v.vi) & 1u);
    let kind = (v.flags >> 1u) & 3u;

    // Whole bases first: both are integers, so the difference is exact and small near the view.
    var x0 = ((v.x.x - u.offset_bp) - u.offset_frac) * u.ppb;
    var x1 = ((v.x.y - u.offset_bp) - u.offset_frac) * u.ppb;
    if (kind == 1u) {
        x1 = x0 + TIP;
    } else if (kind == 2u) {
        x0 = x1 - TIP;
    } else if ((v.flags & FLAG_MARKER) != 0u) {
        x0 = x0 - 1.0;
        x1 = x0 + 2.0;
    } else {
        if ((v.flags & FLAG_GAP) != 0u && x1 - x0 >= 4.0) {
            x1 = x1 - 1.0;
        }
        if ((v.flags & FLAG_MIN3) != 0u && x1 - x0 < 3.0) {
            let c = (x0 + x1) * 0.5;
            x0 = c - 1.5;
            x1 = c + 1.5;
        }
        if (x1 - x0 < u.min_px) {
            x1 = x0 + u.min_px;
        }
    }

    var y0 = v.y.x;
    var y1 = v.y.y;
    if ((v.flags & FLAG_COV) != 0u) {
        y0 = u.cov_bottom - (u.cov_bottom - y0) * u.cov_scale;
        y1 = u.cov_bottom - (u.cov_bottom - y1) * u.cov_scale;
    }
    let scrolls = (v.flags & FLAG_SCROLL) != 0u;
    if (scrolls) {
        y0 = y0 - u.scroll_y;
        y1 = y1 - u.scroll_y;
    }

    // Clamped to just outside the canvas, not merely to a large number: a long thin quad (an
    // intron line, a deep read) reaching millions of points out was cut by the GPU's clipper, and
    // the two triangles' cut diagonals no longer met, so at some zooms the one pixel row a 1 pt
    // line covers fell in the gap and the whole line vanished.
    let px = clamp(mix(x0, x1, cx), -16.0, u.size.x + 16.0);
    let py = mix(y0, y1, cy);
    if (scrolls) {
        out.clip_y = py - u.scroll_top;
    }

    out.pos = vec4<f32>(px / u.size.x * 2.0 - 1.0, 1.0 - py / u.size.y * 2.0, 0.0, 1.0);

    let c = unpack4x8unorm(v.color);
    var rgb = c.rgb;
    if (u.linear_out != 0u) {
        rgb = srgb_to_linear(rgb);
    }
    out.color = vec4<f32>(rgb * c.a, c.a); // premultiplied, as egui blends
    // Where this corner is across the instance as drawn before clamping (rounded corners and
    // arrowheads are shaped from it).
    out.uv = vec2<f32>(select(cx, (px - x0) / (x1 - x0), x1 - x0 > 1.0e-6), cy);
    out.kind = kind;
    if ((v.flags & FLAG_ROUND) != 0u && kind == 0u) {
        out.size = vec2<f32>(x1 - x0, y1 - y0);
    }
    return out;
}

@fragment
fn fs_main(i: VsOut) -> @location(0) vec4<f32> {
    if (i.clip_y < 0.0) {
        discard;
    }
    let half = abs(i.uv.y - 0.5) * 2.0;
    if (i.kind == 1u && half > 1.0 - i.uv.x) {
        discard;
    }
    if (i.kind == 2u && half > i.uv.x) {
        discard;
    }
    if (i.size.x > 0.0) {
        // Signed distance to a rounded rectangle filling the quad, faded over one pixel.
        let r = min(RADIUS, 0.5 * min(i.size.x, i.size.y));
        let p = abs((i.uv - vec2<f32>(0.5)) * i.size) - (0.5 * i.size - vec2<f32>(r));
        let d = length(max(p, vec2<f32>(0.0))) + min(max(p.x, p.y), 0.0) - r;
        let a = clamp(0.5 - d / u.min_px, 0.0, 1.0);
        if (a <= 0.0) {
            discard;
        }
        return i.color * a;
    }
    return i.color;
}
