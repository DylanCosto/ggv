//! Scripted runs for verification and CI: `--screenshot` saves the real GPU frame once data has
//! loaded; `--bench` drives panning and zooming and reports frame times. Both then exit.

use std::{
    path::PathBuf,
    sync::{Arc, Mutex},
    time::Instant,
};

use eframe::egui;

/// Set when something a scripted run asked for failed (a screenshot or report not written, a
/// file that didn't open): the process then exits with status 1 rather than 0.
static FAILED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Reports a scripted run's failure on stderr and marks the run failed.
pub fn fail(message: &str) {
    eprintln!("ggv-desktop: {message}");
    FAILED.store(true, std::sync::atomic::Ordering::Relaxed);
}

pub fn failed() -> bool {
    FAILED.load(std::sync::atomic::Ordering::Relaxed)
}

pub struct Automation {
    pub screenshot: Option<PathBuf>,
    pub bench: bool,
    /// Zoom all the way out and back in, then report (reproduces zoom-related load bugs).
    pub zoom_cycle: bool,
    /// Set by the app once the zoom-in phase is back at the starting zoom.
    pub cycle_back: Option<Instant>,
    cycle: u32,
    pub report: bool,
    settled: u32,
    shot: Arc<Mutex<Option<Arc<egui::ColorImage>>>>,
    shot_requested: bool,
    frame: u32,
    samples: Vec<(f32, f64)>,
    /// Per sample: whether new data arrived in the previous frame.
    arrived: Vec<bool>,
    last_loads: u64,
    last: Option<Instant>,
    loads_at_start: (u64, u64),
    finished: bool,
    /// When the window started, and when the view was first ready (the initial load).
    started: Instant,
    pub first_ready: Option<std::time::Duration>,
    pub report_loads_done: bool,
}

/// What the app should do to the view this frame.
pub enum Step {
    Wait,
    Pan(f32),
    Zoom(f64),
    /// One frame of the zoom cycle (1-based); the app zooms and calls back when done.
    Cycle(u32),
    Done,
}

const PAN_FRAMES: u32 = 200;
const ZOOM_FRAMES: u32 = 200;
const JUMP_FRAMES: u32 = 200;

impl Automation {
    pub fn new(screenshot: Option<PathBuf>, bench: bool, zoom_cycle: bool) -> Option<Self> {
        (screenshot.is_some() || bench || zoom_cycle).then(|| Self {
            screenshot,
            bench,
            zoom_cycle,
            cycle_back: None,
            cycle: 0,
            report: false,
            settled: 0,
            shot: Arc::default(),
            shot_requested: false,
            frame: 0,
            samples: Vec::new(),
            arrived: Vec::new(),
            last_loads: 0,
            last: None,
            loads_at_start: (0, 0),
            finished: false,
            started: Instant::now(),
            first_ready: None,
            report_loads_done: false,
        })
    }

    /// Advances the script by one frame. `ready` is true once data for the view has loaded.
    /// `loads`: loads so far, and how many of them came from memory.
    pub fn step(&mut self, ctx: &egui::Context, ready: bool, cpu_usage: Option<f32>, loads: (u64, u64)) -> Step {
        if self.finished {
            return Step::Done;
        }
        ctx.request_repaint();
        if self.frame == 0 && (!ready || self.settled < 12) {
            if ready && self.first_ready.is_none() {
                self.first_ready = Some(self.started.elapsed());
            }
            if ready {
                self.settled += 1;
            }
            return Step::Wait;
        }

        if self.zoom_cycle {
            match self.cycle_back {
                None => {
                    self.cycle += 1;
                    return Step::Cycle(self.cycle);
                }
                // Back at the start: give loads up to 15 s to settle, then report.
                Some(t) if !ready && t.elapsed().as_secs() < 15 => return Step::Wait,
                Some(_) => self.report = true,
            }
        }

        if let Some(path) = &self.screenshot {
            if !self.shot_requested {
                let slot = self.shot.clone();
                let ctx2 = ctx.clone();
                ctx.request_screenshot(move |image| {
                    *slot.lock().unwrap() = Some(image);
                    ctx2.request_repaint();
                });
                self.shot_requested = true;
                return Step::Wait;
            }
            let Some(image) = self.shot.lock().unwrap().take() else { return Step::Wait };
            match save_png(&image, path) {
                Ok(()) => println!("screenshot: {} ({}x{})", path.display(), image.size[0], image.size[1]),
                Err(e) => fail(&format!("screenshot {}: {e}", path.display())),
            }
            self.screenshot = None;
        }

        if !self.bench {
            self.finished = true;
            return Step::Done;
        }

        let now = Instant::now();
        if self.frame == 0 {
            self.loads_at_start = loads;
        } else if let (Some(last), Some(cpu)) = (self.last, cpu_usage) {
            self.samples.push((cpu * 1000.0, (now - last).as_secs_f64() * 1000.0));
            self.arrived.push(loads.0 != self.last_loads);
        }
        self.last_loads = loads.0;
        self.last = Some(now);

        let f = self.frame;
        self.frame += 1;
        if f < PAN_FRAMES {
            Step::Pan(0.02) // 2% of the width per frame: stays inside the loaded span
        } else if f < PAN_FRAMES + ZOOM_FRAMES {
            let g = f - PAN_FRAMES;
            Step::Zoom(if (g / 25) % 2 == 0 { 1.04 } else { 1.0 / 1.04 })
        } else if f < PAN_FRAMES + ZOOM_FRAMES + JUMP_FRAMES {
            Step::Pan(0.25) // fast flinging: forces background loads while frames keep coming
        } else {
            self.report(loads.0 - self.loads_at_start.0, loads.1 - self.loads_at_start.1);
            self.finished = true;
            Step::Done
        }
    }

    fn report(&self, loads: u64, from_memory: u64) {
        let phase = |name: &str, range: std::ops::Range<usize>| {
            let s = &self.samples[range.start.min(self.samples.len())..range.end.min(self.samples.len())];
            if s.is_empty() {
                return;
            }
            let mut cpu: Vec<f32> = s.iter().map(|x| x.0).collect();
            let mut wall: Vec<f64> = s.iter().map(|x| x.1).collect();
            cpu.sort_by(f32::total_cmp);
            wall.sort_by(f64::total_cmp);
            let pct = |v: &[f64], p: f64| v[((v.len() - 1) as f64 * p) as usize];
            let cpu64: Vec<f64> = cpu.iter().map(|x| f64::from(*x)).collect();
            println!(
                "{name:<14} frames {:>3}  cpu ms p50 {:>5.2} p95 {:>5.2} max {:>6.2}   wall ms p50 {:>5.2} p95 {:>5.2} max {:>6.2}  (~{:.0} fps)",
                s.len(),
                pct(&cpu64, 0.5),
                pct(&cpu64, 0.95),
                cpu64[cpu64.len() - 1],
                pct(&wall, 0.5),
                pct(&wall, 0.95),
                wall[wall.len() - 1],
                1000.0 / pct(&wall, 0.5),
            );
        };
        let (p, z, j) = (PAN_FRAMES as usize, ZOOM_FRAMES as usize, JUMP_FRAMES as usize);
        phase("pan", 0..p);
        phase("zoom", p..p + z);
        phase("fling+load", p + z..p + z + j);
        println!("background loads during bench: {loads} ({from_memory} served from memory: cache or decoded ahead)");
        if let Some(t) = self.first_ready {
            println!("first view ready after {:.0} ms (window start to everything loaded)", t.as_secs_f64() * 1000.0);
        }
        // Spikes: frames over twice the median CPU time, and how many came with new data.
        let mut cpu: Vec<f32> = self.samples.iter().map(|s| s.0).collect();
        cpu.sort_by(f32::total_cmp);
        if let Some(&median) = cpu.get(cpu.len() / 2) {
            let spikes: Vec<usize> = (0..self.samples.len()).filter(|&i| self.samples[i].0 > 2.0 * median).collect();
            // cpu_usage reports the previous frame, so data arriving one frame earlier counts too.
            let with_data = spikes
                .iter()
                .filter(|&&i| self.arrived[i] || (i > 0 && self.arrived[i - 1]) || self.arrived.get(i + 1).copied().unwrap_or(false))
                .count();
            println!("spikes (> 2x median cpu): {} of {} frames, {} next to a data arrival", spikes.len(), self.samples.len(), with_data);
        }
    }
}

/// The part of a screenshot inside `rect` (points), e.g. just the tracks.
pub fn crop(image: &egui::ColorImage, rect: &egui::Rect, ppp: f32) -> egui::ColorImage {
    let [w, h] = image.size;
    let x0 = ((rect.left() * ppp).round().max(0.0) as usize).min(w);
    let y0 = ((rect.top() * ppp).round().max(0.0) as usize).min(h);
    let x1 = ((rect.right() * ppp).round() as usize).clamp(x0, w);
    let y1 = ((rect.bottom() * ppp).round() as usize).clamp(y0, h);
    let mut cropped = egui::ColorImage::new([x1 - x0, y1 - y0], vec![egui::Color32::TRANSPARENT; (x1 - x0) * (y1 - y0)]);
    for y in y0..y1 {
        let src = &image.pixels[y * w + x0..y * w + x1];
        cropped.pixels[(y - y0) * (x1 - x0)..(y - y0 + 1) * (x1 - x0)].copy_from_slice(src);
    }
    cropped
}

/// Saves the part of a screenshot inside `rect` (points).
#[cfg(test)]
pub fn save_png_cropped(image: &egui::ColorImage, rect: &egui::Rect, ppp: f32, path: &std::path::Path) -> Result<(), String> {
    save_png(&crop(image, rect, ppp), &path.to_path_buf())
}

pub fn encode_png(image: &egui::ColorImage) -> Result<Vec<u8>, String> {
    let [w, h] = image.size;
    let bytes: Vec<u8> = image.pixels.iter().flat_map(|c| c.to_array()).collect();
    let size = tiny_skia::IntSize::from_wh(w as u32, h as u32).ok_or("empty image")?;
    let pixmap = tiny_skia::Pixmap::from_vec(bytes, size).ok_or("bad image")?;
    pixmap.encode_png().map_err(|e| e.to_string())
}

fn save_png(image: &egui::ColorImage, path: &PathBuf) -> Result<(), String> {
    std::fs::write(path, encode_png(image)?).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn export_crops_to_the_tracks_in_physical_pixels() {
        // A 10x6 image at 2 pixels per point; crop the rect (1,1)-(3,2) points = (2,2)-(6,4) px.
        let px = |x: usize, y: usize| egui::Color32::from_rgb(x as u8, y as u8, 0);
        let pixels = (0..6).flat_map(|y| (0..10).map(move |x| px(x, y))).collect();
        let image = egui::ColorImage::new([10, 6], pixels);
        let dir = std::env::temp_dir().join(format!("ggv-export-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("crop.png");
        let rect = egui::Rect::from_min_max(egui::pos2(1.0, 1.0), egui::pos2(3.0, 2.0));
        save_png_cropped(&image, &rect, 2.0, &path).unwrap();
        let png = tiny_skia::Pixmap::load_png(&path).unwrap();
        assert_eq!((png.width(), png.height()), (4, 2));
        let first = png.pixel(0, 0).unwrap();
        assert_eq!((first.red(), first.green()), (2, 2), "top-left is source pixel (2, 2)");
        let last = png.pixel(3, 1).unwrap();
        assert_eq!((last.red(), last.green()), (5, 3));
        let _ = std::fs::remove_dir_all(dir);
    }
}
