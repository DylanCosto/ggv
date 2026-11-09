//! GPU side of the track view: one instanced draw per track.
//!
//! All tracks share one pipeline. Each track (and the reference row) has its own slot: a
//! uniform buffer and an instance buffer. Every callback writes only its own slot in
//! `prepare`, so several tracks can draw in one frame without clobbering each other.

use std::{collections::HashMap, sync::Arc};

use eframe::wgpu::util::DeviceExt as _;
use eframe::{egui, egui_wgpu, wgpu};

/// One rectangle or arrowhead. 32 bytes; x in bases from the data origin, y in points.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Instance {
    pub x0: f32,
    pub x1: f32,
    pub y0: f32,
    pub y1: f32,
    /// RGBA, red in the lowest byte.
    pub color: u32,
    pub flags: u32,
    /// Shown only while `lod_min <= points per base < lod_max`.
    pub lod_min: f32,
    pub lod_max: f32,
}

const _: () = assert!(size_of::<Instance>() == 32);

pub const FLAG_SCROLL: u32 = 1;
pub const KIND_ARROW_RIGHT: u32 = 1 << 1;
pub const KIND_ARROW_LEFT: u32 = 2 << 1;
pub const FLAG_GAP: u32 = 8;
pub const FLAG_MARKER: u32 = 16;
/// At least 3 px wide, centred (variants, which are often 1 bp).
pub const FLAG_MIN3: u32 = 32;
/// Coverage: y in depth units above the coverage baseline, scaled by `Uniforms::cov_scale`.
pub const FLAG_COV: u32 = 64;
/// Rounded corners with anti-aliased edges (variant sites, reference bases).
pub const FLAG_ROUND: u32 = 128;

/// A GPU slot id, unique for the process: each track (and each pane's reference row) draws
/// into its own slot, so several panes and tracks never overwrite each other's buffers.
pub fn new_slot() -> u64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(1);
    NEXT.fetch_add(1, Ordering::Relaxed)
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Uniforms {
    pub size: [f32; 2],
    /// Whole bases from the data origin to the left edge; `offset_frac` is the rest. Split so
    /// the view lines up with egui to a fraction of a pixel even when the origin is millions of
    /// bases away (one f32 there only resolves a quarter or half base).
    pub offset_bp: f32,
    pub ppb: f32,
    pub scroll_y: f32,
    pub scroll_top: f32,
    pub min_px: f32,
    pub linear_out: u32,
    pub cov_scale: f32,
    pub cov_bottom: f32,
    pub offset_frac: f32,
    pub _pad: f32,
}

const _: () = assert!(size_of::<Uniforms>() % 16 == 0);

/// Instances tagged with a generation, so the GPU buffer is rewritten only when they change.
#[derive(Clone)]
pub struct InstanceSet {
    pub generation: u64,
    pub instances: Arc<Vec<Instance>>,
}

struct Slot {
    uniform_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    /// Instance buffers, each within the device's largest buffer, drawn in order: a track with
    /// more instances than one buffer can hold (millions of reads) is split, not refused.
    chunks: Vec<Chunk>,
    generation: u64,
}

struct Chunk {
    buffer: wgpu::Buffer,
    capacity: u64,
    count: u32,
}

struct GpuState {
    pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    linear_out: bool,
    slots: HashMap<u64, Slot>,
}

impl GpuState {
    fn slot(&mut self, device: &wgpu::Device, id: u64) -> &mut Slot {
        let layout = &self.bind_group_layout;
        self.slots.entry(id).or_insert_with(|| {
            let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ggv_uniforms"),
                contents: bytemuck::bytes_of(&Uniforms::default()),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });
            let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("ggv_uniforms"),
                layout,
                entries: &[wgpu::BindGroupEntry { binding: 0, resource: uniform_buffer.as_entire_binding() }],
            });
            Slot { uniform_buffer, bind_group, chunks: Vec::new(), generation: u64::MAX }
        })
    }
}

/// Creates the pipeline once, at startup.
pub fn init(render_state: &egui_wgpu::RenderState) {
    let device = &render_state.device;
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("ggv_tracks"),
        source: wgpu::ShaderSource::Wgsl(include_str!("tracks.wgsl").into()),
    });

    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("ggv_uniforms"),
        entries: &[wgpu::BindGroupLayoutEntry {
            binding: 0,
            visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }],
    });

    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("ggv_tracks"),
        bind_group_layouts: &[Some(&bind_group_layout)],
        immediate_size: 0,
    });

    let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("ggv_tracks"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_main"),
            buffers: &[Some(wgpu::VertexBufferLayout {
                array_stride: size_of::<Instance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &wgpu::vertex_attr_array![
                    0 => Float32x2, 1 => Float32x2, 2 => Uint32, 3 => Uint32, 4 => Float32x2
                ],
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        primitive: wgpu::PrimitiveState::default(),
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: Some("fs_main"),
            targets: &[Some(wgpu::ColorTargetState {
                format: render_state.target_format,
                blend: Some(wgpu::BlendState::PREMULTIPLIED_ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        multiview_mask: None,
        cache: None,
    });

    render_state.renderer.write().callback_resources.insert(GpuState {
        pipeline,
        bind_group_layout,
        linear_out: render_state.target_format.is_srgb(),
        slots: HashMap::new(),
    });
}

/// Draws `set` into `rect` using slot `slot`. `retire` lists slots of closed tracks to free.
pub fn paint(painter: &egui::Painter, rect: egui::Rect, slot: u64, uniforms: Uniforms, set: InstanceSet, retire: Vec<u64>) {
    painter.add(egui_wgpu::Callback::new_paint_callback(rect, TrackCallback { slot, uniforms, set, retire }));
}

struct TrackCallback {
    slot: u64,
    uniforms: Uniforms,
    set: InstanceSet,
    retire: Vec<u64>,
}

impl egui_wgpu::CallbackTrait for TrackCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        _screen: &egui_wgpu::ScreenDescriptor,
        _encoder: &mut wgpu::CommandEncoder,
        resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        let Some(gpu) = resources.get_mut::<GpuState>() else { return Vec::new() };
        for id in &self.retire {
            gpu.slots.remove(id);
        }

        let mut uniforms = self.uniforms;
        uniforms.linear_out = u32::from(gpu.linear_out);
        let slot = gpu.slot(device, self.slot);
        queue.write_buffer(&slot.uniform_buffer, 0, bytemuck::bytes_of(&uniforms));

        if slot.generation != self.set.generation {
            let size = size_of::<Instance>() as u64;
            // Chunks of at most 256 MB (or the device's largest buffer), as whole instances: the
            // headroom below can then waste at most one chunk's slack.
            let max_bytes = (device.limits().max_buffer_size.min(256 << 20) / size) * size;
            let per_chunk = (max_bytes / size) as usize;
            let parts: Vec<&[Instance]> = self.set.instances.chunks(per_chunk.max(1)).collect();
            slot.chunks.truncate(parts.len().max(1));
            for (k, part) in parts.iter().enumerate() {
                let bytes: &[u8] = bytemuck::cast_slice(part);
                let needed = (bytes.len() as u64).max(size);
                // Reallocated when too small, or four times too large (after zooming in from a
                // dense view), so a buffer doesn't keep its largest size for good.
                if slot.chunks.get(k).is_none_or(|c| c.capacity < needed || c.capacity / 4 > needed) {
                    // Grow with headroom so panning into slightly denser data doesn't reallocate.
                    let capacity = needed.next_power_of_two().min(max_bytes).max(needed);
                    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
                        label: Some("ggv_instances"),
                        size: capacity,
                        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                        mapped_at_creation: false,
                    });
                    let chunk = Chunk { buffer, capacity, count: 0 };
                    if k < slot.chunks.len() {
                        slot.chunks[k] = chunk;
                    } else {
                        slot.chunks.push(chunk);
                    }
                }
                queue.write_buffer(&slot.chunks[k].buffer, 0, bytes);
                slot.chunks[k].count = part.len() as u32;
            }
            if parts.is_empty() {
                for c in &mut slot.chunks {
                    c.count = 0;
                }
            }
            slot.generation = self.set.generation;
        }
        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(gpu) = resources.get::<GpuState>() else { return };
        let Some(slot) = gpu.slots.get(&self.slot) else { return };
        if slot.chunks.iter().all(|c| c.count == 0) {
            return;
        }
        render_pass.set_pipeline(&gpu.pipeline);
        render_pass.set_bind_group(0, &slot.bind_group, &[]);
        for c in slot.chunks.iter().filter(|c| c.count > 0) {
            render_pass.set_vertex_buffer(0, c.buffer.slice(..));
            render_pass.draw(0..6, 0..c.count);
        }
    }
}
