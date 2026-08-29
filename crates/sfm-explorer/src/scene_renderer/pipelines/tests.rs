// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! Pipeline-creation tests against a headless `wgpu` device.
//!
//! Every pipeline in this module compiles a WGSL shader. wgpu-core runs naga's
//! full validation on the noop backend, so building them here is a real check
//! that the shaders parse and type-check — the one thing an edit to a `.wgsl`
//! file cannot otherwise fail on until the GUI is launched on a real GPU.

use super::super::gpu_types::{
    gbuffer_targets, PointUniforms, ReconUniforms, COLOR_FORMAT, GBUFFER_DEPTH_STATE,
    LINEAR_DEPTH_FORMAT, PICK_FORMAT,
};

fn device() -> (wgpu::Device, wgpu::Queue) {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::NOOP,
        backend_options: wgpu::BackendOptions {
            noop: wgpu::NoopBackendOptions::enabled(),
            ..Default::default()
        },
        ..wgpu::InstanceDescriptor::new_without_display_handle()
    });
    let adapter =
        pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions::default()))
            .expect("noop adapter");
    pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))
        .expect("noop device")
}

#[test]
fn every_shader_compiles() {
    let (device, _queue) = device();
    super::points::create(&device);
    super::edl::create(&device);
    super::frustum::create(&device);
    let image_quad = super::image_quad::create(&device);
    super::patch::create(&device);
    super::target::create(&device);
    super::track_ray::create(&device);
    let bg_image = super::bg_image::create(&device);
    // These two reuse the bind-group layout of the pipeline they extend.
    super::bg_distorted::create(&device, &bg_image.bind_group_layout);
    super::distorted_quad::create(&device, &image_quad.bind_group_layout);
    device.poll(wgpu::PollType::Poll).expect("device poll");
}

#[test]
fn the_point_uniforms_struct_fits_the_buffer_the_pipeline_allocates() {
    // The buffer is sized from `size_of::<PointUniforms>()`, so a write of the
    // whole struct must land exactly. wgpu validates the bounds, which is what
    // catches a field added on the Rust side without the padding a uniform
    // buffer's 16-byte alignment demands.
    let (device, queue) = device();
    let resources = super::points::create(&device);
    queue.write_buffer(
        &resources.uniform_buffer,
        0,
        bytemuck::bytes_of(&PointUniforms {
            view_proj: [[0.0; 4]; 4],
            view: [[0.0; 4]; 4],
            camera_right: [1.0, 0.0, 0.0],
            _pad0: 0.0,
            camera_up: [0.0, 1.0, 0.0],
            selected_point_index: u32::MAX,
            hovered_point_index: u32::MAX,
            screen_width: 800.0,
            screen_height: 600.0,
            infinity_point_px: 3.0,
            _pad: [0.0; 4],
        }),
    );
    assert_eq!(
        std::mem::size_of::<PointUniforms>() % 16,
        0,
        "a WGSL uniform struct is rounded up to a 16-byte multiple; the Rust \
         side has to match or the tail of the buffer is garbage"
    );
    device.poll(wgpu::PollType::Poll).expect("device poll");
}

#[test]
fn the_recon_uniforms_struct_matches_its_wgsl_layout() {
    // Five shaders declare this block; all five must agree with the Rust
    // definition, and wgpu validates the write against the buffer size.
    let (device, queue) = device();
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("recon uniforms"),
        size: std::mem::size_of::<ReconUniforms>() as u64,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(
        &buffer,
        0,
        bytemuck::bytes_of(&ReconUniforms {
            model: [[0.0; 4]; 4],
            point_size: 1.0,
            point_pick_base: 0,
            image_pick_base: 0,
            pickable: 1,
            tint_color: [0.0; 4],
            show_infinity: 1.0,
            _pad: [0.0; 3],
        }),
    );
    assert_eq!(std::mem::size_of::<ReconUniforms>() % 16, 0);
    // mat4 (64) + four scalars (16) + vec4 (16) + a scalar rounded up to the
    // struct's 16-byte alignment (16). The four scalars before the vec4 are
    // sized to satisfy *its* alignment; the trailing pad satisfies the struct's.
    assert_eq!(std::mem::size_of::<ReconUniforms>(), 112);
    device.poll(wgpu::PollType::Poll).expect("device poll");
}

/// [`GBUFFER_DEPTH_STATE`] is exactly the state the five pass-1 pipelines each
/// wrote out before it existed.
///
/// The hoisted constant spells `stencil` and `bias` out in full, because
/// neither `StencilState::default()` nor `DepthBiasState::default()` is a
/// `const fn`. That makes it a *copy* of two `Default` impls owned by `wgpu`,
/// and a `wgpu` upgrade that changed either one would leave the copy behind
/// without failing to compile — the same hazard the entry-name and
/// camera-registry pins were added for. This is that pin.
#[test]
fn gbuffer_depth_state_matches_the_defaults_it_replaced() {
    assert_eq!(
        GBUFFER_DEPTH_STATE,
        wgpu::DepthStencilState {
            format: wgpu::TextureFormat::Depth32Float,
            depth_write_enabled: Some(true),
            depth_compare: Some(wgpu::CompareFunction::Greater),
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState::default(),
        }
    );
}

/// [`gbuffer_targets`] emits the three attachments in the order the pass binds
/// them, and blends only the colour one.
///
/// Order is the part worth pinning: the formats are distinct, so a *dropped*
/// attachment fails pipeline creation loudly, but a **reordered** one is three
/// well-formed targets in the wrong slots — the shader's `@location(2)` pick ID
/// would be written into the linear-depth attachment, and nothing in the type
/// system objects.
#[test]
fn gbuffer_targets_are_in_pass_attachment_order() {
    let blend = wgpu::BlendState::PREMULTIPLIED_ALPHA_BLENDING;
    let targets = gbuffer_targets(Some(blend));
    let formats: Vec<_> = targets.iter().map(|t| t.as_ref().unwrap().format).collect();
    assert_eq!(
        formats,
        vec![COLOR_FORMAT, LINEAR_DEPTH_FORMAT, PICK_FORMAT]
    );

    // Only @location(0) blends; a blended pick ID is not an ID.
    assert_eq!(targets[0].as_ref().unwrap().blend, Some(blend));
    assert!(targets[1].as_ref().unwrap().blend.is_none());
    assert!(targets[2].as_ref().unwrap().blend.is_none());

    // The opaque pipelines pass None and must still get the same three formats.
    let opaque = gbuffer_targets(None);
    assert!(opaque[0].as_ref().unwrap().blend.is_none());
    let opaque_formats: Vec<_> = opaque.iter().map(|t| t.as_ref().unwrap().format).collect();
    assert_eq!(opaque_formats, formats);
}
