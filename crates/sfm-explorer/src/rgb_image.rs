// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! RGB → RGBA expansion for egui texture uploads.
//!
//! egui uploads textures as `ColorImage`, which is RGBA; every RGB source in
//! this crate — embedded thumbnails, decoded full-resolution images — has to
//! be widened with an opaque alpha channel first. The expansion was written
//! out three times, and one copy hard-coded the thumbnail edge length instead
//! of reading it, which made it the only copy that would panic if that length
//! ever changed (`ColorImage::from_rgba_unmultiplied` asserts that the
//! declared size matches the buffer). Both helpers here take the dimensions
//! from the data, so there is no size to keep in sync.

use ndarray::ArrayView3;

/// Expand packed 3-channel RGB into an opaque `ColorImage`.
///
/// `rgb` must hold `width * height` pixels in row-major order; any trailing
/// bytes that do not complete a pixel are ignored, which
/// [`ColorImage::from_rgba_unmultiplied`](egui::ColorImage::from_rgba_unmultiplied)
/// would otherwise turn into a panic on the length assertion.
pub(crate) fn rgb_to_color_image(rgb: &[u8], width: usize, height: usize) -> egui::ColorImage {
    let mut rgba = Vec::with_capacity(width * height * 4);
    for pixel in rgb.chunks_exact(3).take(width * height) {
        rgba.extend_from_slice(&[pixel[0], pixel[1], pixel[2], 255]);
    }
    // Short sources would otherwise trip the size assertion inside egui.
    rgba.resize(width * height * 4, 0);
    egui::ColorImage::from_rgba_unmultiplied([width, height], &rgba)
}

/// Expand one `(y, x, rgb)` thumbnail slice into an opaque `ColorImage`.
///
/// The width and height come from `slice` itself rather than from the
/// `THUMBNAIL_SIZE` the reconstruction reader happens to use today, so a
/// change to the stored thumbnail size cannot desynchronize the panels from
/// the data.
pub(crate) fn rgb_view_to_color_image(slice: ArrayView3<'_, u8>) -> egui::ColorImage {
    let (height, width) = (slice.shape()[0], slice.shape()[1]);
    match slice.as_slice() {
        Some(contiguous) => rgb_to_color_image(contiguous, width, height),
        None => {
            let owned: Vec<u8> = slice.iter().copied().collect();
            rgb_to_color_image(&owned, width, height)
        }
    }
}

#[cfg(test)]
mod tests;
