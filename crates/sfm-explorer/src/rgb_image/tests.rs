// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

use ndarray::{Array3, Axis};

use super::{rgb_to_color_image, rgb_view_to_color_image};

/// Build a `(h, w, 3)` thumbnail whose pixels encode their own coordinates.
fn ramp(height: usize, width: usize) -> Array3<u8> {
    Array3::from_shape_fn((height, width, 3), |(y, x, c)| (y * 10 + x * 3 + c) as u8)
}

#[test]
fn rgb_to_color_image_widens_with_opaque_alpha() {
    let rgb = [1_u8, 2, 3, 4, 5, 6];
    let image = rgb_to_color_image(&rgb, 2, 1);

    assert_eq!(image.size, [2, 1]);
    assert_eq!(image.pixels.len(), 2);
    assert_eq!(image.pixels[0].to_array(), [1, 2, 3, 255]);
    assert_eq!(image.pixels[1].to_array(), [4, 5, 6, 255]);
}

#[test]
fn rgb_view_to_color_image_reads_its_size_from_the_data() {
    // The regression this module exists for: a non-square thumbnail. The old
    // image_browser copy declared [128, 128] regardless, which trips the size
    // assertion inside ColorImage::from_rgba_unmultiplied.
    let thumbnails = ramp(4, 7);
    let image = rgb_view_to_color_image(thumbnails.view());

    assert_eq!(image.size, [7, 4]);
    assert_eq!(image.pixels.len(), 28);
}

#[test]
fn rgb_view_to_color_image_preserves_row_major_order() {
    let thumbnails = ramp(3, 5);
    let image = rgb_view_to_color_image(thumbnails.view());

    for y in 0..3 {
        for x in 0..5 {
            let expected = [
                thumbnails[[y, x, 0]],
                thumbnails[[y, x, 1]],
                thumbnails[[y, x, 2]],
                255,
            ];
            assert_eq!(
                image.pixels[y * 5 + x].to_array(),
                expected,
                "at ({y}, {x})"
            );
        }
    }
}

#[test]
fn rgb_view_to_color_image_handles_a_non_contiguous_view() {
    // `index_axis` on a 4-D thumbnail stack yields a contiguous view, but a
    // sliced view does not; both must produce the same pixels.
    let stack = Array3::from_shape_fn((3, 4, 3), |(y, x, c)| (y * 20 + x * 5 + c) as u8);
    let strided = stack.slice(ndarray::s![.., ..;2, ..]);
    assert!(strided.as_slice().is_none(), "expected a strided view");

    let image = rgb_view_to_color_image(strided);

    assert_eq!(image.size, [2, 3]);
    assert_eq!(image.pixels[0].to_array(), [0, 1, 2, 255]);
    assert_eq!(image.pixels[1].to_array(), [10, 11, 12, 255]);
}

#[test]
fn thumbnail_stack_slice_round_trips() {
    // The shape the panels actually pass: one image out of a (n, h, w, 3) stack.
    let stack = ndarray::Array4::from_shape_fn((2, 4, 6, 3), |(i, y, x, c)| {
        (i * 100 + y * 10 + x * 2 + c) as u8
    });
    let image = rgb_view_to_color_image(stack.index_axis(Axis(0), 1));

    assert_eq!(image.size, [6, 4]);
    assert_eq!(image.pixels[0].to_array(), [100, 101, 102, 255]);
}
