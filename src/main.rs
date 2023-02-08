// use ndarray::array;
// use ndarray::prelude::*;
// use ndarray_rand::RandomExt;
// use ndarray_rand::rand_distr::Uniform;
// extern crate raster;

use attractor::{Attractor, AttractorResult};

use colorgrad::Color;

#[allow(dead_code)]
fn g() -> colorgrad::Gradient {
    return colorgrad::CustomGradient::new()
        .colors(&[
            Color::from_rgba8(255, 0, 0, 255),
            Color::from_rgba8(255, 255, 255, 255),
        ])
        .build()
        .unwrap();
}

fn main() {
    let now = std::time::Instant::now();
    // let a = Array::random(16, Uniform::new(0., 1.));
    // let a = array!(-0.589143526796629, 0.674281478192589, 0.738887963849125, -0.303279513858981, 0.777634792527011, 0.0712083741377903, 0.0940309871410827, -1.17350599974869, -0.780243673606325, 0.536210512392793, 1.60723417066195, -1.2742777716953, -0.0660871384554965, 0.448090598863179);

    // set scaling parameter - higher means larger image with more iterations
    // let scale:i32 = 1;

    // define the attractor configuration
    // let sprott = Attractor {
    //     a,
    //     n: 4_000_000 * scale.pow(2) as usize,
    //     filename: String::from("attractor.png"),
    //     xy: (1.0, 1.0),
    //     img_size: (1920 * scale as usize, 1080 * scale as usize),
    //     perc_threshold: 0.83,
    //     gradient_fn: colorgrad::turbo
    //     // gradient_fn: colorgrad::reds
    //     // gradient_fn: g
    // };

    let mut sprott = Attractor::default();
    sprott.n = 1_000_000;
    println!("Computing attractor with {} iterations", sprott.n);

    let xy: AttractorResult = sprott.compute();

    println!("Completed in {:.2} seconds", now.elapsed().as_secs_f32());

    let nn: usize = xy.n;
    if xy.finite {
        println!("Saved to {}", xy.filename);
    } else {
        println!("Encountered infinity");
        println!("nn: {}", nn);
    }
}
