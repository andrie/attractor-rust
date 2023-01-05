use ndarray::array;
// use ndarray::prelude::*;
// use ndarray_rand::RandomExt;
// use ndarray_rand::rand_distr::Uniform;
// extern crate raster;

pub mod attractor;
use crate::attractor::{AttractorResult, AttractorInput};


fn main() {
    println!("Running code");
    // let a = Array::random(16, Uniform::new(0., 1.));
    let a = array!(-0.589143526796629, 0.674281478192589, 0.738887963849125, -0.303279513858981, 0.777634792527011, 0.0712083741377903, 0.0940309871410827, -1.17350599974869, -0.780243673606325, 0.536210512392793, 1.60723417066195, -1.2742777716953, -0.0660871384554965, 0.448090598863179);
    let config = AttractorInput {
        a: a,
        n: 2_000_000,
        filename: String::from("attractor.png"),
        xy: (1.0, 1.0),
        img_size: (1024 * 4, 768 * 4), 
        perc_threshold: 0.7
    };
    let xy:AttractorResult = attractor::attractor(&config);

    let nn: usize = xy.n;
    if xy.finite {
        println!("saved to {}", xy.filename);
    } else {
        println!("Encountered infinity");
        println!("nn: {}", nn);
    }

}