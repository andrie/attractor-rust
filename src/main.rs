use ndarray::prelude::*;
// use ndarray_rand::{RandomExt, SamplingStrategy};
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;

pub mod attractor;
use crate::attractor::Xy;

fn main() {
    let n:usize = 1000000;
    let a = Array::random(16, Uniform::new(0., 1.));
    let xy:Xy = attractor::attractor(a, n, 0., 0.);

    println!("n: {}", n);
    println!("x: {}", xy.x[n]);
    println!("y: {}", xy.y[n]);
}
