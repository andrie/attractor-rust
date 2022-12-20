// use ndarray::prelude::*;
// use ndarray::Array;
// use ndarray_rand::RandomExt;
// use ndarray_rand::rand_distr::Uniform;
use ndarray::Array1;

pub struct Xy {
    pub x: Vec<f64>,
    pub y: Vec<f64>
}

 pub fn attractor(a:Array1<f64>, n:usize, x0:f64, y0:f64) -> Xy { 
    // let a = Array::random(16, Uniform::new(-5., 5.));
    let a1 = a[0];
    let a2 = a[1];
    let a3 = a[2];
    let a4 = a[3];
    let a5 = a[4];
    let a6 = a[5];
    let a7 = a[6];
    let a8 = a[7];
    let a9 = a[8];
    let a10 = a[9];
    let a11 = a[10];
    let a12 = a[11];
    let a13 = a[12];
    let a14 = a[13];

    // println!("a1: {}, n: {}", a1, n);

    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();

    x.push(x0);
    y.push(y0);

    let mut x1: f64;
    let mut y1: f64;

    let mut i:usize = 0;
    while i < n {
        // x1 = a1 + a2*x[i] +  a3*y[i] +  a4*pow(fabs(x[i]), a5)  +  a6*pow(fabs(y[i]),  a7);
        // y1 = a8 + a9*x[i] + a10*y[i] + a11*pow(fabs(x[i]), a12) + a13*pow(fabs(y[i]), a14);
        x1 = a1 + a2*x[i] +  a3*y[i] +  a4*f64::powf(f64::abs(x[i]), a5)  +  a6*f64::powf(f64::abs(y[i]), a7);
        y1 = a8 + a9*x[i] + a10*y[i] + a11*f64::powf(f64::abs(x[i]), a12) + a13*f64::powf(f64::abs(y[i]), a14);
        x.push(x1);
        y.push(y1);
        i += 1;
    }
    
    // return a new data frame
    return Xy {x:x, y:y};
 }