#![warn(missing_docs)]


//! # attractor crate
//!
//! Create and Plot Strange Attractors
//! 
//! Create and plot attractors, specifically an attractor described by Julian
//! Sprott in his book "_Strange Attractors: Creating Patterns in Chaos_".
//! <https://sprott.physics.wisc.edu/fractals/booktext/sabook.pdf>
//! 
//! The resulting attractor is discretized (binned) into matrix, making plotting fast.
//! 


use ndarray::prelude::*;
// use ndarray::Array;
// use ndarray_rand::RandomExt;
// use ndarray_rand::rand_distr::Uniform;
use ndarray::Array1;
use inc_stats::Percentiles as Percentiles;
use image::{RgbaImage, ImageBuffer};
use indicatif::ProgressBar;

/// Holds attractor results.
/// 
pub struct AttractorResult {
    /// filename
    pub filename: String,
    /// number of points computed
    pub n: usize,
    /// does the attractor remain in finite bounds
    pub finite: bool
}

/// Attractor input configuration
/// 
pub struct AttractorInput {
    /// coefficient array
    pub a: Array1<f64>,
    /// number of iterations
    pub n: usize,
    /// starting x and y value,
    pub filename: String,
    /// initial values
    pub xy: (f64, f64),
    /// image size
    pub img_size: (usize, usize),
    /// percentile threshold to limit sampling points
    pub perc_threshold: f64
}


/// Scales (interpolates) value and maps onto a different output range
/// 
fn remap(x:&f64, xmin:&f64, xmax:&f64, ymin:&f64, ymax:&f64) -> f64 {
    if xmin == xmax { return 0.5; }
    return ymin + (x - xmin) * (ymax - ymin) / (xmax - xmin);
}

/// Generate strange attractor using equation 7E of Sprott.
/// 
/// This attractor is based on Equation 7E of Sprott (see references).
///
/// $ x_{i+1} = a_{1} + a_{2} * x_{i} +  a_{3} * y_{i} +  a_{4} * |x_{i}|^{a5}  +  a6 * |y_{i}|^{a7} $
///
/// \deqn{y_{i+1} = a_{8} + a_{9} * x_{i} +  a_{10} * y_{i} +  a_{11} * |x_{i}|^{a12}  +  a13 * |y_{i}|^{a14}}
///
/// # references
/// 
/// Julien C. Sprott, "Strange Attractors: Creating Patterns in Chaos", page 418, Equation 7e, 
/// <https://sprott.physics.wisc.edu/fractals/booktext/sabook.pdf>
///
fn sprott_7e(a:&Array1<f64>, &x:&f64, &y:&f64) -> (f64, f64) {
    let a1  = a[0];
    let a2  = a[1];
    let a3  = a[2];
    let a4  = a[3];
    let a5  = a[4];
    let a6  = a[5];
    let a7  = a[6];
    let a8  = a[7];
    let a9  = a[8];
    let a10 = a[9];
    let a11 = a[10];
    let a12 = a[11];
    let a13 = a[12];
    let a14 = a[13];

    let x1:f64 = a1 + a2*x +  a3*y +  a4*(x.abs().powf(a5))  +  a6*((y.abs()).powf(a7));
    let y1:f64 = a8 + a9*x + a10*y + a11*(x.abs().powf(a12)) + a13*((y.abs()).powf(a14));

    return (x1, y1)
}


/// 
// pub fn attractor(a:&Array1<f64>, n:&usize, x0:f64, y0:f64, 
    // img_size: (usize, usize), perc_threshold: f64) -> AttractorResult { 
pub fn attractor(config: &AttractorInput) -> AttractorResult { 
    let a = &config.a;
    let n = config.n;
    let img_size = config.img_size;
    let perc_threshold = config.perc_threshold;
    let (x0, y0) = config.xy;

    let mut xy:Vec<(f64, f64)> = Vec::with_capacity(n);

    xy.push((x0, y0));

    let mut x1: f64 = x0;
    let mut y1: f64 = y0;

    let mut nn: usize = n;
    
    // compute a burn-in of 1000 samples and discard
    let mut finite = true;
    for i in 1..1000 {
        (x1, y1) = sprott_7e(&a, &x1, &y1);
        if x1.is_infinite() || y1.is_infinite() {
            nn = i;
            finite = false;
            break;      
        }
    }

    if !finite {
        return AttractorResult{filename: String::from(""), n: nn, finite: false};
    }

    // compute nn new samples
    let progress_bar = ProgressBar::new((n / 100_000) as u64);
    for i in 1..1_000_000 {
        if i % 100_000 == 0 { progress_bar.inc(1); }
        (x1, y1) = sprott_7e(&a, &x1, &y1);
        if x1.is_infinite() || y1.is_infinite() {
            finite = false;
            nn = i;
            break;      
        }
        xy.push((x1, y1));
    }

    if !finite {
        return AttractorResult{filename: String::from(""), n: nn, finite: false};
    }

    // compute percentiles
    let xrange: Vec<f64>;
    let yrange: Vec<f64>;
    (xrange, yrange) = compute_percentiles(&xy, &perc_threshold);

    // discretize into image
    let mut image: Array2<i32> = Array2::zeros(img_size);
    discretize_image(&mut image, &xy, img_size, &xrange, &yrange);

    let (mut xsize, mut ysize) = img_size; 
    xsize -= 1;
    ysize -= 1;

    for _i in 1_000_000..n {
        if _i % 100_000 == 0 { progress_bar.inc(1); }
        (x1, y1) = sprott_7e(&a, &x1, &y1);
        remap_element(&mut image, &(x1, y1), &xrange[0], &xrange[1], 
            &yrange[0], &yrange[1], &xsize, &ysize)
    }
    progress_bar.finish();

    // convert array to image with colour gradient
    println!("Converting array to image...");
    let imgbuf:RgbaImage = convert_to_image(image);

    // save png image
    println!("Saving image...");
    save_attractor(imgbuf, &config.filename);
    
    let filename = config.filename.to_string();
    return AttractorResult {filename: filename, n:nn, finite: finite};
}



/// Remaps an xy element into a matrix location and increments the pixel count
fn remap_element(image: &mut Array2<i32>, elem: &(f64, f64), 
    xmin: &f64, xmax: &f64, ymin: &f64, ymax: &f64,
    xsize: &usize, ysize: &usize) {
    let (x, y) = elem; 
    let xx:isize = remap(&x, &xmin, &xmax, &0.0, &(*xsize as f64)) as isize;
    let yy:isize = remap(&y, &ymin, &ymax, &0.0, &(*ysize as f64)) as isize;
    if xx >= 0 && yy >= 0 && xx <= *xsize as isize && yy <= *ysize as isize {
        image[[xx as usize, yy as usize]] += 1;
    }
}

/// Discretizes an xy array into an image
fn discretize_image(image: &mut Array2<i32>, xy: &Vec<(f64, f64)>, img_size: (usize, usize), 
    xrange: &Vec<f64>, yrange: &Vec<f64>){
    let (mut xsize, mut ysize) = img_size; 
    xsize -= 1;
    ysize -= 1;
    for elem in xy {
        remap_element(image, &elem, &xrange[0], &xrange[1], 
            &yrange[0], &yrange[1], &xsize, &ysize);
    }
}

/// Convert array to RgbaImage
fn convert_to_image(image: Array2<i32>) -> RgbaImage {
    let shape = image.shape();
    let imgx: u32 = shape[0] as u32;
    let imgy: u32 = shape[1] as u32;
    
    println!("Computing max value...");
    let img_max = *image.iter().max().unwrap() as f32;
    println!("Rescaling values...");
    let raw_img:Array2<f32> = image.mapv(|x| ((x as f32 / img_max) * 1.0_f32.exp()).ln_1p());
    println!("Converting to standard layout...");
    let raw_img = raw_img.reversed_axes();
    
    let grad = colorgrad::rainbow();
    let mut imgbuf: RgbaImage = ImageBuffer::new(imgx, imgy);
    
    println!("Converting to rgba image...");
    
    for (p, v) in imgbuf.pixels_mut()
        .zip(raw_img.iter()) {
            let t = (*v as f64) / 1.0_f64.exp().ln_1p();
            let mut rgba = [0,0,0,255];
            if t != 0.0 {
                rgba = grad.at(t).to_rgba8();
            }
            *p = image::Rgba(rgba);
        }
        
    return imgbuf;
}

/// Save RgbaImage to png file
fn save_attractor(imgbuf: RgbaImage, filename: &String) {
    // let filename = String::from("attractor.png");
    imgbuf.save(filename).unwrap();
}

/// Given a vector and a threshold, compute percentile cutofss
fn compute_percentiles(xy: &Vec<(f64, f64)>, perc_threshold: &f64) -> (Vec<f64>, Vec<f64>) {
    // compute percentiles for xrange and yrange

    let mut xperc = Percentiles::new();
    let mut yperc = Percentiles::new();
    for elem in xy {
        let (x, y) = elem;
        xperc.add(x);
        yperc.add(y);
    }
    let percs = [(1.0 - perc_threshold) / 2.0, perc_threshold + (1.0 - perc_threshold) / 2.0];
    let xrange = xperc.percentiles(&percs).unwrap().unwrap();
    let yrange = yperc.percentiles(&percs).unwrap().unwrap();
    return (xrange, yrange)
}

