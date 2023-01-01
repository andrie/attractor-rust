use ndarray::prelude::*;
// use ndarray::Array;
// use ndarray_rand::RandomExt;
// use ndarray_rand::rand_distr::Uniform;
use ndarray::Array1;
use inc_stats::Percentiles as Percentiles;
// use colorgrad::{rainbow};
use image::{RgbaImage, ImageBuffer};

pub struct Xy {
    pub image: Array2<f32>,
    pub n: usize
}

// fn scale01(x:&f64, xmin:&f64, xmax:&f64) -> f64 {
//     if xmin == xmax { return 0.5; }
//     return (x - xmin) / (xmax - xmin);
// }

fn remap(x:&f64, xmin:&f64, xmax:&f64, ymin:&f64, ymax:&f64) -> f64 {
    if xmin == xmax { return 0.5; }
    return ymin + (x - xmin) * (ymax - ymin) / (xmax - xmin);
}

fn sprott(a:&Array1<f64>, &x:&f64, &y:&f64) -> (f64, f64) {
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

pub fn attractor(a:&Array1<f64>, n:&usize, x0:f64, y0:f64, 
    img_size: (usize, usize), perc_threshold: f64) -> Xy { 

    let mut xy:Vec<(f64, f64)> = Vec::new();

    xy.push((x0, y0));

    let mut x1: f64 = x0;
    let mut y1: f64 = y0;

    let mut nn: usize = *n;
    
    // compute a burn-in of 1000 samples and discard
    let mut finite = true;
    for _i in 1..1000 {
        (x1, y1) = sprott(&a, &x1, &y1);
        if x1.is_infinite() || y1.is_infinite() {
            nn = _i;
            finite = false;
            break;      
        }
    }

    // compute nn new samples
    let mut i:usize = 0;
    while finite && i < *n {
        (x1, y1) = sprott(&a, &x1, &y1);
        if x1.is_infinite() || y1.is_infinite() {
            nn = i;
            finite = false;
            break;      
        }
        xy.push((x1, y1));
        i += 1;
    }

    // compute percentiles for xrange and yrange

    let mut xperc = Percentiles::new();
    let mut yperc = Percentiles::new();
    for elem in &xy {
        let (x, y) = elem;
        xperc.add(x);
        yperc.add(y);
    }
    let percs = [(1.0 - perc_threshold) / 2.0, perc_threshold + (1.0 - perc_threshold) / 2.0];
    let xrange = xperc.percentiles(&percs).unwrap().unwrap();
    let yrange = yperc.percentiles(&percs).unwrap().unwrap();

    // discretize into image
    let mut image: Array2<i32> = Array2::zeros(img_size);
    if finite {
        nn = *n;
        let (mut xsize, mut ysize) = img_size; 
        xsize -= 1;
        ysize -= 1;
        // println!("x and y size: {}, {}", &xsize, &ysize);
        for elem in xy {
            let (x, y) = elem; 
            let xx:isize = remap(&x, &xrange[0], &xrange[1], &0.0, &(xsize as f64)) as isize;
            let yy:isize = remap(&y, &yrange[0], &yrange[1], &0.0, &(ysize as f64)) as isize;
            if xx >= 0 && yy >= 0 && xx <= xsize as isize && yy <= ysize as isize {
                image[[xx as usize, yy as usize]] += 1;
            }
        }
    }

    let image = save_attractor(image);
    
    return Xy {image:image, n:nn};
 }

pub fn save_attractor(image: Array2<i32>) -> Array2<f32> {
    let shape = image.shape();
    let imgx: u32 = shape[0] as u32;
    let imgy: u32 = shape[1] as u32;

    let img_max = *image.iter().max().unwrap() as f32;
    let raw_img:Array2<f32> = image.mapv(|x| ((x as f32 / img_max) * 1.0_f32.exp()).ln_1p());
    raw_img.as_standard_layout();

    let grad = colorgrad::rainbow();
    let mut imgbuf: RgbaImage = ImageBuffer::new(imgx, imgy);

    let mut maxt = 0.0;
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let t = raw_img[[x as usize, y as usize]] / 1.0_f32.exp().ln_1p();
        if t > maxt { maxt = t; }
        let rgba = grad.at(t as f64).to_rgba8();
        *pixel = image::Rgba(rgba);
    }

    println!("maxt = {}", maxt);

    imgbuf.save("attractor.png").unwrap();

    return raw_img;
}
