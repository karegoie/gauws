use std::f64::consts::PI;
use rayon::prelude::*;

use ndarray::{Array, Array1, s, Array2, Axis};
use rustfft::{FftPlanner, num_complex::Complex};

use std::sync::{Arc, Mutex};

pub struct Params {
    num: usize,
    tradeoff: f64,
    t_values: Vec<usize>,
}

fn linspace(start: f64, stop: f64, num: usize) -> Array1<f64> {
    let delta = (stop - start) / ((num - 1) as f64);
    Array::from_shape_fn(num, |i| start + (i as f64) * delta)
}

fn psi(t: Array1<f64>, s: f64) -> Array1<Complex<f64>> {
    let kappa_sigma = (-0.5 * s.powf(2.0)).exp();
    let c_sigma = (1.0 + (-s.powf(2.0)).exp() - 2.0*(-0.75*s.powf(2.0)).exp()).powf(-0.5);
    let i_s = Complex::new(0.0, s as f64);
    // Psi = c_sigma * (pi^(-1/4))* e^(-1/2 * t^2) * (e^(i*s*t) - kappa_sigma)
    let coeff = Complex::new(c_sigma * (PI.powf(-0.25)), 0.0);
    let part1 = coeff * t.mapv(|x| ((-0.5 * x * x) + i_s*x).exp());
    let part2 = c_sigma * (PI.powf(-0.25)) * kappa_sigma * t.mapv(|x| (-0.5 * x * x).exp());
    part1 - part2
}
fn wavelet_convolution(tup: (&Array1<f64>, usize), tradeoff: f64) -> Array1<f64> {
    let f = tup.0;
    let tt = tup.1;
    let f_len = f.len();

    let mut f_hat = Array1::zeros(f_len + tt);
    f_hat.slice_mut(s![..f_len]).assign(f);
    let h = psi(linspace(0.0, std::cmp::min(tt, f_len) as f64, tt), tradeoff);
    let mut h_hat = Array1::zeros(f_len + tt);
    h_hat.slice_mut(s![..h.len()]).assign(&h);

    let mut planner = FftPlanner::new();
    let fft_len = f_len + tt;
    let fft = planner.plan_fft_forward(fft_len);
    let ifft = planner.plan_fft_inverse(fft_len);

    let mut f_hat_complex: Vec<Complex<f64>> = f_hat.iter().map(|&val| Complex::new(val, 0.0)).collect();
    let mut h_hat_complex: Vec<Complex<f64>> = h_hat.to_vec();

    fft.process(&mut f_hat_complex);
    fft.process(&mut h_hat_complex);

    let mut result_complex: Vec<Complex<f64>> = f_hat_complex.iter().zip(h_hat_complex.iter()).map(|(&a, &b)| a * b).collect();

    ifft.process(&mut result_complex);

    let result_real: Vec<f64> = result_complex.iter().map(|&val| (val.re*val.re + val.im*val.im).sqrt()).collect();
    let result_view = Array1::from_shape_vec(f_len + tt, result_real).unwrap();
    let start_index = tt / 2;
    let end_index = start_index + f_len;
    result_view.slice(s![start_index..end_index]).to_owned()
}

fn cwt_perform(f: &Array1<f64>, opt: &Params) -> Array2<f64> {
    let f_len = f.len();
    let t_values: Vec<usize> = opt.t_values.clone();
    // linspace(opt.from, opt.to, opt.num).mapv(|x| x as usize).to_vec();

    // Initialize result_2d array
    let result_2d = Arc::new(Mutex::new(Array2::zeros((t_values.len(), f_len))));

    // Perform wavelet convolution and assign results directly to result_2d
    t_values.par_iter().enumerate().for_each(|(i, &t)| {
        let row = wavelet_convolution((&f, t), opt.tradeoff);
        result_2d.lock().unwrap().slice_mut(s![i, ..]).assign(&row);
    });

    let result_cwt_perform = result_2d.lock().unwrap();
    result_cwt_perform.to_owned()
}

pub fn cwt(sig_seqs: &mut Vec<Vec<bool>>, opt: &Params) -> Array2<f64> {
    let sig_seqs: Array2<bool> = Array::from_shape_vec((sig_seqs.len(), sig_seqs[0].len()), sig_seqs.par_iter_mut().flatten().map(|x| *x).collect()).unwrap();

    let mut cwt_result = Array2::<f64>::zeros((opt.num, sig_seqs.dim().0));
    sig_seqs.axis_iter(Axis(1)).for_each( |f| {
        let one_base_result = cwt_perform(&f.map(|&x| if x { 1.0 } else { 0.0 }).to_owned(), opt);
        cwt_result = &cwt_result + &one_base_result;
    });

    cwt_result = cwt_result / 4.0;

    cwt_result
}

pub fn normalize(matrix: &mut Array2<f64>) {
    let mut min = f64::MAX;
    let mut max = f64::MIN;

    for row in matrix.axis_iter(Axis(0)) {
        for value in row.iter() {
            if *value < min {
                min = *value;
            }
            if *value > max {
                max = *value;
            }
        }
    }

    let range = max - min;

    matrix.mapv_inplace(|x| (x - min) / range);
}