use ndarray::{Array1, Array2, Axis};
use rand::SeedableRng as _;
use hmmm::HMM;

fn thresholding(matrix: &mut Array2<f64>, threshold: f64) {
    matrix.mapv_inplace(|x| if x > threshold { 1.0 } else { 0.0 });
}

fn sub_to_number(matrix: &Array2<f64>) -> Vec<usize> {
    let mut result = Vec::new();
    for row in matrix.axis_iter(Axis(0)) {
        let mut number = 0;
        for (i, value) in row.iter().enumerate() {
            number += (*value as i32 * 2i32.pow(i as u32)) as usize;
        }
        result.push(number);
    }
    result
}

pub fn hmm_train(cwt_matrix: &mut Array2<f64>) -> HMM {
    thresholding(cwt_matrix, 0.5);
    let training_ys = sub_to_number(&cwt_matrix);
    let mut rng = rand::rngs::StdRng::seed_from_u64(1337);
    let hmm = HMM::train(&training_ys, 3, 2, &mut rng);

    hmm
}