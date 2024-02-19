use ndarray::{Array1, Array2, Axis, array};
use rand::SeedableRng;
use hmmm::HMM;

fn thresholding(matrix: &mut Array2<f64>, threshold: f64) {
    matrix.mapv_inplace(|x| if x > threshold { 1.0 } else { 0.0 });
}

fn sub_to_number(matrix: &Array2<f64>) -> Array1<usize> {
    let mut result = Vec::new();
    for row in matrix.axis_iter(Axis(0)) {
        let mut number = 0;
        for (i, value) in row.iter().enumerate() {
            number += (*value as i32 * 2i32.pow(i as u32)) as usize;
        }
        result.push(number);
    }
    Array1::from_shape_vec(result.len(), result).unwrap()
}

pub fn hmm_train(cwt_matrix: &mut Array2<f64>) -> HMM {
    // Initialize HMM
    // 0 is non-coding region
    // 1 is start codon
    // 2 is exon
    // 3 is splicing donor
    // 4 is intron
    // 5 is splicing acceptor
    // 6 is stop codon
    let init_sequence_answer = array![0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 2, 2, 2,
    6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6];
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    let hmm = HMM::train(&init_sequence, 7, 8, &mut rng);
    // 8 for 2^3

    thresholding(cwt_matrix, 0.5);
    let training_ys = sub_to_number(&cwt_matrix);

    while convergence(&training_ys) {

        let training_ys = hmm.most_likely_sequence(ys)
    }

    
    
    hmm
}

fn convergence(ys: &Array1<usize>) -> bool {
    // criteria one: 
    // TODO: implement criteria one
    false
}