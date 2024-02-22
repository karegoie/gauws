use ndarray::{Array1, Array2, Axis};
use super::hmm;

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

pub fn hmm_training_func(cwt_matrix: &mut Array2<f64>, max_iter: usize) -> () {
    // Initialize HMM
    // 0 is non-coding region
    // 1 is start codon
    // 2 is exon
    // 3 is splicing donor
    // 4 is intron
    // 5 is splicing acceptor
    // 6 is stop codon
    let init_sequence_answer = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 2, 2, 2,
    6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 2, 2, 2, 2, 
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6];
    // 0~7 values
    let init_sequence_obsveration = vec![0, 1, 2, 1, 2, 0,
    0, 1, 2, 0, 4, 5, 6, 7, 2, 3, 4, 2, 1, 4, 5, 6, 1, 4,
    0, 3, 2, 3, 1, 5, 6, 7, 5, 5, 4, 1, 2, 3, 4, 5, 6, 7, 
    0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1,
    2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3,
    4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5,
    6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
    0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1,
    2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3];

    //let mut hmm = hmm::Hmm::train(&init_sequence_answer, &init_sequence_obsveration, 0.001, 1000, 0.0001);

    //thresholding(cwt_matrix, 0.5);
    //let training_ys = sub_to_number(&cwt_matrix);

    //let mut old_convergence = hmm.evaluate(&training_ys).log10() / training_ys.len() as f64;
    //for _ in 0..max_iter {
        //hmm.learn(&training_ys, 0.0);
        //let new_convergence = hmm.evaluate(&training_ys).log10() / training_ys.len() as f64;

        //if new_convergence - old_convergence < 0.001{
        //    break;
        //}
        //old_convergence = new_convergence;
    //}
    //hmm
}
