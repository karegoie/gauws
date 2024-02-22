use ndarray::{Array1, Array2, Axis};

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