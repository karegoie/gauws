use rayon::prelude::*;

struct Converter;
trait Convert {
    fn convert(&self, value: u8) -> Vec<bool>;
}

impl Convert for Converter {
    fn convert(&self, value: u8) -> Vec<bool> {
        match value {
            b'A' | b'a' => vec![true, false, false, false],
            b'C' | b'c' => vec![false, true, false, false],
            b'G' | b'g' => vec![false, false, true, false],
            b'T' | b't' => vec![false, false, false, true],
            _ => vec![false, false, false, false],
        }
    }
}

pub fn convert_to_signal(sequence: &mut Vec<u8>) -> Vec<Vec<bool>> {
    let converter = Converter;
    //sequence.par_iter_mut().map(|x| converter.convert(*x)).collect()
    // parallely allocate converted sequence to stack
    let mut converted_sequence = Vec::with_capacity(sequence.len());
    converted_sequence.par_extend(sequence.par_iter_mut().map(|x| converter.convert(*x)));
    converted_sequence
}