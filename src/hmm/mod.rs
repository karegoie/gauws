use std::collections::{HashMap, VecDeque};

pub struct Hmm {
    states: Vec<usize>,
    symbols: Vec<usize>,
    start_prob: HashMap<usize, f64>,
    trans_prob: HashMap<usize, HashMap<usize, f64>>,
    emit_prob: HashMap<usize, HashMap<usize, f64>>,
}

impl Hmm {
    fn _normalize_prob(prob: &HashMap<usize, f64>, item_set: &Vec<usize>) -> HashMap<usize, f64> {
        let mut result = HashMap::new();
        //if prob is empty
        if prob.is_empty() {
            let number = item_set.len();
            for item in item_set.iter() {
                result.insert(*item, 1.0 / number as f64);
            }
        }
        let mut prob_sum = 0.0;
        for item in item_set.iter() {
            prob_sum += prob[item];
        }
        if prob_sum > 0.0 {
            for item in item_set.iter() {
                result.insert(*item, prob[item] / prob_sum);
            }
        }
        else {
            for item in item_set.iter() {
                result.insert(*item, 0.0);
            }
        }
        result
    }

    fn _normalize_prob_two_dim(prob: &HashMap<usize, HashMap<usize, f64>>, item_set1: &Vec<usize>, item_set2: &Vec<usize>) -> HashMap<usize, HashMap<usize, f64>>{
        let mut result = HashMap::new();
        if prob.is_empty() {
            for item in item_set1.iter() {
                result.insert(*item, Hmm::_normalize_prob(&HashMap::new(), item_set2));
            }
        }
        else {
            for item in item_set1.iter() {
                result.insert(*item, Hmm::_normalize_prob(&prob[item], item_set2));
            }
        }
        
        result
    }

    fn _count(item: usize, count: &mut HashMap<usize, f64>) {
        if !count.contains_key(&item) {
            count.insert(item, 0.0);
            
        }
        *count.get_mut(&item).unwrap() += 1.0;
    }

    fn _count_two_dim(item1: usize, item2: usize, count: &mut HashMap<usize, HashMap<usize, f64>>) {
        if !count.contains_key(&item1) {
            count.insert(item1, HashMap::new());
        }
        Hmm::_count(item2, count.get_mut(&item1).unwrap());
    }

    pub fn new(states: Vec<usize>, symbols: Vec<usize>, start_prob: HashMap<usize, f64>,
        trans_prob: HashMap<usize, HashMap<usize, f64>>, emit_prob: HashMap<usize, HashMap<usize, f64>>) -> Self {
            // unique values
            let mut states = states;
            states.sort();
            states.dedup();
            let mut symbols = symbols;
            symbols.sort();
            symbols.dedup();
            let start_prob = Hmm::_normalize_prob(&start_prob, &states);
            let trans_prob = Hmm::_normalize_prob_two_dim(&trans_prob, &states, &states);
            let emit_prob = Hmm::_normalize_prob_two_dim(&emit_prob, &states, &symbols);
        Hmm {
            states,
            symbols,
            start_prob,
            trans_prob,
            emit_prob,
        }
    }

    pub fn forward(&self, ys: &Vec<usize>) -> Vec<HashMap<usize, f64>> {
        let seq_len = ys.len();
        if seq_len ==  0 {
            return Vec::new();
        }
        let mut alpha = Vec::new();
        alpha.push(HashMap::new());
        for state in self.states.iter() {
            alpha[0].insert(*state, self.start_prob[state] * self.emit_prob[state][&ys[0]] as f64);
        }

        for index in 1..seq_len {
            alpha.push(HashMap::new());
            for state_to in self.states.iter() {
                let mut prob = 0.0;
                for state_from in self.states.iter() {
                    prob += alpha[index - 1][state_from] * self.trans_prob[state_from][state_to]
                }
                alpha[index].insert(*state_to, prob * self.emit_prob[state_to][&ys[index]] as f64);
            }
        }
        return alpha;
    }

    pub fn backward(&self, ys: &Vec<usize>) -> VecDeque<HashMap<usize, f64>> {
        let seq_len = ys.len();
        if seq_len == 0 {
            return VecDeque::new();
        }

        let mut beta = VecDeque::new();
        beta.push_back(HashMap::new());
        for state in self.states.iter() {
            beta[0].insert(*state, 1.0);
        }

        for index in (1..seq_len).rev() {
            beta.push_front(HashMap::new());
            for state_from in self.states.iter() {
                let mut prob = 0.0;
                for state_to in self.states.iter() {
                    prob += beta[1][state_to] * self.trans_prob[state_from][state_to] * 
                    self.emit_prob[state_to][&ys[index]];
                }
                beta[0].insert(*state_from, prob);
            }
        }
        return beta;
    }
    
    pub fn decode(&self, ys: &Vec<usize>) -> VecDeque<usize> {
        let seq_len = ys.len();
        if seq_len == 0 {
            return VecDeque::new();
        }

        // delta is viterbi probability
        let mut delta = HashMap::new();

        // calculate delta in start position
        for state in self.states.iter() {
            delta.insert(*state, self.start_prob[state] * self.emit_prob[state][&ys[0]] as f64);
        }

        // pre is backtracer
        let mut pre = Vec::new();
        for index in 1..seq_len {
            let mut delta_bar = HashMap::new();
            let mut pre_state = HashMap::new();
            for state_to in self.states.iter() {
                let mut max_prob = 0.0;
                let mut max_state = 0;
                for state_from in self.states.iter() {
                    let prob: f64 = delta[state_from] * self.trans_prob[state_from][state_to];
                    if prob > max_prob {
                        max_prob = prob;
                        max_state = *state_from;
                    }
                }
                delta_bar.insert(*state_to, max_prob * self.emit_prob[state_to][&ys[index]] as f64);
                pre_state.insert(*state_to, max_state);
            }
            delta = delta_bar;
            pre.push(pre_state);
        }

        // backtrace
        let mut max_prob = 0.0;
        let mut max_state = 0;
        for state in self.states.iter() {
            if delta[state] > max_prob {
                max_prob = delta[state];
                max_state = *state;
            }
        }

        if max_state == 0 {
            return VecDeque::new();
        }

        let mut result = VecDeque::new();
        result.push_back(max_state);
        for index in (1..seq_len).rev() {
            max_state = pre[index-1][&max_state];
            result.push_front(max_state);
        }
        return result;
    }

    pub fn learn(&mut self, ys: &Vec<usize>, smoothing: f64) -> &mut Self {
        let seq_len = ys.len();

        let alpha = self.forward(ys);
        let beta = self.backward(ys);

        let mut gamma:Vec<Vec<f64>> = Vec::new();
        for index in 0..seq_len {
            let mut prob_sum = 0.0;
            gamma.push(Vec::new());
            for state in self.states.iter() {
                let prob = alpha[index][state] * beta[index][state];
                gamma[index][*state] = prob;
                prob_sum += prob;
            }

            if prob_sum == 0.0 {
                continue;
            }

            for state in self.states.iter() {
                gamma[index][*state] /= prob_sum;
            }
        }

        let mut xi: Vec<Vec<Vec<f64>>> = Vec::new();
        for index in 0..seq_len - 1 {
            let mut prob_sum = 0.0;
            xi.push(Vec::new());
            for state_from in self.states.iter() {
                xi.push(Vec::new());
                for state_to in self.states.iter() {
                    let prob = alpha[index][state_from] *
                    beta[index + 1][state_to] *
                    self.trans_prob[state_from][state_to] *
                    self.emit_prob[state_to][&ys[index + 1]];
                    xi[index][*state_from][*state_to] = prob;
                    prob_sum += prob;
                }
            }
            if prob_sum == 0.0 {
                continue;
            }

            for state_from in self.states.iter() {
                for state_to in self.states.iter() {
                    xi[index][*state_from][*state_to] /= prob_sum;
                }
            }
        }

        let states_len = self.states.len();
        let symbols_len = self.symbols.len();
        for state in self.states.iter() {
            self.start_prob.insert(*state, 
            (gamma[0][*state] + smoothing) / (1.0 + states_len as f64 * smoothing));
            
            let mut gamma_sum = 0.0;
            for index in 0..seq_len - 1 {
                gamma_sum += gamma[index][*state];
            }

            if gamma_sum > 0.0 {
                let denominator = gamma_sum + smoothing * states_len as f64;
                for state_to in self.states.iter() {
                    let mut xi_sum = 0.0;
                    for index in 0..seq_len - 1 {
                        xi_sum += xi[index][*state][*state_to];
                    }
                    self.trans_prob.get_mut(state).unwrap().insert(*state_to, (xi_sum + smoothing) / denominator);
                }
            } 
            else {
                for state_to in self.states.iter() {
                    self.trans_prob.get_mut(state).unwrap().insert(*state_to, 0.0);
                }
            }
            
            // update emit_prob
            gamma_sum += gamma[seq_len - 1][*state];
            let mut emit_gamma_sum = HashMap::new();
            for symbol in self.symbols.iter() {
                emit_gamma_sum.insert(*symbol, 0.0);
            }

            for index in 0..seq_len {
                if let Some(value) = emit_gamma_sum.get_mut(&ys[index]) {
                    *value += gamma[index][*state];
                }
            }

            if gamma_sum > 0.0 {
                let denominator = gamma_sum + smoothing * symbols_len as f64;
                for symbol in self.symbols.iter() {
                    self.emit_prob.get_mut(state).unwrap().insert(*symbol, 
                    (emit_gamma_sum[symbol] + smoothing) / denominator);
                }
            } 
            else {
                for symbol in self.symbols.iter() {
                    self.emit_prob.get_mut(state).unwrap().insert(*symbol, 0.0);
                }
            }
        }
        self
    }

    pub fn evaluate(&self, ys: &Vec<usize>) -> f64 {
        let seq_len = ys.len();
        if seq_len == 0 {
            return 0.0;
        }

        let alpha = self.forward(ys);
        let mut prob = 0.0;
        for state in self.states.iter() {
            prob += alpha[seq_len - 1][state];
        }
        return prob;
    }

    pub fn get_init_model(xs: &Vec<Vec<usize>>, ys: &Vec<Vec<usize>>) -> Self {
        let mut symbol_count = HashMap::new();
        let mut state_count = HashMap::new();
        let mut state_start_count = HashMap::new();
        let mut state_symbol_count = HashMap::new();
        let mut state_trans_count = HashMap::new();

        for (state_list, symbol_list) in xs.iter().zip(ys.iter()) {
            let mut pre_state = 0;
            for (state, symbol) in state_list.iter().zip(symbol_list.iter()) {
                Hmm::_count(*state, &mut state_count);
                Hmm::_count(*symbol, &mut symbol_count);
                Hmm::_count_two_dim(*state, *symbol, &mut state_symbol_count);
                if pre_state == 0 {
                    Hmm::_count(*state, &mut state_start_count);
                }
                else {
                    Hmm::_count_two_dim(pre_state, *state, &mut state_trans_count)
                }
                pre_state = *state;
            }
        }
        Hmm::new(
        state_count.keys().map(|x| *x).collect(),
        symbol_count.keys().map(|x| *x).collect(),
        state_start_count, state_trans_count, state_symbol_count
        )
}

    pub fn train(xs: &Vec<Vec<usize>>, ys: &Vec<Vec<usize>>, tol: f64, max_iter: usize, smoothing: f64) -> Self {
        let mut model = Hmm::get_init_model(xs, ys);
        
        let mut old_likelihood = 0.0;

        for y in ys.iter() {
            old_likelihood += model.evaluate(y).log10();
        }

        old_likelihood = old_likelihood / ys.len() as f64;

        for _ in 0..max_iter {
            let mut new_likelihood = 0.0;
            for y in ys.iter() {
                model.learn(y, smoothing);
                new_likelihood += model.evaluate(y).log10();
            }
            new_likelihood = new_likelihood / ys.len() as f64;

            if (new_likelihood - old_likelihood).abs() < tol {
                break;
            }
            old_likelihood = new_likelihood;
        }
        return model
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_train() {
        let xs = vec![vec![0, 1, 0, 1, 0, 1], vec![1, 0, 1, 0]];
        let ys = vec![vec![1, 2, 3, 1, 2, 3], vec![1, 2, 3, 1]];
        
        let tol = 0.01;
        let max_iter = 100;
        let smoothing = 0.001;

        let model = Hmm::train(&xs, &ys, tol, max_iter, smoothing);

        // Add assertions based on your expected output
        // For example, if you expect the model to have certain states, you can assert like this:
        // assert_eq!(model.states, vec![1, 2, 3, 4]);
        println!("{:?}", model.states);
    }

    #[test]
    fn test_evaluate() {
        let xs = vec![vec![0, 1, 0, 1, 0, 1], vec![1, 0, 1, 0]];
        let ys = vec![vec![1, 2, 3, 1, 2, 3], vec![1, 2, 3, 1]];
        
        let tol = 0.01;
        let max_iter = 10;
        let smoothing = 0.1;

        let model = Hmm::train(&xs, &ys, tol, max_iter, smoothing);

        let prob = model.evaluate(&vec![1, 2, 3]);

        // Add assertions based on your expected output
        // For example, if you expect the probability to be a certain value, you can assert like this:
        // assert_eq!(prob, expected_prob);
        println!("{:?}", prob);
    }

    #[test]
    fn test_get_init_model() {
        let xs = vec![vec![0, 1, 0, 1, 0, 1], vec![1, 0, 1, 0]];
        let ys = vec![vec![1, 2, 3, 1, 2, 3], vec![1, 2, 3, 1]];

        let model = Hmm::get_init_model(&xs, &ys);

        // Add assertions based on your expected output
        // For example, if you expect the model to have certain states, you can assert like this:
        // assert_eq!(model.states, vec![1, 2, 3, 4]);
        println!("{:?}", model.states);
    }
}