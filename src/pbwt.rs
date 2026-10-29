//! Generalized positional BWT over alignment columns: each column is permuted by the
//! (stable) order of the rows' suffixes seen so far, so similar sequences become adjacent.

pub struct Pbwt {
    prev: Vec<u32>,
    curr: Vec<u32>,
}

impl Pbwt {
    pub fn new(n_rows: usize) -> Self {
        Pbwt { prev: (0..n_rows as u32).collect(), curr: vec![0; n_rows] }
    }

    /// Column in row order -> column in current PBWT order. Symbols must be < 128.
    pub fn forward(&mut self, src: &[u8], dest: &mut Vec<u8>) {
        dest.clear();
        dest.extend(self.prev.iter().map(|&r| src[r as usize]));
        self.reorder(dest);
    }

    /// Column in PBWT order -> column in row order.
    pub fn reverse(&mut self, src: &[u8], dest: &mut Vec<u8>) {
        dest.clear();
        dest.resize(src.len(), 0);
        for (&r, &c) in self.prev.iter().zip(src) {
            dest[r as usize] = c;
        }
        self.reorder(src);
    }

    /// Stable counting sort of the current order by this column's symbols (given in current order).
    fn reorder(&mut self, permuted: &[u8]) {
        let mut start = [0u32; 129];
        for &c in permuted {
            start[(c & 127) as usize + 1] += 1;
        }
        for i in 1..129 {
            start[i] += start[i - 1];
        }
        for (&r, &c) in self.prev.iter().zip(permuted) {
            let s = &mut start[(c & 127) as usize];
            self.curr[*s as usize] = r;
            *s += 1;
        }
        std::mem::swap(&mut self.prev, &mut self.curr);
    }
}
