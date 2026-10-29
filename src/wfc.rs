//! Weighted frequency count (WFC) transform, bit-exact with CoMSA's `CWFCCore`
//! (weight function "Deo w9", p = 4.0, history window 16384).
//!
//! Each column starts from the same initial symbol list. A symbol's weight is the sum,
//! over its recent occurrences, of a step function of the occurrence's age; the
//! list is kept sorted by weight and the output is the symbol's position in it.

use std::sync::OnceLock;

const MAX_TIME: usize = 4096 * 4;
const P: f64 = 4.0;

/// (age, weight delta) pairs: when an occurrence reaches `age`, its weight changes by `delta`.
fn updates() -> &'static [(usize, i32)] {
    static TABLE: OnceLock<Vec<(usize, i32)>> = OnceLock::new();
    TABLE.get_or_init(build_updates)
}

fn build_updates() -> Vec<(usize, i32)> {
    // Weight divisors; the factor order matches the C++ expression for identical rounding.
    let mut divs = vec![(1usize, 1.0f64)];
    for i in 2..=MAX_TIME {
        let f = |t: usize, v: f64| if i >= t { v } else { 1.0 };
        let div = P * i as f64 * f(4, 1.4) * f(8, 1.1) * f(16, 1.1) * f(32, 1.1) * f(64, 1.1)
            * f(1024, 1.1) * f(2048, 1.1) * f(4096, 4.0) * f(16384, 4.0);
        divs.push((i, div));
    }

    // Round down to a power of 2
    let round_pow2 = |x: f64| -> i64 {
        let mut lx = x as i64;
        while lx & (lx - 1) != 0 {
            lx &= lx - 1;
        }
        lx
    };

    let max_div = round_pow2(divs.last().unwrap().1) as u32;
    let mut ups = vec![(1usize, max_div as i32)];
    let mut p_div = round_pow2(divs[0].1);
    for &(age, div) in &divs[1..] {
        let c_div = round_pow2(div);
        if c_div != p_div {
            ups.push((age, (max_div as i64 / c_div - max_div as i64 / p_div) as i32));
            p_div = c_div;
        }
    }
    let sum = ups.iter().fold(0u32, |s, &(_, d)| s.wrapping_add(d as u32));
    ups.push((MAX_TIME + 1, (sum as i32).wrapping_neg()));
    ups
}

pub struct Wfc {
    init_list: [(u8, i32); 128],
    init_pos: [u8; 128],
    list: [(u8, i32); 128],
    pos: [u8; 128],
    history: Vec<u8>,
    history_size: usize,
    n_active: usize,
    updates: &'static [(usize, i32)],
}

impl Wfc {
    pub fn new() -> Self {
        let order = crate::initial_symbol_order();
        let mut init_list = [(0u8, 0i32); 128];
        let mut init_pos = [0u8; 128];
        for (i, &c) in order.iter().enumerate() {
            init_list[i] = (c, 0);
            init_pos[c as usize] = i as u8;
        }
        Wfc {
            init_list,
            init_pos,
            list: init_list,
            pos: init_pos,
            history: Vec::new(),
            history_size: 0,
            n_active: 0,
            updates: updates(),
        }
    }

    fn reset(&mut self) {
        self.list = self.init_list;
        self.pos = self.init_pos;
        self.history.clear();
        self.history_size = 0;
        self.n_active = 0;
    }

    /// Symbols -> ranks. Symbols must be < 128.
    pub fn forward(&mut self, src: &[u8], out: &mut Vec<u8>) {
        self.reset();
        out.clear();
        out.reserve(src.len());
        for &c in src {
            out.push(self.pos[c as usize]);
            self.insert(c);
        }
    }

    /// Ranks -> symbols. Ranks must be < 128.
    pub fn reverse(&mut self, src: &[u8], out: &mut Vec<u8>) {
        self.reset();
        out.clear();
        out.reserve(src.len());
        for &r in src {
            let c = self.list[r as usize & 127].0;
            out.push(c);
            self.insert(c);
        }
    }

    #[inline]
    fn insert(&mut self, x: u8) {
        self.history.push(x);
        if self.history_size < MAX_TIME {
            self.history_size += 1;
            if self.n_active < self.updates.len() && self.history_size == self.updates[self.n_active].0 {
                self.n_active += 1;
            }
        }

        // Apply the weight changes of all occurrences crossing an age boundary,
        // merging consecutive changes to the same symbol.
        let ups: &'static [(usize, i32)] = self.updates;
        let hist = &self.history[..];
        let hp = hist.len();
        let (list, pos) = (&mut self.list, &mut self.pos);
        let mut p_sym = x;
        let mut p_val = ups[0].1;
        for &(age, delta) in &ups[1..self.n_active] {
            let c_sym = hist[hp - age];
            if c_sym == p_sym {
                p_val += delta;
            } else {
                apply(list, pos, p_sym, p_val);
                p_sym = c_sym;
                p_val = delta;
            }
        }
        apply(list, pos, p_sym, p_val);
    }
}

/// Add `val` to the weight of `sym` and restore the list order.
#[inline]
fn apply(list: &mut [(u8, i32); 128], pos: &mut [u8; 128], sym: u8, val: i32) {
    let mut i = pos[sym as usize & 127] as usize & 127;
    list[i].1 += val;
    if val < 0 {
        while i < 127 && list[i].1 < list[i + 1].1 {
            list.swap(i, i + 1);
            pos[list[i].0 as usize & 127] = i as u8;
            i += 1;
        }
    } else {
        while i > 0 && list[i].1 >= list[i - 1].1 {
            list.swap(i, i - 1);
            pos[list[i].0 as usize & 127] = i as u8;
            i -= 1;
        }
    }
    pos[sym as usize & 127] = i as u8;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Table produced by CoMSA 1.1's CWFCCore(9, 16384, 4.0).
    #[test]
    fn matches_cpp_table() {
        let expect: &[(usize, i32)] = &[(1, 2097152), (2, -1835008), (4, -131072), (6, -65536), (11, -32768), (19, -16384), (35, -8192), (64, -4096), (125, -2048), (250, -1024), (500, -512), (1000, -256), (1817, -128), (3303, -64), (4096, -48), (6606, -8), (13212, -4), (16384, -3), (16385, -1)];
        assert_eq!(updates(), expect);
    }

    #[test]
    fn roundtrip() {
        let src: Vec<u8> = (0..50000u32).map(|i| b"ACDEFGHIK--..~*"[((i * 2654435761) >> 7) as usize % 15]).collect();
        let (mut ranks, mut back) = (Vec::new(), Vec::new());
        let mut w = Wfc::new();
        w.forward(&src, &mut ranks);
        assert!(ranks.iter().all(|&r| r < 56));
        w.reverse(&ranks, &mut back);
        assert_eq!(src, back);
    }
}
