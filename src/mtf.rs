//! MTF-1 transform ("fast" variant): a symbol at position 1 moves to the front, any symbol
//! further back moves to position 1. Bit-exact with CoMSA's `CMTFCore`.

pub struct Mtf {
    init_list: [u8; 128],
    init_pos: [u8; 128],
    list: [u8; 128],
    pos: [u8; 128],
}

impl Mtf {
    pub fn new() -> Self {
        let list = crate::initial_symbol_order();
        let mut pos = [0u8; 128];
        for (i, &c) in list.iter().enumerate() {
            pos[c as usize] = i as u8;
        }
        Mtf { init_list: list, init_pos: pos, list, pos }
    }

    pub fn forward(&mut self, src: &[u8], out: &mut Vec<u8>) {
        self.list = self.init_list;
        self.pos = self.init_pos;
        out.clear();
        out.reserve(src.len());
        for &c in src {
            out.push(self.pos[c as usize]);
            self.insert(c);
        }
    }

    pub fn reverse(&mut self, src: &[u8], out: &mut Vec<u8>) {
        self.list = self.init_list;
        self.pos = self.init_pos;
        out.clear();
        out.reserve(src.len());
        for &r in src {
            let c = self.list[r as usize & 127];
            out.push(c);
            self.insert(c);
        }
    }

    #[inline]
    fn insert(&mut self, sym: u8) {
        let p = self.pos[sym as usize] as usize;
        if p == 0 {
            return;
        }
        let to = if p == 1 { 0 } else { 1 };
        self.list.copy_within(to..p, to + 1);
        self.list[to] = sym;
        for i in to..=p {
            self.pos[self.list[i] as usize] = i as u8;
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn roundtrip() {
        let src: Vec<u8> = (0..20000u32).map(|i| b"ACDEFGHIK--..~*"[((i * 2654435761) >> 7) as usize % 15]).collect();
        let (mut ranks, mut back) = (Vec::new(), Vec::new());
        let mut m = super::Mtf::new();
        m.forward(&src, &mut ranks);
        assert!(ranks.iter().all(|&r| r < 56));
        m.reverse(&ranks, &mut back);
        assert_eq!(src, back);
    }
}
