//! Carryless range coder (Subbotin) with simple adaptive frequency models, bit-exact with
//! CoMSA's `CRangeEncoder` / `CRangeDecoder` / `CSimpleModel`.

use crate::error::{Result, corrupt};

const TOP: u32 = 0x00ff_ffff;
const MASK64: u64 = 0xff00_0000_0000_0000;

pub struct Encoder {
    low: u64,
    range: u32,
    out: Vec<u8>,
}

impl Encoder {
    pub fn new() -> Self {
        Encoder { low: 0, range: u32::MAX, out: Vec::new() }
    }

    #[inline]
    fn encode(&mut self, sym_freq: u32, cum_freq: u32, total: u32) {
        self.range /= total;
        self.low = self.low.wrapping_add(self.range.wrapping_mul(cum_freq) as u64);
        self.range = self.range.wrapping_mul(sym_freq);
        while self.range <= TOP {
            if (self.low ^ self.low.wrapping_add(self.range as u64)) & MASK64 != 0 {
                let r = self.low as u32;
                self.range = (r | TOP) - r;
            }
            self.out.push((self.low >> 56) as u8);
            self.low <<= 8;
            self.range <<= 8;
        }
    }

    pub fn finish(mut self) -> Vec<u8> {
        for _ in 0..8 {
            self.out.push((self.low >> 56) as u8);
            self.low <<= 8;
        }
        self.out
    }
}

pub struct Decoder<'a> {
    buf: &'a [u8],
    pos: usize,
    low: u64,
    range: u32,
    code: u64,
}

impl<'a> Decoder<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        let mut d = Decoder { buf, pos: 0, low: 0, range: u32::MAX, code: 0 };
        for _ in 0..8 {
            d.code = (d.code << 8) | d.next_byte() as u64;
        }
        d
    }

    #[inline]
    fn next_byte(&mut self) -> u8 {
        let b = self.buf.get(self.pos).copied().unwrap_or(0);
        self.pos += 1;
        b
    }

    /// True if the decoder had to read past the end of its input (damaged data).
    pub fn overrun(&self) -> bool {
        self.pos > self.buf.len()
    }

    #[inline]
    fn cumulative_freq(&mut self, total: u32) -> u32 {
        self.range /= total;
        (self.code / self.range as u64) as u32
    }

    #[inline]
    fn update(&mut self, sym_freq: u32, low_end: u32) {
        let r = low_end.wrapping_mul(self.range);
        self.code = self.code.wrapping_sub(r as u64);
        self.low = self.low.wrapping_add(r as u64);
        self.range = self.range.wrapping_mul(sym_freq);
        while self.range <= TOP {
            if (self.low ^ self.low.wrapping_add(self.range as u64)) & MASK64 != 0 {
                let r = self.low as u32;
                self.range = (r | TOP) - r;
            }
            self.code = (self.code << 8).wrapping_add(self.next_byte() as u64);
            self.low <<= 8;
            self.range <<= 8;
        }
    }
}

/// Adaptive frequency model; counts start at 1 and are halved when the total reaches `max_total`.
pub struct Model {
    stats: Vec<u32>,
    total: u32,
    max_total: u32,
}

impl Model {
    pub fn new(n_symbols: usize, max_total: u32) -> Self {
        Model { stats: vec![1; n_symbols], total: n_symbols as u32, max_total }
    }

    #[inline]
    fn bump(&mut self, x: usize) {
        self.stats[x] += 1;
        self.total += 1;
        if self.total >= self.max_total {
            self.total = 0;
            for s in &mut self.stats {
                *s = s.div_ceil(2);
                self.total += *s;
            }
        }
    }

    #[inline]
    pub fn encode(&mut self, enc: &mut Encoder, x: usize) {
        let left: u32 = self.stats[..x].iter().sum();
        enc.encode(self.stats[x], left, self.total);
        self.bump(x);
    }

    #[inline]
    pub fn decode(&mut self, dec: &mut Decoder) -> Result<usize> {
        let cum = dec.cumulative_freq(self.total);
        let mut left = 0u32;
        for (x, &f) in self.stats.iter().enumerate() {
            if left + f > cum {
                dec.update(f, left);
                self.bump(x);
                return Ok(x);
            }
            left += f;
        }
        corrupt("range decoder out of sync")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let syms: Vec<usize> = (0..20000u32).map(|i| ((i * 7919) % 13 % 5) as usize).collect();
        let mut enc = Encoder::new();
        let mut m = Model::new(5, 256);
        for &s in &syms {
            m.encode(&mut enc, s);
        }
        let bytes = enc.finish();
        let mut dec = Decoder::new(&bytes);
        let mut m = Model::new(5, 256);
        for &s in &syms {
            assert_eq!(m.decode(&mut dec).unwrap(), s);
        }
        assert!(!dec.overrun());
    }
}
