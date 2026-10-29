//! Context-modelled range coding of the RLE-0 token stream, column by column.
//!
//! Each token is coded as a prefix (run bit 0, run bit 1, rank 1, other); "other" ranks are
//! further coded as a selector (bit length) and suffix (remaining bits).

use crate::error::{Result, corrupt};
use crate::rangecoder::{Decoder, Encoder, Model};
use crate::rle::{RUN_BIT0, RUN_BIT1};

/// Context sizes, chosen from the alignment size (rows * columns).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CtxLength {
    Tiny = 0,
    Small = 1,
    Medium = 2,
    Large = 3,
    Huge = 4,
}

// (prefix, selector, suffix) context counts
const CONTEXTS: [[usize; 3]; 5] = [[25, 8, 8], [125, 64, 8], [625, 64, 64], [3125, 64, 64], [3125, 512, 64]];

impl CtxLength {
    pub fn for_cells(cells: u64) -> Self {
        match cells {
            0..10_000 => CtxLength::Tiny,
            10_000..200_000 => CtxLength::Small,
            200_000..5_000_000 => CtxLength::Medium,
            5_000_000..20_000_000 => CtxLength::Large,
            _ => CtxLength::Huge,
        }
    }

    pub fn from_u8(x: u8) -> Option<Self> {
        Some(match x {
            0 => CtxLength::Tiny,
            1 => CtxLength::Small,
            2 => CtxLength::Medium,
            3 => CtxLength::Large,
            4 => CtxLength::Huge,
            _ => return None,
        })
    }
}

struct Models {
    prefix: Vec<Model>,
    selector: Vec<Model>,
    suffix: Vec<Model>,
}

impl Models {
    fn new(ctx: CtxLength) -> Self {
        let [np, ns, nsuf] = CONTEXTS[ctx as usize];
        Models {
            prefix: (0..np).map(|_| Model::new(4, 1 << 8)).collect(),
            selector: (0..ns).map(|_| Model::new(5, 1 << 8)).collect(),
            suffix: (0..nsuf).map(|i| Model::new(1 << (i % 8 + 1), 1 << 10)).collect(),
        }
    }
}

pub struct EntropyEncoder {
    enc: Encoder,
    m: Models,
    /// Total number of tokens coded (CoMSA's `pre_entropy_sequences_size`).
    pub n_tokens: u64,
}

impl EntropyEncoder {
    pub fn new(ctx: CtxLength) -> Self {
        EntropyEncoder { enc: Encoder::new(), m: Models::new(ctx), n_tokens: 0 }
    }

    /// Code the RLE-0 tokens of one column.
    pub fn encode_column(&mut self, tokens: &[u8]) {
        let (np, ns, nsuf) = (self.m.prefix.len(), self.m.selector.len(), self.m.suffix.len());
        let mut cp = np - 1;
        let mut cs = ns - 1;
        for &x in tokens {
            let prefix = match x {
                RUN_BIT0 => 0,
                RUN_BIT1 => 1,
                1 => 2,
                _ => 3,
            };
            self.m.prefix[cp].encode(&mut self.enc, prefix);
            cp = (cp * 5 + prefix) % np;
            if prefix < 3 {
                continue;
            }
            // WFC/MTF ranks of valid symbols are < 64, so the selector is 0..=4.
            assert!((2..64).contains(&x), "rank {x} out of range for the entropy coder");
            let bits = 8 - x.leading_zeros() as usize;
            let suffix = x as usize - (1 << (bits - 1));
            self.m.selector[cs].encode(&mut self.enc, bits - 2);
            cs = ((cs << 3) + bits - 2) % ns;
            self.m.suffix[cs % nsuf].encode(&mut self.enc, suffix);
        }
        self.n_tokens += tokens.len() as u64;
    }

    pub fn finish(self) -> Vec<u8> {
        self.enc.finish()
    }
}

pub struct EntropyDecoder<'a> {
    dec: Decoder<'a>,
    m: Models,
    pub n_tokens: u64,
}

impl<'a> EntropyDecoder<'a> {
    pub fn new(ctx: CtxLength, data: &'a [u8]) -> Self {
        EntropyDecoder { dec: Decoder::new(data), m: Models::new(ctx), n_tokens: 0 }
    }

    /// Decode one column of `n` ranks (zero runs already expanded) into `out`.
    pub fn decode_column(&mut self, n: usize, out: &mut Vec<u8>) -> Result<()> {
        let (np, ns, nsuf) = (self.m.prefix.len(), self.m.selector.len(), self.m.suffix.len());
        let mut cp = np - 1;
        let mut cs = ns - 1;
        let mut run_code = 0usize;
        let mut run_bits = 0u32;
        out.clear();

        while out.len() < n {
            let prefix = self.m.prefix[cp].decode(&mut self.dec)?;
            cp = (cp * 5 + prefix) % np;
            self.n_tokens += 1;

            if prefix > 1 && run_bits > 0 {
                // A literal ends the pending zero run.
                let len = run_code + (1 << run_bits) - 1;
                if out.len() + len >= n {
                    return corrupt("zero run overflows column");
                }
                out.resize(out.len() + len, 0);
                run_bits = 0;
            }

            match prefix {
                0 | 1 => {
                    if run_bits == 0 {
                        run_code = 0;
                    }
                    if run_bits >= 40 {
                        return corrupt("zero run too long");
                    }
                    run_bits += 1;
                    run_code += prefix << (run_bits - 1);
                    // A run that exactly completes the column ends it (no literal follows).
                    let len = run_code + (1 << run_bits) - 1;
                    if out.len() + len == n {
                        out.resize(n, 0);
                        run_bits = 0;
                    }
                }
                2 => out.push(1),
                _ => {
                    let sel = self.m.selector[cs].decode(&mut self.dec)?;
                    cs = ((cs << 3) + sel) % ns;
                    let suffix = self.m.suffix[cs % nsuf].decode(&mut self.dec)?;
                    out.push(((1usize << (sel + 1)) + suffix) as u8);
                }
            }
            if self.dec.overrun() {
                return corrupt("alignment stream truncated");
            }
        }
        Ok(())
    }
}
