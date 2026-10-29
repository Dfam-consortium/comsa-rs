//! Zero-run-length coding (RLE-0). Runs of rank 0 are written as the bits of (len + 1),
//! least significant first and without the leading 1, using symbols 125 (bit 0) and 126 (bit 1).

pub const RUN_BIT0: u8 = 125;
pub const RUN_BIT1: u8 = 126;

pub fn encode(src: &[u8], out: &mut Vec<u8>) {
    out.clear();
    let mut run = 0u64;
    for &x in src {
        if x == 0 {
            run += 1;
        } else {
            if run > 0 {
                emit_run(out, run);
                run = 0;
            }
            out.push(x);
        }
    }
    if run > 0 {
        emit_run(out, run);
    }
}

#[inline]
fn emit_run(out: &mut Vec<u8>, len: u64) {
    let mut c = len + 1;
    while c != 1 {
        out.push(RUN_BIT0 + (c & 1) as u8);
        c >>= 1;
    }
}
