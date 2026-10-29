//! CoMSA's variable-length integer: one length byte (0..=8) followed by that many
//! little-endian bytes.

use crate::error::{Result, corrupt};

pub fn put(v: &mut Vec<u8>, x: u64) {
    let n = (64 - x.leading_zeros()).div_ceil(8);
    v.push(n as u8);
    v.extend_from_slice(&x.to_le_bytes()[..n as usize]);
}

pub fn get(buf: &[u8], pos: &mut usize) -> Result<u64> {
    let Some(&n) = buf.get(*pos) else {
        return corrupt("truncated integer");
    };
    let n = n as usize;
    if n > 8 || *pos + 1 + n > buf.len() {
        return corrupt("bad integer encoding");
    }
    let mut bytes = [0u8; 8];
    bytes[..n].copy_from_slice(&buf[*pos + 1..*pos + 1 + n]);
    *pos += 1 + n;
    Ok(u64::from_le_bytes(bytes))
}

#[cfg(test)]
mod tests {
    #[test]
    fn roundtrip() {
        let vals = [0u64, 1, 255, 256, 65535, 1 << 40, u64::MAX];
        let mut v = Vec::new();
        for &x in &vals {
            super::put(&mut v, x);
        }
        assert_eq!(&v[..3], &[0, 1, 1]);
        let mut pos = 0;
        for &x in &vals {
            assert_eq!(super::get(&v, &mut pos).unwrap(), x);
        }
        assert_eq!(pos, v.len());
        assert!(super::get(&v, &mut pos).is_err());
    }
}
