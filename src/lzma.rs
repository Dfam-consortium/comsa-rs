//! LZMA (.xz container, CRC64) for sequence names and annotation text.

use liblzma::stream::{Action, CONCATENATED, Check, Filters, LzmaOptions, PRESET_EXTREME, Status, Stream};

use crate::error::{Result, corrupt};

/// Dictionary size of presets 7-9; larger than any useful dictionary for a single family.
const MAX_DICT: u32 = 64 << 20;

pub fn compress(input: &[u8], preset: u32, extreme: bool) -> Result<Vec<u8>> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut opts = LzmaOptions::new_preset(preset | if extreme { PRESET_EXTREME } else { 0 })?;
    // A dictionary larger than the input gains nothing but costs a large, zeroed allocation.
    let mut dict = 4096u32;
    while (dict as usize) < input.len() && dict < MAX_DICT {
        dict <<= 1;
    }
    opts.dict_size(dict);
    let mut filters = Filters::new();
    filters.lzma2(&opts);
    let mut stream = Stream::new_stream_encoder(&filters, Check::Crc64)?;

    let mut out = Vec::with_capacity(input.len() / 4 + 128);
    let mut consumed = 0;
    loop {
        if out.capacity() == out.len() {
            out.reserve(out.len().max(4096));
        }
        let before = stream.total_in();
        let status = stream.process_vec(&input[consumed..], &mut out, Action::Finish)?;
        consumed += (stream.total_in() - before) as usize;
        if status == Status::StreamEnd {
            return Ok(out);
        }
    }
}

pub fn decompress(input: &[u8]) -> Result<Vec<u8>> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut stream = Stream::new_stream_decoder(u64::MAX, CONCATENATED)?;
    let mut out = Vec::with_capacity(input.len() * 4);
    let mut consumed = 0;
    loop {
        if out.capacity() == out.len() {
            out.reserve(out.len().max(4096));
        }
        let (before_in, before_out) = (stream.total_in(), stream.total_out());
        let status = stream.process_vec(&input[consumed..], &mut out, Action::Finish)?;
        consumed += (stream.total_in() - before_in) as usize;
        if status == Status::StreamEnd {
            return Ok(out);
        }
        if consumed == input.len() && stream.total_out() == before_out && out.len() < out.capacity() {
            return corrupt("truncated LZMA stream");
        }
    }
}
