//! Ordered parallel map: a producer feeds numbered items to a pool of workers and the
//! results are consumed in production order on the calling thread.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::Mutex;
use std::thread;

use crate::error::Result;

/// `produce` gets a `send` callback, which returns false once the pipeline is shutting
/// down (the producer should then stop). Memory is bounded by `capacity` items queued on
/// each side plus one item per worker.
pub fn ordered_map<T, U, P, W, C>(workers: usize, capacity: usize, produce: P, work: W, mut consume: C) -> Result<()>
where
    T: Send,
    U: Send,
    P: FnOnce(&mut dyn FnMut(T) -> bool) -> Result<()> + Send,
    W: Fn(T) -> U + Sync,
    C: FnMut(U) -> Result<()>,
{
    let (tx_in, rx_in) = sync_channel::<(usize, T)>(capacity.max(1));
    let (tx_out, rx_out) = sync_channel::<(usize, U)>(capacity.max(1));
    let rx_in = Mutex::new(rx_in);
    let stop = AtomicBool::new(false);

    thread::scope(|s| {
        let producer = s.spawn(|| {
            let mut next = 0usize;
            let mut send = |t: T| {
                if stop.load(Ordering::Relaxed) || tx_in.send((next, t)).is_err() {
                    return false;
                }
                next += 1;
                true
            };
            let r = produce(&mut send);
            drop(tx_in);
            r
        });

        for _ in 0..workers.max(1) {
            let tx = tx_out.clone();
            let (rx, work) = (&rx_in, &work);
            s.spawn(move || {
                loop {
                    let msg = rx.lock().unwrap().recv();
                    let Ok((i, t)) = msg else { break };
                    if tx.send((i, work(t))).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx_out);

        let mut pending = BTreeMap::new();
        let mut next = 0usize;
        let mut result = Ok(());
        for (i, u) in rx_out.iter() {
            pending.insert(i, u);
            while let Some(u) = pending.remove(&next) {
                next += 1;
                if result.is_ok() {
                    result = consume(u);
                    if result.is_err() {
                        stop.store(true, Ordering::Relaxed);
                    }
                }
            }
        }
        let produced = producer.join().expect("producer thread panicked");
        result.and(produced)
    })
}
