"""Download a random sample of Dfam seed alignments (Stockholm) for round-trip testing.

usage: fetch.py <out_dir> <n_total> <n_curated>
Samples curated (DF) accessions from 1..31000 and uncurated (DR) accessions from
1..7650000, skipping accessions that return 404, until n_total files are saved.
Six requests at a time; a manifest is appended to <out_dir>/../manifest.tsv.
"""
import os, random, sys, time, threading
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

out_dir = sys.argv[1]
target = int(sys.argv[2])
n_df = int(sys.argv[3])          # how many curated DF families among the target
URL = 'https://www.dfam.org/api/families/{}/seed?format=stockholm'

rng = random.Random(20261029)
lock = threading.Lock()
ok = {'DF': 0, 'DR': 0}
log = open(os.path.join(out_dir, '..', 'manifest.tsv'), 'a')

def fetch(acc):
    path = os.path.join(out_dir, acc + '.sto')
    if os.path.exists(path):
        return acc, 200, os.path.getsize(path)
    for attempt in range(6):
        try:
            req = urllib.request.Request(URL.format(acc), headers={'User-Agent': 'comsa-rs validation (https://github.com/Dfam-consortium/comsa-rs)'})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            with open(path + '.part', 'wb') as f:
                f.write(data)
            os.rename(path + '.part', path)
            return acc, 200, len(data)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return acc, 404, 0
            if e.code not in (429, 500, 502, 503, 504):
                return acc, e.code, 0
        except Exception:
            pass
        time.sleep(2 ** attempt)
    return acc, -1, 0

def candidates(prefix, hi, want):
    seen = set()
    while True:
        with lock:
            if ok[prefix] >= want:
                return
        n = rng.randint(1, hi)
        if n in seen:
            continue
        seen.add(n)
        yield f'{prefix}{n:09d}'

def run(prefix, hi, want):
    with ThreadPoolExecutor(6) as ex:
        it = candidates(prefix, hi, want)
        pending = set()
        while True:
            while len(pending) < 12:
                try:
                    pending.add(ex.submit(fetch, next(it)))
                except StopIteration:
                    break
            if not pending:
                return
            finished, _ = wait(pending, return_when=FIRST_COMPLETED)
            for done in finished:
                pending.discard(done)
                record(prefix, *done.result())

def record(prefix, acc, code, size):
    log.write(f'{acc}\t{code}\t{size}\n'); log.flush()
    if code == 200:
        with lock:
            ok[prefix] += 1
            if (ok['DF'] + ok['DR']) % 500 == 0:
                print(ok, flush=True)

run('DF', 31000, n_df)
run('DR', 7650000, target - n_df)
print('done', ok, flush=True)
