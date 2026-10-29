"""Per-file round-trip verification of Dfam Stockholm files with the Rust and C++ CoMSA tools.

usage: verify.py <files_dir> <work_dir> <rust_bin> <cpp_fixed_bin> <cpp_orig_bin> <jobs>
"""
import os, sys, subprocess, json, collections, random
from concurrent.futures import ThreadPoolExecutor

files_dir, work, RS, CF, CO, jobs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])
os.makedirs(work, exist_ok=True)

VALID = set(b'-.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*~')


def normalize(data):
    """Expected decompressed text: non-empty lines, LF endings, every family ends with //."""
    out, fam = [], []
    for line in data.replace(b'\r', b'\n').split(b'\n'):
        if not line:
            continue
        if line == b'//':
            if fam:
                out += fam + [b'//']
            fam = []
        else:
            fam.append(line)
    if fam:
        out += fam + [b'//']
    return b''.join(l + b'\n' for l in out)


def to_fasta(norm):
    """First family as aligned FASTA (name, sequence) or None."""
    recs = []
    for line in norm.split(b'\n'):
        if line == b'//':
            break
        if line and not line.startswith(b'#'):
            parts = line.split(None, 1)
            if len(parts) == 2:
                recs.append((parts[0], parts[1]))
    return b''.join(b'>' + n + b'\n' + s + b'\n' for n, s in recs) if recs else None


def run(args, timeout=600):
    try:
        r = subprocess.run(args, capture_output=True, timeout=timeout)
        return r.returncode, r.stderr.decode('utf-8', 'replace').strip().replace('\n', ' ')[-300:]
    except subprocess.TimeoutExpired:
        return 'TIMEOUT', ''


def read(p):
    try:
        with open(p, 'rb') as f:
            return f.read()
    except FileNotFoundError:
        return None


def check(name):
    src = os.path.join(files_dir, name)
    base = os.path.join(work, name[:-4])
    raw = read(src)
    norm = normalize(raw)
    r = {'file': name, 'size': len(raw), 'issues': []}
    bad = r['issues'].append

    rc_rs, err_rs = run([RS, 'Sc', '-t', '1', src, base + '.rs.smsac'])
    rc_cf, err_cf = run([CF, 'Sc', '-t', '1', src, base + '.cf.smsac'])
    r['rs_rc'], r['cf_rc'] = rc_rs, rc_cf
    if rc_rs != 0 or rc_cf != 0:
        r['rejected'] = err_rs or err_cf
        if (rc_rs == 0) != (rc_cf == 0):
            bad(f'rust/C++ disagree on accepting input: rust={rc_rs} {err_rs} | c++={rc_cf} {err_cf}')
        for p in (base + '.rs.smsac', base + '.cf.smsac'):
            if os.path.exists(p):
                bad(f'output left behind after failure: {p}')
                os.remove(p)
        return r

    a_rs, a_cf = read(base + '.rs.smsac'), read(base + '.cf.smsac')
    r['comp_size'] = len(a_rs)
    if a_rs != a_cf:
        bad('rust and fixed C++ archives differ')

    def decode_check(tool, label, archive):
        out = base + '.' + label + '.dec'
        # the original binary has no -t option (and exits 0 on bad options)
        args = [tool, 'Sd', archive, out] if tool == CO else [tool, 'Sd', '-t', '1', archive, out]
        rc, err = run(args)
        if rc != 0:
            bad(f'{label}: decode failed rc={rc} {err}')
        elif read(out) != norm:
            bad(f'{label}: decoded text differs from normalized input')
        if os.path.exists(out):
            os.remove(out)

    decode_check(RS, 'rust_sd', base + '.rs.smsac')
    decode_check(CF, 'cppfixed_sd_of_rust', base + '.rs.smsac')
    if b'~' not in norm:
        decode_check(CO, 'cpporig_sd_of_rust', base + '.rs.smsac')

    # Fast variant
    rc, err = run([RS, 'Sc', '-f', '-t', '1', src, base + '.f.smsac'])
    if rc != 0:
        bad(f'fast compress failed {err}')
    else:
        decode_check(RS, 'rust_fast_sd', base + '.f.smsac')
        rc2, _ = run([CF, 'Sc', '-f', '-t', '1', src, base + '.cff.smsac'])
        if rc2 == 0 and read(base + '.f.smsac') != read(base + '.cff.smsac'):
            bad('fast: rust and fixed C++ archives differ')

    # Random access / listing
    rc, err = run([RS, 'Sl', base + '.rs.smsac'])
    if rc != 0:
        bad(f'Sl failed {err}')

    # FASTA path on the first family
    fa = to_fasta(norm)
    if fa is not None:
        with open(base + '.fa', 'wb') as f:
            f.write(fa)
        rc, err = run([RS, 'Fc', '-t', '1', base + '.fa', base + '.msac'])
        if rc != 0:
            bad(f'Fc failed {err}')
        else:
            for tool, label in ((RS, 'rust_fd'), (CF, 'cppfixed_fd')):
                rc, err = run([tool, 'Fd', base + '.msac', base + '.fa.dec'])
                if rc != 0 or read(base + '.fa.dec') != fa:
                    bad(f'{label}: FASTA round trip failed rc={rc} {err}')

    for ext in ('.rs.smsac', '.cf.smsac', '.f.smsac', '.cff.smsac', '.fa', '.msac', '.fa.dec'):
        if os.path.exists(base + ext):
            os.remove(base + ext)
    return r


names = sorted(n for n in os.listdir(files_dir) if n.endswith('.sto'))
with ThreadPoolExecutor(jobs) as ex:
    results = list(ex.map(check, names))

with open(os.path.join(work, 'results.json'), 'w') as f:
    json.dump(results, f)

n = len(results)
rejected = [r for r in results if 'rejected' in r]
issues = [r for r in results if r['issues']]
ok = [r for r in results if 'rejected' not in r and not r['issues']]
print(f'files: {n}  passed: {len(ok)}  rejected (both tools): {len(rejected) - len([r for r in rejected if r["issues"]])}  with issues: {len(issues)}')
tot_raw = sum(r['size'] for r in ok)
tot_c = sum(r['comp_size'] for r in ok)
if tot_c:
    print(f'passed files: {tot_raw/1e6:.1f} MB -> {tot_c/1e6:.2f} MB (per-file archives), ratio {tot_raw/tot_c:.2f}')
reasons = collections.Counter(r['rejected'].split(': ', 2)[-1][:110] for r in rejected)
for k, v in reasons.most_common(15):
    print(f'  rejected x{v}: {k}')
for r in issues[:30]:
    print('ISSUE', r['file'], r['issues'][:3])
