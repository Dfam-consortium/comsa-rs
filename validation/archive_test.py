"""Whole-archive test: all accepted files concatenated (and via @list) into one archive.

usage: archive_test.py <files_dir> <work_dir> <results.json> <rust> <cpp_fixed> <cpp_orig>
"""
import os, sys, json, random, subprocess, time

files_dir, work, results, RS, CF, CO = sys.argv[1:7]

def normalize(data):
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

res = json.load(open(results))
accepted = sorted(r['file'] for r in res if 'rejected' not in r)
print(f'accepted files: {len(accepted)}')

big = os.path.join(work, 'all.sto')
lst = os.path.join(work, 'all.list')
norm_parts = {}
with open(big, 'wb') as out, open(lst, 'w') as l:
    for name in accepted:
        d = open(os.path.join(files_dir, name), 'rb').read()
        if d and not d.endswith(b'\n'):
            d += b'\n'
        out.write(d)
        l.write(os.path.join(files_dir, name) + '\n')
        norm_parts[name] = normalize(d)
expected = b''.join(norm_parts[n] for n in accepted)
raw_size = os.path.getsize(big)
print(f'concatenated input: {raw_size/1e6:.1f} MB')

def timed(args):
    t = time.time()
    r = subprocess.run(args, capture_output=True)
    return r.returncode, time.time() - t, r.stderr.decode()[-300:]

def check_dec(label, tool, archive, extra=()):
    out = os.path.join(work, 'dec.sto')
    rc, t, err = timed([tool, 'Sd', *extra, archive, out])
    same = rc == 0 and open(out, 'rb').read() == expected
    print(f'{label:45s} rc={rc} {t:7.2f}s  {"OK" if same else "MISMATCH " + err}')
    os.remove(out)

a_rs, a_cf, a_co, a_list = (os.path.join(work, x) for x in ('rs.smsac', 'cf.smsac', 'co.smsac', 'list.smsac'))
for label, args in [
    ('rust Sc (default threads)', [RS, 'Sc', big, a_rs]),
    ('fixed C++ Sc (default threads)', [CF, 'Sc', big, a_cf]),
    ('original C++ Sc', [CO, 'Sc', big, a_co]),
    ('rust Sc @list', [RS, 'Sc', '@' + lst, a_list]),
]:
    rc, t, err = timed(args)
    size = os.path.getsize(args[-1]) if os.path.exists(args[-1]) else -1
    print(f'{label:45s} rc={rc} {t:7.2f}s  size={size}  ratio={raw_size/size:.1f}' + (f'  {err}' if rc else ''))

rs_b = open(a_rs, 'rb').read()
print('rust archive == fixed C++ archive :', rs_b == open(a_cf, 'rb').read())
print('rust archive == rust @list archive:', rs_b == open(a_list, 'rb').read())
print('original C++ archive size vs rust :', os.path.getsize(a_co), len(rs_b))

check_dec('rust Sd of rust archive', RS, a_rs)
check_dec('rust Sd -t 1 of rust archive', RS, a_rs, ('-t', '1'))
check_dec('fixed C++ Sd of rust archive', CF, a_rs)
check_dec('original C++ Sd of rust archive', CO, a_rs)
check_dec('rust Sd of original C++ archive', RS, a_co)

# Listing and random access
r = subprocess.run([RS, 'Sl', a_rs], capture_output=True)
rows = r.stdout.decode().strip().split('\n')[1:]
print(f'Sl families: {len(rows)} (expected {len(accepted)})')
acs = [row.split('\t')[1].strip() for row in rows]
rng = random.Random(1)
bad = 0
for i in rng.sample(range(len(accepted)), 200):
    out = os.path.join(work, 'one.sto')
    subprocess.run([RS, 'Se', '-eAC', acs[i], a_rs, out], capture_output=True)
    if open(out, 'rb').read() != norm_parts[accepted[i]]:
        bad += 1
        print('  Se mismatch', accepted[i], acs[i])
print(f'Se random access: 200 checked, {bad} mismatches')

# General-purpose compressors for reference
for label, cmd in [('xz -9 -T16', ['xz', '-9', '-T16', '-c', big]), ('gzip -9', ['gzip', '-9', '-c', big])]:
    t = time.time()
    r = subprocess.run(cmd, capture_output=True)
    print(f'{label:45s} {time.time()-t:7.2f}s  size={len(r.stdout)}  ratio={raw_size/len(r.stdout):.1f}')
