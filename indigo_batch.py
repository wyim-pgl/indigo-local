#!/usr/bin/env python3
"""Run indigo on many chromatograms in parallel and summarize the results.

Input is either a folder of traces with one shared reference:
  indigo_batch.py -i traces/ -g ref.fa -o results/ -j 8

or a sample sheet (TSV or CSV) for per-sample references / primers / trims:
  indigo_batch.py -s samples.tsv -o results/ -j 8

Sample sheet columns (header required; only `ab1` and `reference` are mandatory):
  sample  ab1  reference  fwd  rev  ltrim  rtrim  pratio
Relative paths in the sheet are resolved against the sheet's folder.

Each sample is written to <outdir>/<sample>/<sample>.* with a log in
<outdir>/logs/<sample>.log. Samples that already finished successfully are
skipped, so an interrupted run can simply be restarted. At the end
<outdir>/summary.tsv and <outdir>/index.html are written (indigo_summary.py).
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.realpath(__file__))
TRACE_EXT = ('.ab1', '.abi', '.ab!', '.ab', '.scf')


def read_sheet(path):
    with open(path, newline='', encoding='utf-8-sig') as f:  # utf-8-sig drops Excel's BOM
        dialect = 'excel-tab' if '\t' in f.readline() else 'excel'
        f.seek(0)
        rows = list(csv.DictReader(f, dialect=dialect))
    base = os.path.dirname(os.path.abspath(path))
    samples = []
    for i, row in enumerate(rows, 1):
        row = {k.strip().lower(): (v or '').strip() for k, v in row.items() if k}
        if not row.get('ab1') or not row.get('reference'):
            sys.exit(f'{path}: row {i} needs "ab1" and "reference" columns')
        for key in ('ab1', 'reference'):
            if not os.path.isabs(row[key]):
                row[key] = os.path.join(base, row[key])
        row.setdefault('sample', '')
        if not row['sample']:
            row['sample'] = os.path.splitext(os.path.basename(row['ab1']))[0]
        samples.append(row)
    return samples


def read_folder(folder, reference):
    folder = os.path.abspath(folder)
    files = sorted(f for f in os.listdir(folder) if f.lower().endswith(TRACE_EXT))
    return [{'sample': os.path.splitext(f)[0], 'ab1': os.path.join(folder, f),
             'reference': os.path.abspath(reference)} for f in files]


def safe_name(name):
    name = re.sub(r'[^\w.-]', '_', name)
    return '_' + name if name.startswith('.') or not name else name


def build_cmd(s, args, prefix):
    ltrim = s.get('ltrim') or str(args.ltrim)
    rtrim = s.get('rtrim') or str(args.rtrim)
    pratio = s.get('pratio') or str(args.pratio)
    fwd = s.get('fwd') or args.fwd
    rev = s.get('rev') or args.rev
    if fwd or rev:
        cmd = [sys.executable, os.path.join(HERE, 'indigo_primer.py'), s['ab1'], s['reference'],
               '-o', prefix, '-p', pratio, '--ltrim', ltrim, '--rtrim', rtrim, '-m', str(args.mismatches)]
        if fwd:
            cmd += ['-f', fwd]
        if rev:
            cmd += ['-r', rev]
        return cmd
    return [os.path.join(HERE, 'indigo.sh'), s['ab1'], s['reference'], ltrim, rtrim, pratio, prefix]


def run_one(s, args):
    name = s['sample']
    outdir = os.path.join(args.outdir, name)
    prefix = os.path.join(outdir, name)
    done = prefix + '.done'  # written only after indigo finished successfully
    if os.path.exists(done) and not args.force:
        return name, 'skipped'
    try:
        for key in ('ab1', 'reference'):
            if not os.path.exists(s[key]):
                return name, f'missing {key}: {s[key]}'
        os.makedirs(outdir, exist_ok=True)
        for ext in ('.done', '.json', '.html'):  # no stale results from an earlier run
            if os.path.exists(prefix + ext):
                os.remove(prefix + ext)
        log = os.path.join(args.outdir, 'logs', name + '.log')
        with open(log, 'w') as fh:
            rc = subprocess.call(build_cmd(s, args, prefix), stdout=fh, stderr=subprocess.STDOUT, cwd=outdir)
        if rc != 0:
            return name, f'failed (exit {rc}, see logs/{name}.log)'
        open(done, 'w').close()
        return name, 'ok'
    except Exception as e:
        return name, f'failed ({e})'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('-i', '--input', help='folder with .ab1/.abi/.scf files')
    src.add_argument('-s', '--sheet', help='sample sheet (TSV/CSV)')
    ap.add_argument('-g', '--genome', help='reference for all files (with -i)')
    ap.add_argument('-o', '--outdir', default='results', help='output folder [results]')
    ap.add_argument('-j', '--jobs', type=int, default=4, help='parallel jobs [4]')
    ap.add_argument('--ltrim', type=int, default=50, help='default left trim [50]')
    ap.add_argument('--rtrim', type=int, default=50, help='default right trim [50]')
    ap.add_argument('-p', '--pratio', type=float, default=0.33, help='default peak ratio [0.33]')
    ap.add_argument('-f', '--fwd', help='default forward primer (enables primer trimming)')
    ap.add_argument('-r', '--rev', help='default reverse primer (enables primer trimming)')
    ap.add_argument('-m', '--mismatches', type=int, default=2, help='max. primer mismatches [2]')
    ap.add_argument('--force', action='store_true', help='re-run samples that already have results')
    args = ap.parse_args()

    if args.input:
        if not args.genome:
            ap.error('-i needs -g <reference>')
        samples = read_folder(args.input, args.genome)
    else:
        samples = read_sheet(args.sheet)
    if not samples:
        sys.exit('No samples found.')
    for s in samples:
        s['sample'] = safe_name(s['sample'])
    names = [s['sample'] for s in samples]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        sys.exit(f'Duplicate sample names: {", ".join(dups)}')

    args.outdir = os.path.abspath(args.outdir)
    os.makedirs(os.path.join(args.outdir, 'logs'), exist_ok=True)
    print(f'{len(samples)} samples, {args.jobs} parallel jobs -> {args.outdir}', file=sys.stderr)

    status = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, s, args) for s in samples]
        for n, fut in enumerate(as_completed(futures), 1):
            name, st = fut.result()
            status[name] = st
            print(f'[{n}/{len(samples)}] {name}: {st}', file=sys.stderr)

    with open(os.path.join(args.outdir, 'status.tsv'), 'w') as f:
        f.write('sample\tstatus\n')
        for name in names:
            f.write(f'{name}\t{status[name]}\n')

    subprocess.call([sys.executable, os.path.join(HERE, 'indigo_summary.py'), args.outdir])
    failed = [n for n in names if not status[n] in ('ok', 'skipped')]
    print(f'Done: {len(names) - len(failed)} ok/skipped, {len(failed)} failed.', file=sys.stderr)
    if failed:
        print('Failed: ' + ', '.join(failed[:20]) + (' ...' if len(failed) > 20 else ''), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
