#!/usr/bin/env python3
"""Run indigo with primer-based trimming.

The trace is basecalled with tracy, the primer(s) and their reverse complements
are located in the primary basecalls (mismatch-tolerant, IUPAC aware), and the
left/right trim sizes are set (only hits in the outer third of the read count) so that everything up to and including a primer
near the 5' end, and from a primer near the 3' end onwards, is removed.
If no primer is found on a side, the fallback trim size for that side is used.

Example:
  indigo_primer.py sample.ab1 ref.fa -f ACGTACGTACGTACGTAC -r TTGCATGCATGCATGCA -o out
"""
import argparse
import os
import subprocess
import sys
import tempfile

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMP = str.maketrans('ACGTURYSWKMBDHVN', 'TGCAAYRSWMKVHDBN')


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def find_best(read, primer, max_mm):
    """Return (mismatches, start, end) of the best hit of primer in read, or None."""
    best = None
    n = len(primer)
    for i in range(len(read) - n + 1):
        mm = 0
        for j in range(n):
            if read[i + j] not in IUPAC.get(primer[j], primer[j]):
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm and (best is None or mm < best[0]):
            best = (mm, i, i + n)
    return best


def basecall(tracy, ab1, pratio):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'bc.fa')
        subprocess.run([tracy, 'basecall', '-p', str(pratio), '-f', 'fasta', '-o', out, ab1],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with open(out) as f:
            return ''.join(line.strip() for line in f if not line.startswith('>')).upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ab1', help='experiment chromatogram (.ab1/.abi/.scf)')
    ap.add_argument('genome', help='reference (.fa, .fa.gz or wildtype .ab1)')
    ap.add_argument('-f', '--fwd', help='forward primer (5\'->3\')')
    ap.add_argument('-r', '--rev', help='reverse primer (5\'->3\')')
    ap.add_argument('-o', '--out', default='out', help='output prefix [out]')
    ap.add_argument('-p', '--pratio', type=float, default=0.33, help='peak ratio [0.33]')
    ap.add_argument('-m', '--mismatches', type=int, default=2, help='max. primer mismatches [2]')
    ap.add_argument('--pad', type=int, default=0, help='extra bases to trim beyond a primer [0]')
    ap.add_argument('--ltrim', type=int, default=50, help='left trim if no primer found [50]')
    ap.add_argument('--rtrim', type=int, default=50, help='right trim if no primer found [50]')
    ap.add_argument('--tracy', default='tracy', help='tracy executable [tracy]')
    ap.add_argument('--dry-run', action='store_true', help='only report trim sizes')
    args = ap.parse_args()

    args.fwd = (args.fwd or '').upper().replace(' ', '') or None
    args.rev = (args.rev or '').upper().replace(' ', '') or None
    if not (args.fwd or args.rev):
        ap.error('give at least one primer with -f and/or -r')
    for p in (args.fwd, args.rev):
        if p and (len(p) < 10 or set(p) - set(IUPAC)):
            ap.error(f'primer "{p}" must be at least 10 nt of A/C/G/T/IUPAC codes')

    read = basecall(args.tracy, args.ab1, args.pratio)
    L = len(read)
    print(f'Read length: {L} bp', file=sys.stderr)

    left, right = [], []
    for name, primer in zip(['fwd', 'rev'], [args.fwd, args.rev]):
        if not primer:
            continue
        for orient, seq in (('+', primer), ('-', revcomp(primer))):
            hit = find_best(read, seq, args.mismatches)
            if hit is None:
                continue
            mm, s, e = hit
            mid = (s + e) / 2
            side = 'left' if mid < L / 3 else 'right' if mid > 2 * L / 3 else None
            print(f'  {name} primer ({orient}) at {s + 1}-{e}, {mm} mismatch(es) -> '
                  f'{side or "ignored (middle of the read)"}', file=sys.stderr)
            if side is None:
                continue
            if side == 'left':
                left.append(e + args.pad)
            else:
                right.append(L - s + args.pad)

    ltrim = max(left) if left else args.ltrim
    rtrim = max(right) if right else args.rtrim
    print(f'ltrim = {ltrim} ({"primer" if left else "fallback"}), '
          f'rtrim = {rtrim} ({"primer" if right else "fallback"})', file=sys.stderr)
    if ltrim + rtrim >= L - 20:
        sys.exit(f'Error: trimming leaves only {L - ltrim - rtrim} bp of {L} bp.')
    if args.dry_run:
        return

    indigo = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'indigo.sh')
    sys.exit(subprocess.call([indigo, args.ab1, args.genome, str(ltrim), str(rtrim), str(args.pratio), args.out]))


if __name__ == '__main__':
    main()
