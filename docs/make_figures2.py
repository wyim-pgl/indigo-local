# Teaching figures: python docs/make_figures2.py example/out.json docs/img  (needs matplotlib)
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

d = json.load(open(sys.argv[1]))
out = sys.argv[2]
COL = {'A': '#4daf4a', 'C': '#377eb8', 'G': '#212121', 'T': '#e41a1c'}
bc = {int(k): v.split(':') for k, v in d['basecalls'].items()}
ks = sorted(bc)


def clean(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


# ---------- fig6: base quality along the read and the trimmed ends
q = d['basecallQual']
n = len(q)
lt, rt = d['meta']['arguments']['trimLeft'], d['meta']['arguments']['trimRight']
fig, ax = plt.subplots(figsize=(11, 3.6), dpi=150)
ax.bar(range(1, n + 1), q, width=1.0, color=['#adb5bd' if (i < lt or i >= n - rt) else '#4c6ef5' for i in range(n)])
ax.axvspan(0, lt + 0.5, color='#ffe3e3', zorder=0)
ax.axvspan(n - rt + 0.5, n + 1, color='#ffe3e3', zorder=0)
ax.axhline(20, color='#e8590c', ls='--', lw=1)
ax.text(n + 2, 20, 'Q20 (99% accurate)', color='#e8590c', va='center', fontsize=9)
ax.text(lt / 2, max(q) * 1.05, f'left trim\n{lt} bases', ha='center', fontsize=10, color='#c92a2a')
ax.text(n - rt / 2, max(q) * 1.05, f'right trim\n{rt} bases', ha='center', fontsize=10, color='#c92a2a')
ax.text(n / 2, max(q) * 1.05, 'analysed region', ha='center', fontsize=11, color='#364fc7')
ax.set_xlim(0, n + 1)
ax.set_ylim(0, max(q) * 1.25)
ax.set_xlabel('Base position in the read')
ax.set_ylabel('Base quality (Phred)')
ax.set_title('Sanger reads are noisy at both ends - trimming removes them before the analysis', fontsize=11)
clean(ax)
fig.tight_layout()
fig.savefig(f'{out}/fig6_trim_quality.png')

# ---------- fig7: peak ratio at a single position
pos = d['pos']
peak = {b: d['peak' + b] for b in 'ACGT'}
# a clean and a double position close together
dbl = [k for k in ks if '|' in bc[k][1] and int(bc[k][0]) > 300]
sgl = [k for k in ks if '|' not in bc[k][1] and int(bc[k][0]) < 230 and int(bc[k][0]) > 150]
picks = [('single peak', sgl[5]), ('double peak', dbl[3])]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), dpi=150, sharey=False)
for ax, (title, k) in zip(axes, picks):
    lo, hi = k - 18, k + 18
    for b in 'ACGT':
        ax.plot(pos[lo:hi], peak[b][lo:hi], color=COL[b], lw=2, label=b)
    hs = sorted(((peak[b][k - 1], b) for b in 'ACGT'), reverse=True)
    (h1, b1), (h2, b2) = hs[0], hs[1]
    thr = 0.33 * h1
    ax.axhline(thr, color='#e8590c', ls='--', lw=1.2)
    ax.text(hi - 0.5, thr, ' 0.33 x highest', color='#e8590c', fontsize=9, ha='left', va='center')
    ax.axvline(k, color='#999', lw=0.8, ls=':')
    ratio = h2 / h1 if h1 else 0
    call = f'{b1}' if ratio < 0.33 else f'{b1} | {b2}'
    ax.set_title(f'{title}: 2nd/1st = {h2}/{h1} = {ratio:.2f}  ->  call {call}', fontsize=10)
    ax.set_xlabel('Trace signal position')
    clean(ax)
axes[0].set_ylabel('Signal')
axes[0].legend(ncol=4, frameon=False, fontsize=9, loc='upper right')
fig.suptitle('How pratio decides whether a second base is called (secondary peak >= pratio x primary peak)', fontsize=11)
fig.tight_layout()
fig.savefig(f'{out}/fig7_peak_ratio.png')


# ---------- fig8: from Cas9 cut to double peaks (schematic)
fig, ax = plt.subplots(figsize=(12, 6.2), dpi=150)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')
ref = 'CATCCAGGCTTCCCTCTGGAGGGTCCAGG'
x0, w = 6, 2.4


def seq_row(y, seq, label, colors=None, gaps=()):
    ax.text(x0 - 1, y, label, ha='right', va='center', fontsize=11, fontweight='bold')
    for i, c in enumerate(seq):
        fc = colors[i] if colors else 'white'
        ax.add_patch(Rectangle((x0 + i * w, y - 2.2), w * 0.95, 4.4, fc=fc, ec='#dee2e6', lw=0.5))
        ax.text(x0 + i * w + w * 0.47, y, c, ha='center', va='center', fontsize=10, family='monospace',
                color='white' if fc == '#e03131' else 'black')


cols = ['#cfe2ff'] * 20 + ['#ffc078'] * 3 + ['white'] * 6
seq_row(88, ref, 'Reference', cols)
cutx = x0 + 17 * w
ax.plot([cutx, cutx], [83, 93], color='#e03131', lw=2.5)
ax.text(x0 + 10 * w, 95, 'guide RNA (20 nt)', ha='center', color='#1c7ed6', fontsize=10)
ax.text(x0 + 21.5 * w, 95, 'PAM', ha='center', color='#e8590c', fontsize=10, fontweight='bold')
ax.text(cutx, 79.5, 'Cas9 cuts 3 bp before the PAM', ha='center', color='#e03131', fontsize=9)
ax.annotate('', xy=(50, 66), xytext=(50, 76), arrowprops=dict(arrowstyle='->', lw=2))
ax.text(52, 71, 'cell repairs the cut (NHEJ) - often small deletions/insertions', fontsize=10, va='center')
a2 = ref[:17] + '-' * 7 + ref[24:]
seq_row(58, ref, 'Allele 1')
seq_row(50, a2, 'Allele 2', ['#e03131' if c == '-' else 'white' for c in a2])
ax.text(x0 + 30 * w, 58, 'wild type', va='center', fontsize=10, color='#2b8a3e')
ax.text(x0 + 30 * w, 50, '7 bp deletion', va='center', fontsize=10, color='#c92a2a')
ax.annotate('', xy=(50, 33), xytext=(50, 44), arrowprops=dict(arrowstyle='->', lw=2))
ax.text(52, 39, 'PCR amplifies BOTH alleles together -> one Sanger read', fontsize=10, va='center')
# mini traces
import math
seqA, seqB = ref, ref[:17] + ref[24:] + 'TTTCC'
xs = [i / 4 for i in range(4 * 29)]
for b in 'ACGT':
    ys = []
    for x in xs:
        v = 0
        for i in range(29):
            ga = math.exp(-((x - i) ** 2) / 0.08)
            if i < len(seqA) and seqA[i] == b:
                v += ga * (1.0 if i < 17 else 0.5)
            if 17 <= i < len(seqB) and seqB[i] == b:
                v += ga * 0.5
        ys.append(v)
    ax.plot([x0 + x * w + w * 0.47 for x in xs], [6 + y * 18 for y in ys], color=COL[b], lw=1.5)
ax.plot([cutx, cutx], [4, 27], color='#e03131', lw=1.5, ls='--')
ax.text(x0 + 8 * w, 28, 'one peak per position', ha='center', fontsize=10, color='#2b8a3e')
ax.text(x0 + 23 * w, 28, 'two peaks per position', ha='center', fontsize=10, color='#c92a2a')
ax.text(x0 - 1, 15, 'Chromatogram', ha='right', va='center', fontsize=11, fontweight='bold')
fig.savefig(f'{out}/fig8_cut_to_double_peaks.png', bbox_inches='tight')


# ---------- fig9: what different genotypes look like (simulated)
def trace(alleles, n=26):
    xs = [i / 5 for i in range(5 * n)]
    res = {}
    for b in 'ACGT':
        res[b] = [sum(f * math.exp(-((x - i) ** 2) / 0.07) for s, f in alleles for i in range(min(n, len(s))) if s[i] == b)
                  for x in xs]
    return xs, res


base = 'GCTTCCCTCTGGAGGGTCCAGGTTTCCAGATATATGC'
cut = 11
dele = base[:cut] + base[cut + 7:]
ins = base[:cut] + 'A' + base[cut:]
cases = [('Wild type (unedited)', [(base, 1.0)], 'one peak per position everywhere'),
         ('Heterozygous: WT + 7 bp deletion', [(base, 0.5), (dele, 0.5)], 'double peaks from the cut site on'),
         ('Homozygous: 7 bp deletion on both alleles', [(dele, 1.0)], 'single peaks again - the deletion is only seen in the alignment to the reference'),
         ('Biallelic: -7 bp and +1 bp', [(dele, 0.5), (ins, 0.5)], 'double peaks from the cut site on, no wild-type allele')]
fig, axes = plt.subplots(4, 1, figsize=(11, 9), dpi=150, sharex=True)
for ax, (title, al, note) in zip(axes, cases):
    xs, res = trace(al)
    for b in 'ACGT':
        ax.plot(xs, res[b], color=COL[b], lw=1.4)
    ax.axvline(cut - 0.5, color='#e03131', ls='--', lw=1.2)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_title(note, loc='right', fontsize=9, color='#495057')
    ax.set_yticks([])
    ax.set_ylim(0, 1.15)
    clean(ax)
axes[0].text(cut - 0.4, 1.0, 'cut site', color='#e03131', fontsize=9)
axes[-1].set_xlabel('Position (bases)   -   simulated traces for illustration, not real data')
fig.tight_layout()
fig.savefig(f'{out}/fig9_genotype_patterns.png')
print('ok')
