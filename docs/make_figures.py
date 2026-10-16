# Figures for the docs: python docs/make_figures.py example/out.json docs/img  (needs matplotlib)
import json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = json.load(open(sys.argv[1]))
out = sys.argv[2]
col = {'A': '#4daf4a', 'C': '#377eb8', 'G': '#212121', 'T': '#e41a1c'}
bc = {int(k): v.split(':') for k, v in d['basecalls'].items()}
pos = d['pos']

# breakpoint: first base after which most calls are double
ks = sorted(bc)
bp = next(k for i, k in enumerate(ks) if sum('|' in bc[x][1] for x in ks[i:i + 15]) >= 9)

# --- Figure 1: chromatogram, clean -> double peaks
lo, hi = bp - 330, bp + 330
fig, ax = plt.subplots(figsize=(13, 4.2), dpi=150)
for b in 'ACGT':
    ys = d['peak' + b]
    ax.plot(pos[lo:hi], ys[lo:hi], color=col[b], lw=1.2, label=b)
top = max(max(d['peak' + b][lo:hi]) for b in 'ACGT')
for k in ks:
    if lo <= k < hi:
        bases = bc[k][1].split('|')
        for j, b in enumerate(bases):
            ax.text(k, top * (1.08 + 0.09 * j), b, color=col.get(b, 'k'), ha='center', fontsize=8, fontweight='bold')
ax.axvline(bp, color='#ff7f00', ls='--', lw=1.5)
ax.axvspan(lo, bp, color='#e8f5e9', alpha=0.6, zorder=0)
ax.axvspan(bp, hi, color='#fff3e0', alpha=0.6, zorder=0)
ax.text((lo + bp) / 2, top * 1.38, 'Clean: one peak per position', ha='center', fontsize=11)
ax.text((bp + hi) / 2, top * 1.38, 'After the indel: two overlapping alleles = double peaks', ha='center', fontsize=11)
ax.text(bp, -top * 0.12, 'break point (edit site)', color='#ff7f00', ha='center', fontsize=9)
ax.set_ylim(-top * 0.18, top * 1.48)
ax.set_xlim(lo, hi)
ax.set_xlabel('Trace signal position')
ax.set_ylabel('Fluorescence signal')
ax.legend(loc='upper left', ncol=4, frameon=False, bbox_to_anchor=(0, 0.93))
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f'{out}/fig1_chromatogram_double_peaks.png')

# --- Figure 2: decomposition
x, y = d['decomposition']['x'], d['decomposition']['y']
fig, ax = plt.subplots(figsize=(9, 3.8), dpi=150)
ax.plot(x, y, '-o', color='#377eb8', ms=4)
order = sorted(zip(y, x))
(e1, i1), (e2, i2) = order[0], order[1]
ax.annotate(f'{i1} bp: error {e1}\n(edited allele: {abs(i1)} bp deletion)', (i1, e1), xytext=(i1 - 12, e1 + 45),
            arrowprops=dict(arrowstyle='->', color='#e41a1c'), color='#e41a1c', fontsize=10)
ax.annotate(f'{i2} bp: error {e2}\n(wild-type allele)', (i2, e2), xytext=(i2 + 4, e2 + 45),
            arrowprops=dict(arrowstyle='->', color='#4daf4a'), color='#4daf4a', fontsize=10)
ax.set_xlabel('InDel length (bp)   negative = deletion, positive = insertion')
ax.set_ylabel('Decomposition error')
ax.set_title('Lowest errors at -7 and 0  ->  heterozygous 7 bp deletion', fontsize=11)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f'{out}/fig2_decomposition.png')
print('breakpoint signal pos', bp, 'base', bc[bp][0])
