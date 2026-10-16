#!/usr/bin/env python3
"""Render a tracy/indigo JSON result as a standalone HTML report (like the Indigo web app).

Usage: indigo_report.py <prefix>.json [output.html] [--grna SEQ] [--pam NGG] [--cut-offset -3]

With --grna (20 nt protospacer, 5'->3', with or without PAM) the report gets a
"CRISPR edit" panel: guide, PAM and Cas9 cut site highlighted on the reference,
each allele aligned underneath with deleted / inserted bases marked, and the
InDel length and distance to the cut site per allele. Without --grna the panel
still lists the InDels of each allele, centred on the first one.
The panel data is also written to <prefix>.crispr.json (used by indigo_summary.py).
"""
import argparse
import html
import json
import os
import sys

TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Indigo report - __TITLE__</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.6.2/dist/css/bootstrap.min.css">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body { padding: 20px; }
  pre { font-size: 12px; background: #f8f9fa; padding: 10px; }
  h4 { margin-top: 30px; }
  td.show { cursor: pointer; color: #007bff; }
  .help { font-size: 13px; color: #555; background: #f6f8fa; border-left: 3px solid #adb5bd; padding: 6px 10px; margin: 6px 0 12px; }
  .help b { color: #333; }
  .edit { font-family: SFMono-Regular, Menlo, Consolas, monospace; font-size: 14px; line-height: 1.9; overflow-x: auto; white-space: nowrap; background: #fbfbfb; border: 1px solid #eee; padding: 8px; }
  .edit .lab { display: inline-block; width: 44px; color: #666; }
  .edit .num { display: inline-block; width: 52px; color: #999; text-align: right; padding-right: 6px; }
  .edit span.b { display: inline-block; width: 0.7em; text-align: center; }
  .g { background: #cfe2ff; } .pam { background: #ffc078; font-weight: bold; }
  .del { background: #e03131; color: #fff; } .ins { background: #2f9e44; color: #fff; }
  .mm { background: #ffe066; } .nc { color: #ccc; }
  .cut { display: inline-block; width: 0; border-left: 2px solid #e03131; height: 1.3em; vertical-align: middle; margin: 0 1px; }
  .alleles td, .alleles th { padding: 3px 10px; }
</style>
</head>
<body>
<div class="container-fluid">
  <h2>Indigo results</h2>
  <p class="text-muted" id="meta"></p>

__CRISPR__

  <h4 id="sec-variants">Variants</h4>
  <div class="help">All differences of the two separated alleles to the reference (VCF style, left-aligned).
  <b>type</b>: SNV / Deletion / Insertion. <b>ref &rarr; alt</b> <code>GGAGGGTC &rarr; G</code> = the 7 bases after the first G are deleted.
  <b>genotype</b> het. = only one allele carries it. <b>filter</b> LowQual calls are usually noise.
  Inside the double-peak region many SNV-like calls can appear; focus on the Deletion/Insertion rows and the decomposition.
  Click &#128200; to jump to the position in the trace.</div>
  <div id="variants-table"></div>

  <h4 id="sec-trace">Trace</h4>
  <div class="help">The raw chromatogram. Lines = fluorescence of A (green), C (blue), G (black), T (red); dots above = base calls.
  One peak per position = a single sequence. From the edit site onward two alleles overlap and you see <b>two peaks per position</b>.
  Drag to zoom, double-click to reset.</div>
  <div id="trace-chart" style="height:500px"></div>

  <h4 id="sec-decomp">Decomposition</h4>
  <div class="help">For every InDel length (x; negative = deletion, positive = insertion) Indigo tests how well
  "allele 2 = allele 1 shifted by this length" explains the double peaks (y = error, lower is better).
  <b>Minima at 0 and &minus;n</b>: wild type + n bp deletion (heterozygous). <b>Two non-zero minima</b>: both alleles edited (biallelic).
  <b>Only 0 low</b>: no InDel. <b>No clear minimum</b>: mosaic or poor quality.</div>
  <div id="decomposition-chart" style="height:400px"></div>

  <h4 id="sec-align1">Alignment of Alt1 to reference</h4>
  <div class="help">Allele 1 (from the primary peaks) aligned to the reference. <code>|</code> = match; a gap <code>-</code> in Alt = deleted bases,
  a gap in Ref = inserted bases. Gaps at the very start/end are only the trimmed ends.
  "Estimated allelic fraction" = share of this allele in the trace. <code>_reverse</code> = the read maps to the reverse strand.</div>
  <div id="alignment-chart-1"></div>

  <h4 id="sec-align2">Alignment of Alt2 to reference</h4>
  <div class="help">Allele 2 (reconstructed from the secondary peaks) aligned to the reference.</div>
  <div id="alignment-chart-2"></div>

  <h4 id="sec-align3">Alignment of Alt1 to Alt2</h4>
  <div class="help">The two alleles aligned to each other; their difference is the edit.</div>
  <div id="alignment-chart-3"></div>
</div>

<script>
const data = __DATA__;

function zip() {
  const ret = []
  for (let i = 0; i < arguments[0].length; i += 1) {
    const record = [arguments[0][i]]
    for (let j = 1; j < arguments.length; j += 1) record.push(arguments[j][i])
    ret.push(record)
  }
  return ret
}
const isDna = seq => /^[acgt]+$/i.test(seq)
const ungapped = seq => seq.replace(/-/g, '')
const escapeHtml = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]))

// Meta
if (data.meta) {
  const a = data.meta.arguments || {}
  document.getElementById('meta').textContent =
    `${data.meta.program} ${data.meta.version} | input: ${a.input} | reference: ${a.genome} | ` +
    `trimLeft: ${a.trimLeft} | trimRight: ${a.trimRight} | pratio: ${a.pratio}`
}

// Trace
function convertTraceData(data) {
  const ret = []
  const baseCallPat = /(\d+):([A-Z|]+)/
  for (const [pos, A, C, G, T] of zip(data.pos, data.peakA, data.peakC, data.peakG, data.peakT)) {
    const record = { position: pos, peaks: { A, C, G, T }, calls: null }
    const baseCall = data.basecalls[pos]
    if (baseCall) {
      const match = baseCallPat.exec(baseCall)
      if (match) {
        const [, p, bases] = match
        if (!isDna(bases.replace(/\|/g, ''))) continue
        record.calls = { pos: +p, bases: bases.split('|') }
      }
    }
    ret.push(record)
  }
  return ret
}

const traceChart = document.getElementById('trace-chart')
;(function () {
  const recs = convertTraceData(data)
  const colors = { A: '#4daf4a', C: '#377eb8', G: '#212121', T: '#e41a1c' }
  const traces = [], calls = []
  for (const base of ['A', 'C', 'G', 'T']) {
    calls.push({ x: [], y: [], xaxis: 'x', yaxis: 'y2', name: base, mode: 'markers',
                 hoverinfo: 'x+text', text: [], marker: { color: colors[base], size: 10 } })
    traces.push({ x: recs.map(r => r.position), y: recs.map(r => r.peaks[base]),
                  name: base, mode: 'lines', line: { color: colors[base] } })
  }
  const idx = { A: 0, C: 1, G: 2, T: 3 }
  for (const r of recs.filter(r => r.calls !== null)) {
    for (const b of r.calls.bases) {
      calls[idx[b]].x.push(r.position)
      calls[idx[b]].y.push(b)
      calls[idx[b]].text.push(`${b} (pos ${r.calls.pos})`)
    }
  }
  let xRange = [0, 500]
  if (data.chartConfig && data.chartConfig.x && data.chartConfig.x.axis && data.chartConfig.x.axis.range)
    xRange = data.chartConfig.x.axis.range
  Plotly.newPlot(traceChart, calls.concat(traces), {
    yaxis: { title: 'signal', domain: [0, 0.6] },
    yaxis2: { title: 'basecalls', domain: [0.7, 1], categoryorder: 'category descending' },
    xaxis: { title: 'Trace signal position', range: xRange, zeroline: false }
  }, { displayModeBar: true })
})()

// Decomposition
;(function () {
  const d = data.decomposition || { x: [], y: [] }
  Plotly.newPlot('decomposition-chart', [{ x: d.x, y: d.y, mode: 'lines+markers' }], {
    title: d.title || '',
    xaxis: { title: 'InDel length (bp)', zeroline: false },
    yaxis: { title: 'Decomposition error' }
  }, { displayModeBar: true })
})()

// Alignments
function chunked(seq, n) {
  const ret = []
  for (let i = 0; i < seq.length; i += n) ret.push(seq.slice(i, i + n))
  return ret
}
function chunkedAlignment(s1, s2, n) {
  return zip(chunked(s1, n), chunked(s2, n)).map(([l1, l2]) => {
    let m = ''
    for (const [c1, c2] of zip(l1, l2)) m += c1 === c2 ? '|' : ' '
    return [l1, m, l2]
  })
}
function alignmentText(alt, ref, n) {
  const labelWidth = Math.max(alt.label.length, ref.label.length)
  const numberWidth = Math.max(String(alt.sequence.length).length,
                               String(ref.startPosition + ref.sequence.length - 1).length)
  let pos1 = 1
  let pos2 = ref.isReverseComplement ? ref.startPosition + ref.sequence.length - 1 : ref.startPosition
  let body = ''
  for (const [s1, m, s2] of chunkedAlignment(alt.alignmentString, ref.alignmentString, n)) {
    body += `${alt.label.padStart(labelWidth)}  ${String(pos1).padStart(numberWidth)} ${s1}\n` +
            `${' '.repeat(labelWidth + numberWidth + 2)} ${m}\n` +
            `${ref.label.padStart(labelWidth)}  ${String(pos2).padStart(numberWidth)} ${s2}\n\n`
    pos1 += ungapped(s1).length
    pos2 += (ref.isReverseComplement ? -1 : 1) * ungapped(s2).length
  }
  const header = s => `>${s.chromosome}:${s.startPosition}-${s.startPosition + s.sequence.length - 1}` +
    `${s.isReverseComplement ? '_reverse' : '_forward'}` +
    `${s.alleleFraction ? ` (Estimated allelic fraction: ${s.alleleFraction})` : ''}`
  return `${header(alt)}\n${chunked(alt.sequence, n + 20).join('\n')}\n\n` +
         `${header(ref)}\n${chunked(ref.sequence, n + 20).join('\n')}\n\n${body}`
}
function showAlignment(id, alt, ref, score) {
  const el = document.getElementById(id)
  if (!alt.alignmentString || !ref.alignmentString) { el.innerHTML = '<p class="text-muted">n/a</p>'; return }
  el.innerHTML = `<pre>${score ? `Alignment score: ${score}\n\n` : ''}${escapeHtml(alignmentText(alt, ref, 80))}</pre>`
}

const alt1 = { sequence: ungapped(data.alt1align || ''), alignmentString: data.alt1align, isReverseComplement: false,
               chromosome: 'Alt1', startPosition: 1, label: 'Alt1', alleleFraction: data.allele1fraction }
const ref1 = { sequence: ungapped(data.ref1align || ''), alignmentString: data.ref1align, isReverseComplement: data.ref1forward === 0,
               chromosome: data.ref1chr, startPosition: data.ref1pos, label: 'Ref' }
const alt2 = { sequence: ungapped(data.alt2align || ''), alignmentString: data.alt2align, isReverseComplement: false,
               chromosome: 'Alt2', startPosition: 1, label: 'Alt2', alleleFraction: data.allele2fraction }
const ref2 = { sequence: ungapped(data.ref2align || ''), alignmentString: data.ref2align, isReverseComplement: data.ref2forward === 0,
               chromosome: data.ref2chr, startPosition: data.ref2pos, label: 'Ref' }
showAlignment('alignment-chart-1', alt1, ref1, data.align1score)
showAlignment('alignment-chart-2', alt2, ref2, data.align2score)
showAlignment('alignment-chart-3', { ...alt1, alignmentString: data.allele1align },
                                   { ...alt2, alignmentString: data.allele2align }, data.align3score)

// Variants table
;(function () {
  const v = data.variants
  const el = document.getElementById('variants-table')
  if (!v || !v.rows || v.rows.length === 0) { el.innerHTML = '<p class="text-muted">No variants called.</p>'; return }
  window.showVariant = i => {
    Plotly.relayout(traceChart, { 'xaxis.range': v.xranges[i] })
    traceChart.scrollIntoView()
  }
  el.innerHTML = `<table class="table table-sm table-striped table-hover">
    <thead><tr><th></th>${v.columns.map(c => `<th>${escapeHtml(c)}</th>`).join('')}</tr></thead>
    <tbody>${v.rows.map((row, i) => `<tr>
      <td class="show" title="Show in trace viewer" onclick="showVariant(${i})">&#128200;</td>
      ${row.map((val, j) => `<td title="${escapeHtml(v.columns[j])}">${escapeHtml(val)}</td>`).join('')}
    </tr>`).join('')}</tbody></table>`
})()
</script>
</body>
</html>
"""


COMP = str.maketrans('ACGTNRYSWKMBDHV-', 'TGCANYRSWMKVHDB-')
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
         'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}


def rc(seq):
    return seq.upper().translate(COMP)[::-1]


def matches(seq, pattern):
    return len(seq) == len(pattern) and all(b in IUPAC.get(p, p) for b, p in zip(seq, pattern))


def allele_alignment(d, k):
    """Alignment of allele k in reference-forward orientation: (alt, ref, positions)."""
    alt, ref = d.get(f'alt{k}align'), d.get(f'ref{k}align')
    if not alt or not ref:
        return None
    if d.get(f'ref{k}forward') == 0:
        alt, ref = rc(alt), rc(ref)
    pos, p = [], d.get(f'ref{k}pos', 1)
    for c in ref:
        if c == '-':
            pos.append(None)
        else:
            pos.append(p)
            p += 1
    return alt.upper(), ref.upper(), pos


def find_indels(alt, ref, pos):
    """Internal gap runs: deletions (gap in alt) and insertions (gap in ref)."""
    covered = [i for i, c in enumerate(alt) if c != '-']
    if not covered:
        return []
    first, last = covered[0], covered[-1]
    out, i = [], first
    while i <= last:
        if alt[i] == '-' and ref[i] != '-':
            j = i
            while j <= last and alt[j] == '-' and ref[j] != '-':
                j += 1
            out.append({'type': 'deletion', 'length': j - i, 'start': pos[i], 'end': pos[j - 1], 'seq': ref[i:j]})
            i = j
        elif ref[i] == '-':
            j = i
            while j <= last and ref[j] == '-':
                j += 1
            before = next((pos[x] for x in range(i - 1, -1, -1) if pos[x] is not None), None)
            out.append({'type': 'insertion', 'length': j - i, 'start': before, 'end': before, 'seq': alt[i:j]})
            i = j
        else:
            i += 1
    return out


def locate_guide(seq, start, grna, pam, cut_offset):
    """Find the guide (<=2 mismatches, both strands) in the reference slice."""
    g = grna.upper().replace('U', 'T').replace(' ', '')
    if len(g) == 20 + len(pam) and matches(g[-len(pam):], pam):
        g = g[:-len(pam)]  # 20 nt guide given together with its PAM
    best = None
    for strand, q in (('+', g), ('-', rc(g))):
        for i in range(len(seq) - len(q) + 1):
            mm = sum(a != b for a, b in zip(seq[i:i + len(q)], q))
            if mm <= 2 and (best is None or mm < best[0]):
                best = (mm, strand, i)
    if best is None:
        return None
    mm, strand, i = best
    L = len(g)
    if strand == '+':
        pam_pos = (start + i + L, start + i + L + len(pam) - 1)
        pam_seq = seq[i + L:i + L + len(pam)]
        cut_after = start + i + L - 1 + cut_offset   # last reference base before the cut
    else:
        pam_pos = (start + i - len(pam), start + i - 1)
        pam_seq = rc(seq[max(0, i - len(pam)):i])
        cut_after = start + i - 1 - cut_offset
    return {'strand': strand, 'mismatches': mm, 'guide': (start + i, start + i + L - 1), 'pam': pam_pos,
            'pam_seq': pam_seq, 'pam_ok': matches(pam_seq, pam), 'cut_after': cut_after, 'protospacer': g}


def distance_to_cut(ev, cut):
    if cut is None or ev['start'] is None:
        return None
    if ev['type'] == 'insertion':
        return abs(ev['start'] - cut)
    if ev['start'] - 1 <= cut <= ev['end']:
        return 0
    return min(abs(ev['start'] - 1 - cut), abs(ev['end'] - cut))


def crispr_panel(d, grna=None, pam='NGG', cut_offset=-3, flank=30):
    chrom = d.get('ref1chr') or d.get('ref2chr') or 'ref'
    alleles = []
    for k in (1, 2):
        a = allele_alignment(d, k)
        if a:
            alleles.append((k, a, find_indels(*a)))
    if not alleles:
        return None, ''
    guide = None
    if grna:
        _, (alt, ref, pos), _ = alleles[0]
        seq = ''.join(c for c in ref if c != '-')
        start = next(p for p in pos if p is not None)
        guide = locate_guide(seq, start, grna, pam, cut_offset)
    cut = guide['cut_after'] if guide else None
    centre = cut
    if centre is None:
        evs = [e for _, _, es in alleles for e in es if e['start'] is not None]
        centre = evs[0]['start'] if evs else None
    info = {'chrom': chrom, 'guide': guide, 'pam_pattern': pam, 'alleles': []}
    for k, _, evs in alleles:
        for e in evs:
            e['cut_distance'] = distance_to_cut(e, cut)
        info['alleles'].append({'allele': k, 'fraction': d.get(f'allele{k}fraction'), 'indels': evs})

    h = ['<h4 id="sec-crispr">CRISPR edit</h4>']
    h.append('<div class="help">Reference around the edit with the guide (<span class="g">blue</span>), '
             'PAM (<span class="pam">orange</span>) and the Cas9 cut site (<span class="cut"></span> red line, '
             f'{-cut_offset} bp upstream of the PAM). Below, each separated allele: '
             '<span class="del">deleted</span> bases, <span class="ins">inserted</span> bases, '
             '<span class="mm">mismatches</span>, <span class="nc">&middot;</span> = not covered (trimmed). '
             'Inside a repeat the exact position of an InDel is ambiguous (the alignment may place it a few bases '
             'left or right of the VCF position); its length is not. '
             + ('' if grna else '<b>Give the guide RNA (<code>--grna</code>) to see the guide, PAM and cut site.</b>')
             + '</div>')
    if grna and not guide:
        h.append(f'<div class="alert alert-warning">Guide {html.escape(grna)} was not found in the aligned reference '
                 '(&le;2 mismatches, both strands).</div>')
    if guide:
        g = guide
        pam_note = '' if g['pam_ok'] else f' <b style="color:#e03131">(does not match {html.escape(pam)})</b>'
        h.append('<table class="table table-sm alleles" style="width:auto"><tbody>'
                 f'<tr><th>Guide</th><td><code>{html.escape(g["protospacer"])}</code> ({g["strand"]} strand, '
                 f'{g["mismatches"]} mismatch(es)) &nbsp; {html.escape(chrom)}:{g["guide"][0]}-{g["guide"][1]}</td></tr>'
                 f'<tr><th>PAM</th><td><code>{html.escape(g["pam_seq"])}</code>{pam_note} &nbsp; '
                 f'{html.escape(chrom)}:{g["pam"][0]}-{g["pam"][1]}</td></tr>'
                 f'<tr><th>Cut site</th><td>between {html.escape(chrom)}:{g["cut_after"]} and {g["cut_after"] + 1}</td></tr>'
                 '</tbody></table>')
    rows = []
    for a in info['alleles']:
        frac = a['fraction']
        frac_s = f'{frac:.0%}' if isinstance(frac, (int, float)) else ''
        if not a['indels']:
            rows.append(f'<tr><td>Alt{a["allele"]}</td><td>{frac_s}</td><td>no InDel (same length as reference)</td><td></td><td></td></tr>')
        for e in a['indels']:
            sign = '&minus;' if e['type'] == 'deletion' else '+'
            loc = f'{chrom}:{e["start"]}' + (f'-{e["end"]}' if e['end'] != e['start'] else '')
            dist = '' if e['cut_distance'] is None else ('spans the cut site' if e['cut_distance'] == 0 else f'{e["cut_distance"]} bp')
            rows.append(f'<tr><td>Alt{a["allele"]}</td><td>{frac_s}</td>'
                        f'<td><b style="font-size:16px">{sign}{e["length"]} bp</b> {e["type"]} '
                        f'<code>{html.escape(e["seq"])}</code></td><td>{html.escape(loc)}</td><td>{dist}</td></tr>')
    h.append('<table class="table table-sm table-bordered alleles" style="width:auto"><thead><tr>'
             '<th>Allele</th><th>Fraction</th><th>Edit</th><th>Position</th><th>Distance to cut</th>'
             '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table>')
    v = d.get('variants') or {}
    vrows = [dict(zip(v.get('columns', []), r)) for r in v.get('rows', [])]
    vindels = [r for r in vrows if r.get('type') in ('Deletion', 'Insertion')]
    if vindels:
        items = []
        for r in vindels:
            n = len(str(r.get('ref', ''))) - len(str(r.get('alt', '')))
            items.append(f'<code>{html.escape(str(r.get("chr")))}:{r.get("pos")} {html.escape(str(r.get("ref")))}&rarr;'
                         f'{html.escape(str(r.get("alt")))}</code> ({"&minus;" if n > 0 else "+"}{abs(n)} bp, {html.escape(str(r.get("filter")))})')
        h.append('<p style="font-size:13px">Same InDel in the variant table (VCF, left-aligned): ' + '; '.join(items) + '</p>')
    info['vcf_indels'] = [{k: r.get(k) for k in ('chr', 'pos', 'ref', 'alt', 'filter', 'genotype')} for r in vindels]
    if centre is not None:
        lo, hi = centre - flank, centre + flank
        gset = set(range(guide['guide'][0], guide['guide'][1] + 1)) if guide else set()
        pset = set(range(guide['pam'][0], guide['pam'][1] + 1)) if guide else set()
        lines = []
        for k, (alt, ref, pos), _ in alleles:
            covered = [i for i, c in enumerate(alt) if c != '-']
            first, last = (covered[0], covered[-1]) if covered else (0, -1)
            cols, prev = [], None
            for i, p in enumerate(pos):
                q = p if p is not None else prev
                if q is not None and lo <= q <= hi:
                    cols.append(i)
                if p is not None:
                    prev = p
            if not cols:
                continue
            rref, ralt = [], []
            for i in cols:
                rb, ab = ref[i], alt[i]
                cls = 'g' if pos[i] in gset else 'pam' if pos[i] in pset else ''
                rref.append(f'<span class="b {cls}">{rb}</span>')
                if i < first or i > last:
                    acls, ab = 'nc', '&middot;'
                elif ab == '-' and rb != '-':
                    acls = 'del'
                elif rb == '-':
                    acls = 'ins'
                elif ab != rb:
                    acls = 'mm'
                else:
                    acls = ''
                ralt.append(f'<span class="b {acls}">{ab}</span>')
                if cut is not None and pos[i] == cut:
                    rref.append('<span class="cut"></span>')
                    ralt.append('<span class="cut"></span>')
            p0 = next((pos[i] for i in cols if pos[i] is not None), '')
            p1 = next((pos[i] for i in reversed(cols) if pos[i] is not None), '')
            frac = d.get(f'allele{k}fraction')
            frac_s = f'{frac:.0%}' if isinstance(frac, (int, float)) else ''
            edits = ', '.join(('&minus;' if e['type'] == 'deletion' else '+') + f'{e["length"]} bp'
                              for e in info['alleles'][[a['allele'] for a in info['alleles']].index(k)]['indels']
                              if e['start'] is not None and lo - e['length'] <= e['start'] <= hi) or 'no InDel'
            lines.append(f'<div><span class="lab">Ref</span><span class="num">{p0}</span>{"".join(rref)}'
                         f'<span style="color:#999;padding-left:6px">{p1}</span></div>'
                         f'<div><span class="lab">Alt{k}</span><span class="num"></span>{"".join(ralt)}'
                         f'<span style="padding-left:8px"><b>{edits}</b> <span style="color:#666">{frac_s}</span></span></div>'
                         '<div style="height:10px"></div>')
        h.append('<div class="edit">' + ''.join(lines) + '</div>')
    return info, '\n'.join(h)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('json', help='tracy decompose JSON (<prefix>.json)')
    ap.add_argument('out', nargs='?', help='output HTML [<prefix>.html]')
    ap.add_argument('-g', '--grna', help="guide RNA protospacer 5'->3' (20 nt, PAM optional)")
    ap.add_argument('--pam', default='NGG', help='PAM pattern, IUPAC [NGG = SpCas9]')
    ap.add_argument('--cut-offset', type=int, default=-3, help='cut site relative to the PAM [-3 = Cas9]')
    args = ap.parse_args()
    out_path = args.out or os.path.splitext(args.json)[0] + ".html"
    with open(args.json) as f:
        data = json.load(f)
    title = os.path.basename(str(data.get("meta", {}).get("arguments", {}).get("input", args.json)))
    info, panel = crispr_panel(data, args.grna, args.pam.upper(), args.cut_offset)
    if info is not None:
        with open(os.path.splitext(args.json)[0] + ".crispr.json", "w") as f:
            json.dump(info, f, indent=1)
    payload = json.dumps(data).replace("</", "<\\/")
    page = TEMPLATE.replace("__DATA__", payload).replace("__TITLE__", html.escape(title)).replace("__CRISPR__", panel)
    with open(out_path, "w") as f:
        f.write(page)
    print(out_path)


if __name__ == "__main__":
    main()
