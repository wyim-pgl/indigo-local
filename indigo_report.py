#!/usr/bin/env python3
"""Render a tracy/indigo JSON result as a standalone HTML report (like the Indigo web app).

Usage: indigo_report.py <prefix>.json [output.html]
"""
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
</style>
</head>
<body>
<div class="container-fluid">
  <h2>Indigo results</h2>
  <p class="text-muted" id="meta"></p>

  <h4 id="sec-variants">Variants</h4>
  <div id="variants-table"></div>

  <h4 id="sec-trace">Trace</h4>
  <div id="trace-chart" style="height:500px"></div>

  <h4 id="sec-decomp">Decomposition</h4>
  <div id="decomposition-chart" style="height:400px"></div>

  <h4 id="sec-align1">Alignment of Alt1 to reference</h4>
  <div id="alignment-chart-1"></div>

  <h4 id="sec-align2">Alignment of Alt2 to reference</h4>
  <div id="alignment-chart-2"></div>

  <h4 id="sec-align3">Alignment of Alt1 to Alt2</h4>
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


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    json_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(json_path)[0] + ".html"
    with open(json_path) as f:
        data = json.load(f)
    title = os.path.basename(str(data.get("meta", {}).get("arguments", {}).get("input", json_path)))
    payload = json.dumps(data).replace("</", "<\\/")
    html = TEMPLATE.replace("__DATA__", payload).replace("__TITLE__", title)
    with open(out_path, "w") as f:
        f.write(html)
    print(out_path)


if __name__ == "__main__":
    main()
