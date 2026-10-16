#!/usr/bin/env python3
"""Summarize many indigo/tracy results into one table and one HTML index.

Usage: indigo_summary.py <results folder>

Finds every <name>.json below the folder and writes
  <folder>/summary.tsv   one row per sample
  <folder>/index.html    sortable/filterable table linking each sample's report

Columns:
  best_indel     InDel length with the lowest decomposition error (0 = no InDel,
                 negative = deletion, positive = insertion)
  best_error     its decomposition error
  second_indel   the next best InDel length, second_error its error
  hetindel       1 if tracy detected a breakpoint where the trace becomes mixed
  allele1/2_frac estimated fractions of the two alleles
  call           heuristic: WT (only 0 fits), het (two lengths fit, one is 0),
                 biallelic (two non-zero lengths fit), hom (one non-zero length),
                 check (no clear minimum)
  alt1/2_edit    InDels of each allele from its alignment, e.g. -7bp@1225(cut+0)
                 (cut distance only when a guide RNA was given)
  cut_site       Cas9 cut site (last base before the cut), if a guide RNA was given
  indels         InDel variants from the variant table (chr:pos ref>alt)
"""
import csv
import html
import json
import os
import sys

FIELDS = ['sample', 'status', 'call', 'best_indel', 'best_error', 'second_indel', 'second_error',
          'hetindel', 'allele1_frac', 'allele2_frac', 'alt1_edit', 'alt2_edit', 'cut_site',
          'align1_score', 'n_variants', 'n_pass', 'indels', 'ref_chr', 'ref_pos', 'report']


def classify(points):
    """points: list of (indel, error) sorted by error."""
    if len(points) < 2:
        return 'check'
    (i1, e1), (i2, e2) = points[0], points[1]
    rest = [e for _, e in points[2:]] or [e2]
    background = sorted(rest)[len(rest) // 2]
    if background <= 0 or e1 > 0.5 * background:
        return 'check'
    two = e2 < 0.5 * background
    if not two:
        return 'WT' if i1 == 0 else 'hom'
    if 0 in (i1, i2):
        return 'het'
    return 'biallelic'


def summarize(path, folder):
    with open(path) as f:
        d = json.load(f)
    if not isinstance(d, dict) or 'decomposition' not in d:
        return None  # some other json file
    name = os.path.splitext(os.path.basename(path))[0]
    row = {'sample': name, 'status': 'ok', 'report': os.path.relpath(os.path.splitext(path)[0] + '.html', folder)}
    dec = d.get('decomposition') or {}
    pts = sorted(zip(dec.get('x', []), dec.get('y', [])), key=lambda p: (p[1], abs(p[0])))
    if pts:
        row['best_indel'], row['best_error'] = pts[0]
    if len(pts) > 1:
        row['second_indel'], row['second_error'] = pts[1]
    row['call'] = classify(pts)
    row['hetindel'] = d.get('hetindel', '')
    row['allele1_frac'] = d.get('allele1fraction', '')
    row['allele2_frac'] = d.get('allele2fraction', '')
    row['align1_score'] = d.get('align1score', '')
    row['ref_chr'] = d.get('ref1chr', '')
    row['ref_pos'] = d.get('ref1pos', '')
    v = d.get('variants') or {}
    cols = v.get('columns', [])
    rows = [dict(zip(cols, r)) for r in v.get('rows', [])]
    row['n_variants'] = len(rows)
    row['n_pass'] = sum(1 for r in rows if r.get('filter') == 'PASS')
    cj = os.path.splitext(path)[0] + '.crispr.json'
    if os.path.exists(cj):
        with open(cj) as f:
            c = json.load(f)
        for a in c.get('alleles', []):
            evs = []
            for e in a.get('indels', []):
                tag = ('-' if e['type'] == 'deletion' else '+') + f"{e['length']}bp@{e['start']}"
                if e.get('cut_distance') is not None:
                    tag += f"(cut+{e['cut_distance']})"
                evs.append(tag)
            row[f"alt{a['allele']}_edit"] = ';'.join(evs) or 'none'
        if c.get('guide'):
            row['cut_site'] = c['guide']['cut_after']
    row['indels'] = ';'.join(f"{r.get('chr')}:{r.get('pos')}{r.get('ref')}>{r.get('alt')}"
                             for r in rows if r.get('type') in ('Deletion', 'Insertion'))
    return row


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    folder = os.path.abspath(sys.argv[1])
    # status from indigo_batch.py (if present): a failure there wins over a leftover json
    status = {}
    st = os.path.join(folder, 'status.tsv')
    if os.path.exists(st):
        with open(st) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                status[r['sample']] = r['status']
    out = []
    seen = set()
    for root, _, files in os.walk(folder):
        for fn in sorted(files):
            if not fn.endswith('.json') or fn.endswith('.crispr.json'):
                continue
            name = fn[:-5]
            try:
                r = summarize(os.path.join(root, fn), folder)
            except Exception as e:
                if os.path.basename(root) != name:
                    continue  # not an indigo output, ignore
                r = {'sample': name, 'status': f'unreadable json: {e}'}
            if r is None:
                continue
            if status.get(r['sample'], 'ok') not in ('ok', 'skipped'):
                r = {'sample': r['sample'], 'status': status[r['sample']]}
            out.append(r)
            seen.add(r['sample'])
    for name, st_ in status.items():
        if name not in seen:
            out.append({'sample': name, 'status': st_})
    out.sort(key=lambda r: r['sample'])

    with open(os.path.join(folder, 'summary.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(out)

    counts = {}
    for r in out:
        k = r.get('call') or 'failed'
        counts[k] = counts.get(k, 0) + 1
    head = ''.join(f'<th onclick="sortBy({i})">{html.escape(c)}</th>' for i, c in enumerate(FIELDS[:-1]))
    body = []
    for r in out:
        cells = []
        for c in FIELDS[:-1]:
            val = html.escape(str(r.get(c, '')))
            if c == 'sample' and r.get('report') and os.path.exists(os.path.join(folder, r['report'])):
                val = f'<a href="{html.escape(r["report"])}" target="_blank">{val}</a>'
            cells.append(f'<td>{val}</td>')
        cls = html.escape(r.get('call') or 'failed')
        body.append(f'<tr class="{cls}">' + ''.join(cells) + '</tr>')
    summary = ' · '.join(f'{k}: {v}' for k, v in sorted(counts.items()))
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>Indigo batch summary</title>
<style>
body{{font-family:system-ui,sans-serif;margin:20px}} table{{border-collapse:collapse;font-size:13px}}
th,td{{border:1px solid #ddd;padding:4px 8px;white-space:nowrap}} th{{background:#f3f3f3;cursor:pointer;position:sticky;top:0}}
tr.het td{{background:#fff8e1}} tr.hom td,tr.biallelic td{{background:#fde2e2}} tr.WT td{{background:#e8f5e9}}
tr.check td,tr.failed td{{background:#eeeeee;color:#555}} input{{padding:6px;width:300px;margin:8px 0}}
</style></head><body>
<h2>Indigo batch summary</h2><p>{len(out)} samples — {html.escape(summary)}</p>
<input id="q" placeholder="filter (e.g. het, -7, sample name)" oninput="filt()">
<p style="font-size:12px;color:#666">Click a header to sort, a sample name to open its report.
call is a heuristic from the decomposition curve — always check the trace.</p>
<table id="t"><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>
<script>
function filt(){{const q=document.getElementById('q').value.toLowerCase();
 for(const tr of document.querySelectorAll('#t tbody tr')) tr.style.display=tr.textContent.toLowerCase().includes(q)?'':'none'}}
let dir=1;function sortBy(i){{const tb=document.querySelector('#t tbody');dir=-dir;
 const rows=[...tb.rows].sort((a,b)=>{{const x=a.cells[i].textContent,y=b.cells[i].textContent;
 const nx=parseFloat(x),ny=parseFloat(y);return (isNaN(nx)||isNaN(ny)?x.localeCompare(y):nx-ny)*dir}});
 rows.forEach(r=>tb.appendChild(r))}}
</script></body></html>"""
    with open(os.path.join(folder, 'index.html'), 'w') as f:
        f.write(page)
    print(f'{len(out)} samples -> {folder}/summary.tsv, {folder}/index.html ({summary})', file=sys.stderr)


if __name__ == '__main__':
    main()
