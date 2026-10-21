#!/usr/bin/env python3
"""Plot audited LHT-04/NJ-lan RepeatMasker composition and divergence data."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/repeat_landscape_mpl')
import csv
import hashlib
import html
import io
import json
import re
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from audit_panel_alignment import require_matplotlib_panel_alignment

DELIVERY = Path(__file__).resolve().parents[1]
ROOT = DELIVERY / 'repeat_source'
PDF = DELIVERY / 'figures'
FIG = DELIVERY / 'figures'
DATA = DELIVERY / 'qa' / 'repeat_plot_data'
for p in (PDF, FIG, DATA):
    p.mkdir(parents=True, exist_ok=True)
SAMPLES = ('LHT-04', 'NJ-lan')
STEMS = {'LHT-04': 'lht04_final.shortids.fasta', 'NJ-lan': 'njlan_final.shortids.fasta'}
DIVSUMS = {'LHT-04': 'lht04.divsum', 'NJ-lan': 'njlan.divsum'}
CATEGORIES = ['Unclassified', 'SINEs', 'Rolling-circles', 'LTRs', 'LINEs',
              'DNA transposons', 'Small RNA', 'Other repeats']
DISPLAY_CATEGORIES = ['Unclassified', 'SINEs', 'Rolling-\ncircles', 'LTRs', 'LINEs',
                      'DNA\ntransposons', 'Small\nRNA', 'Other\nrepeats']
SHORT_CATEGORIES = ['Uncl.', 'SINE', 'RC', 'LTR', 'LINE', 'DNA', 'RNA', 'Other']
SPECIES_LABELS = {'LHT-04': 'LHT-04', 'NJ-lan': 'NJ-lan'}
COLORS = ['#aaaaaa', '#cc6b3c', '#367cad', '#ef5225', '#61418a',
          '#287b19', '#c4d8f3', '#e8d944']
TABLE_ROWS = ['Unclassified:', 'SINEs:', 'Rolling-circles', 'LTR elements:',
              'LINEs:', 'DNA transposons', 'Small RNA:']
plt.rcParams.update({'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'], 'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 12, 'axes.linewidth': .75,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})

def write_tsv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def parse_tbl(sample):
    path = ROOT / sample / (STEMS[sample] + '.tbl')
    text = path.read_text()
    span = int(re.search(r'total length:\s+(\d+)', text)[1])
    masked = int(re.search(r'bases masked:\s+(\d+)', text)[1])
    entries = {}
    for line in text.splitlines():
        m = re.match(r'\s*(.*?)\s+(\d+)\s+(\d+) bp\s+([\d.]+) %', line)
        if m:
            entries[m[1]] = {'elements': int(m[2]), 'bp': int(m[3]),
                             'reported_pct': float(m[4])}
    values = [entries[k]['bp'] for k in TABLE_ROWS]
    values.append(sum(entries[k]['bp'] for k in ('Satellites:', 'Simple repeats:', 'Low complexity:')))
    assert sum(values) == masked, (sample, sum(values), masked)
    return {'sample': sample, 'span': span, 'masked': masked,
            'composition': values, 'entries': entries,
            'pct': np.array(values) / span * 100}

official = (DELIVERY / 'scripts' / 'createRepeatLandscape.pl').read_text()
NAME_MAP = dict(re.findall(r'"([^"\n]+)"\s*=>\s*"([^"\n]+)"', official))
EXCLUDED = {'Simple_repeat', 'Low_complexity', 'Satellite', 'Structural_RNA', 'Segmental', 'ARTEFACT'}

def normalize_class(name):
    name = name.replace('?', '')
    return NAME_MAP.get(name, name)

def parse_divsum(sample):
    text = (ROOT / sample / DIVSUMS[sample]).read_text()
    families, coverage = text.split('Coverage for each repeat class and divergence (Kimura)')
    summary = []
    for line in families.splitlines():
        fields = line.split('\t')
        if len(fields) == 5 and fields[2].isdigit():
            summary.append({'class': fields[0], 'repeat': fields[1],
                            'abs_bp': int(fields[2]), 'well_characterized_bp': int(fields[3]),
                            'weighted_kimura_pct': None if fields[4] == '----' else float(fields[4])})
    lines = [line.split() for line in coverage.splitlines() if line.strip()]
    assert lines[0][0] == 'Div'
    raw_classes = lines[0][1:]
    bins, raw = [], []
    for row in lines[1:]:
        assert len(row) == len(raw_classes) + 1
        bins.append(int(row[0])); raw.append([int(x) for x in row[1:]])
    bins = np.array(bins); raw = np.array(raw, dtype=np.int64)
    assert np.array_equal(bins, np.arange(len(bins)))
    fixed = defaultdict(lambda: np.zeros(len(bins), dtype=np.int64))
    for i, cl in enumerate(raw_classes):
        fixed[normalize_class(cl)] += raw[:, i]
    return {'bins': bins, 'raw': raw, 'raw_classes': raw_classes,
            'fixed': dict(fixed), 'family_summary': summary}

samples = {s: parse_tbl(s) for s in SAMPLES}
for s in SAMPLES:
    samples[s].update(parse_divsum(s))

# Retain common superfamilies; pool rare classes within the same repeat type.
# Threshold uses the larger genome fraction across the two samples, not row counts.
allclasses = sorted(set().union(*(d['fixed'] for d in samples.values())))
mapping = {}
for cl in allclasses:
    if cl in EXCLUDED:
        mapping[cl] = None
        continue
    major = cl.split('/')[0]
    max_pct = max(samples[s]['fixed'].get(cl, np.array([0])).sum() / samples[s]['span'] * 100 for s in SAMPLES)
    if cl == 'Unknown':
        label = 'Unknown'
    elif cl == 'RC/Helitron':
        label = cl
    elif major in ('DNA', 'LINE', 'LTR', 'SINE'):
        label = cl if '/' in cl and max_pct >= .15 else major + '/Other'
    else:
        label = cl
    mapping[cl] = label

priority = {'DNA': 0, 'LINE': 1, 'LTR': 2, 'RC': 3, 'SINE': 4, 'PLE': 5,
            'Other': 6, 'Unknown': 7}
labels = sorted({v for v in mapping.values() if v is not None},
    key=lambda x: (priority.get(x.split('/')[0], 6), x.endswith('/Other'), x))
colors = {}
palettes = {'DNA': plt.cm.YlGn, 'LINE': plt.cm.Purples, 'LTR': plt.cm.OrRd}
for major, palette in palettes.items():
    group = [x for x in labels if x.split('/')[0] == major]
    for label, value in zip(group, np.linspace(.40, .90, len(group))):
        colors[label] = matplotlib.colors.to_hex(palette(value))
for i, label in enumerate(labels):
    colors.setdefault(label, {'Unknown': '#aaaaaa', 'RC/Helitron': '#367cad'}.get(label, matplotlib.colors.to_hex(plt.cm.tab20(i % 20))))

composition_rows, raw_rows, grouped_rows, mapping_rows, audit = [], [], [], [], {}
for cl in allclasses:
    mapping_rows.append({'normalized_class': cl, 'plotted_group': mapping[cl] or 'EXCLUDED_NON_TE',
                         'color': colors.get(mapping[cl], '')})
for s, d in samples.items():
    grouped = {label: np.zeros(len(d['bins']), dtype=np.int64) for label in labels}
    for cl, v in d['fixed'].items():
        if mapping[cl] is not None:
            # RepeatMasker overlap adjudication can emit a tiny negative cell
            # after subtraction (LHT-04 has one -6 bp cell).  Retain the raw
            # value in the source table and clip only the plotted aggregate.
            grouped[mapping[cl]] += np.maximum(v, 0)
    d['grouped'] = grouped
    assert all(np.all(v >= 0) for v in grouped.values()), 'A plotted aggregate has negative coverage'
    d['te_bp'] = sum(int(v.sum()) for v in grouped.values())
    d['tail_bp'] = sum(int(v[d['bins'] > 50].sum()) for v in grouped.values())
    for cat, bp in zip(CATEGORIES, d['composition']):
        composition_rows.append({'sample': s, 'category': cat, 'occupied_bp': bp,
                                  'assembly_bp': d['span'], 'genome_pct': bp / d['span'] * 100})
    for i, b in enumerate(d['bins']):
        for j, cl in enumerate(d['raw_classes']):
            raw_rows.append({'sample': s, 'kimura_bin_floor_pct': int(b), 'class': cl,
                'normalized_class': normalize_class(cl), 'well_characterized_bp': int(d['raw'][i, j]),
                'assembly_bp': d['span'], 'genome_pct': d['raw'][i, j] / d['span'] * 100})
        for label, v in grouped.items():
            grouped_rows.append({'sample': s, 'kimura_bin_floor_pct': int(b), 'group': label,
                'well_characterized_bp': int(v[i]), 'assembly_bp': d['span'],
                'genome_pct': v[i] / d['span'] * 100, 'shown_in_figure': int(b <= 50)})
    audit[s] = {'genome_bp': d['span'], 'masked_bp': d['masked'], 'masked_pct': d['masked']/d['span']*100,
                'composition_sum_bp': sum(d['composition']), 'composition_sum_equals_masked': True,
                'divsum_coverage_well_characterized_bp': int(d['raw'].sum()),
                'divsum_family_well_characterized_bp': sum(r['well_characterized_bp'] for r in d['family_summary']),
                'family_minus_coverage_well_characterized_bp': sum(r['well_characterized_bp'] for r in d['family_summary'])-int(d['raw'].sum()),
                'te_landscape_bp_all_emitted_bins': d['te_bp'],
                'te_landscape_bp_bins_0_to_50': d['te_bp']-d['tail_bp'], 'te_landscape_bp_bins_51_to_70': d['tail_bp'],
                'non_te_bp_excluded_from_emitted_coverage': int(d['raw'].sum())-d['te_bp'],
                'negative_raw_cells_retained_in_raw_table': [{'bin': int(d['bins'][i]), 'class': d['raw_classes'][j], 'bp': int(d['raw'][i,j])} for i,j in zip(*np.where(d['raw'] < 0))],
                'negative_bp_clipped_only_for_plot': int(-d['raw'][d['raw'] < 0].sum()),
                'plotting_adjustment': 'Negative divergence cells are retained in kimura_all_classes_all_bins.tsv and replaced by zero only in plotted aggregates.',
                'all_plotted_aggregates_nonnegative': True,
                'max_kimura_bin': int(d['bins'].max()), 'raw_classes': d['raw_classes'],
                'plot_groups': labels, 'CpG_adjustment': 'RepeatMasker default; -noCpGMod was not used',
                'coverage_scope': 'Official calcDivergenceFromAlign emits bins 0 through 70. Simple/low-complexity/artefact families are assigned divergence 100 and omitted from this table.'}
    write_tsv(DATA/(s+'_family_divergence_summary.tsv'), d['family_summary'])

write_tsv(DATA/'repeat_composition.tsv', composition_rows)
write_tsv(DATA/'kimura_all_classes_all_bins.tsv', raw_rows)
write_tsv(DATA/'kimura_plotted_groups_all_bins.tsv', grouped_rows)
write_tsv(DATA/'landscape_class_mapping.tsv', mapping_rows)
(DATA/'validation.json').write_text(json.dumps(audit, indent=2) + '\n')

composition_max = max(40, np.ceil(max(d['pct'].max() for d in samples.values()) / 5) * 5 + 2)
landscape_max = max(np.sum(list(d['grouped'].values()), axis=0)[d['bins'] <= 50].max()/d['span']*100 for d in samples.values())
landscape_step = .5 if landscape_max <= 5 else 1
landscape_max = np.ceil(landscape_max / landscape_step) * landscape_step
tooltip = {}

def axes_style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#e5e5e5', linewidth=.65)
    ax.tick_params(length=3, width=.6)
    for sp in ax.spines.values():
        sp.set_color('#262626')

def draw_row(fig, grid, sample, letters):
    d = samples[sample]
    ax = fig.add_subplot(grid[0]); bx = fig.add_subplot(grid[1])
    axes_style(ax); axes_style(bx)
    bars = ax.bar(np.arange(8), d['pct'], color=COLORS, edgecolor='#252525', linewidth=.65, width=.72)
    for i, (bar, val) in enumerate(zip(bars, d['pct'])):
        key = f'{sample}_composition_{i}'
        bar.set_gid(key)
        tooltip[key] = f'{sample} | {CATEGORIES[i]} | {d["composition"][i]:,} bp | {val:.4f}% of genome'
        label_y = val + .50 if val >= 6 else 6.6 + 1.4 * (i % 2)
        ax.text(bar.get_x()+bar.get_width()/2, label_y, f'{val:.2f}', ha='center', va='bottom', fontsize=7)
    ax.set_ylim(0, composition_max); ax.set_xlim(-.65, 7.65)
    ax.set_xticks(range(8), SHORT_CATEGORIES, ha='center')
    ax.tick_params(axis='x', labelsize=7)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_ylabel('Proportion of the genome (%)')
    ax.set_title(f'{SPECIES_LABELS[sample]} · Repeat composition', loc='left', pad=13, fontweight='bold')
    ax.text(-.17, 1.055, letters[0], transform=ax.transAxes, fontsize=14, fontweight='bold')
    idx = d['bins'] <= 50
    b = d['bins'][idx]
    bottom = np.zeros(len(b))
    for label in labels:
        values = d['grouped'][label][idx] / d['span'] * 100
        stacks = bx.bar(b, values, bottom=bottom, color=colors[label], width=.92, linewidth=0)
        for k, (patch, v) in enumerate(zip(stacks, values)):
            key = f'{sample}_kimura_{labels.index(label)}_{k}'
            patch.set_gid(key)
            tooltip[key] = f'{sample} | Kimura [{int(b[k])}, {int(b[k])+1})% | {label} | {v:.5f}% of genome'
        bottom += values
    bx.set_xlim(-.8, 50.8); bx.set_ylim(0, landscape_max*1.06)
    bx.set_xticks(np.arange(0, 51, 10)); bx.xaxis.set_minor_locator(MultipleLocator(5))
    bx.yaxis.set_major_locator(MultipleLocator(landscape_step))
    bx.grid(axis='x', color='#eeeeee', linewidth=.5)
    bx.set_xlabel('Kimura divergence (%)')
    bx.set_ylabel('Proportion of the genome (%)')
    bx.set_title('Repeat landscape', loc='left', pad=13, fontweight='bold')
    bx.text(-.14, 1.055, letters[1], transform=bx.transAxes, fontsize=14, fontweight='bold')
    if len(grid) > 2:
        lx = fig.add_subplot(grid[2]); lx.axis('off')
        handles = [Patch(facecolor=colors[label], label=label) for label in reversed(labels)]
        lx.legend(handles=handles, loc='upper left', frameon=False, borderaxespad=0,
                  fontsize=7.2, handlelength=1.3, handletextpad=.55, labelspacing=.3)
    return ax, bx

def save(fig, name):
    fig.savefig(PDF/(name+'.pdf'))
    fig.savefig(FIG/(name+'.png'), dpi=600)
    fig.savefig(FIG/(name+'.tiff'), dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    out = io.StringIO(); fig.savefig(out, format='svg')
    svg = out.getvalue()
    for key, title in tooltip.items():
        svg = svg.replace(f'<g id="{key}">', f'<g id="{key}"><title>{html.escape(title)}</title>')
    (FIG/(name+'.svg')).write_text(svg)
    plt.close(fig)

fig = plt.figure(figsize=(7.2, 9.2))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.48], left=.11, right=.97, bottom=.23, top=.86, hspace=.72, wspace=.38)
combined_axes = []
for row, s in enumerate(SAMPLES):
    combined_axes.extend(draw_row(fig, [gs[row,i] for i in range(2)], s, ('a','b') if row == 0 else ('c','d')))
fig.text(.07, .975, 'Repeat composition and divergence landscapes', fontsize=13, fontweight='bold', va='top')
fig.text(.07, .942, 'LHT-04: 38.97% repeats, 115.73 Mb     |     NJ-lan: 46.47% repeats, 291.89 Mb', fontsize=8.5, va='top', color='#444444')
handles = [Patch(facecolor=colors[label], label=label) for label in reversed(labels)]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .073), ncol=5,
           frameon=False, fontsize=6, handlelength=1.2, handletextpad=.35,
           columnspacing=.75, labelspacing=.22)
fig.text(.07, .022, 'Other repeats = satellites + simple repeats + low complexity. Landscape: CpG-adjusted Kimura, 1% bins; non-TE repeats excluded.\nComposition uses non-overlapping RepeatMasker .tbl totals; landscape uses official divergence-alignment counts.', fontsize=6, color='#444444')
require_matplotlib_panel_alignment(
    fig,
    axes=combined_axes,
    panel_ids=['a', 'b', 'c', 'd'],
    row_groups=(('a', 'b'), ('c', 'd')),
    column_groups=(('a', 'c'), ('b', 'd')),
    json_out=FIG/'Figure2_repeat_composition_landscapes.alignment.json',
    overlay_svg=FIG/'Figure2_repeat_composition_landscapes.alignment.svg',
    require_panel_labels=True,
    strict=True,
)
save(fig, 'Figure2_repeat_composition_landscapes')
print(json.dumps(audit, indent=2))
