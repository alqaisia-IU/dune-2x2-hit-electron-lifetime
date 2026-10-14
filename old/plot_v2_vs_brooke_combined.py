#!/usr/bin/env python3
"""Overlay DUNE 2x2 V2 lifetime fits with Brooke Russell's slide reference.

Plotting-only. No packet processing or fitting. Brooke reference CSV consists
of approximately digitized May 16, 2025 slide-4 July 2024 lifetime values,
assigned to *assumed* acquisition timestamps. These are NOT verified matches,
and no Brooke measurement uncertainties were available in that source.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SERIES = [
    ('n2', 'Exact n=2, separate IO amplitudes', 'o', '#2166ac'),
    ('shared', 'Exact n=2+n=3, independent IO/n amplitudes', '^', '#17864a'),
    ('n3', 'Exact n=3, separate IO amplitudes', 's', '#d66d10'),
]


def inputs():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--v2', type=Path, default=Path('hit_based_v2_full/v2_lifetime_timeline.csv'),
                   help='Output of run_hit_based_v2_validation.py')
    p.add_argument('--brooke', type=Path,
                   default=Path(__file__).with_name('brooke_digitized_approx_v2.csv'),
                   help='Approximately digitized Brooke 2x2 slide lifetime points')
    p.add_argument('--outdir', type=Path, default=Path('hit_based_v2_full/brooke_overlay'))
    p.add_argument('--focus-max-ms', type=float, default=2.2,
                   help='Upper y-axis bound in main comparison panel')
    return p.parse_args()


def parse_v2(path):
    if not path.is_file():
        raise FileNotFoundError(f'V2 CSV not found: {path}. Check --v2 path.')
    d = pd.read_csv(path)
    need = {'time_CDT', 'packet_file'}
    for tag, _, _, _ in SERIES:
        need.add(tag + '_tau_ms')
        need.add(tag + '_tau_error_ms')
    absent = need - set(d.columns)
    if absent:
        raise ValueError(f'V2 CSV is missing columns: {sorted(absent)}')
    d['time_CDT'] = pd.to_datetime(d['time_CDT'], errors='coerce')
    d = d[d['time_CDT'].notna()].copy()
    for tag, _, _, _ in SERIES:
        for suffix in ('_tau_ms', '_tau_error_ms'):
            d[tag + suffix] = pd.to_numeric(d[tag + suffix], errors='coerce')
    if d.empty:
        raise ValueError('No readable acquisition times in V2 CSV')
    # One row per actual processed packet-file acquisition. No averaging of repeated rows.
    d = d.sort_values('time_CDT').drop_duplicates('packet_file', keep='last')
    return d


def parse_brooke(path):
    if not path.is_file():
        raise FileNotFoundError(f'Brooke reference CSV not found: {path}. Copy it beside this script.')
    b = pd.read_csv(path)
    need = {'acquisition_start_cdt', 'brooke_lifetime_ms'}
    if not need.issubset(b.columns):
        raise ValueError(f'Brooke CSV is missing columns: {sorted(need - set(b.columns))}')
    b['assumed_time_CDT'] = pd.to_datetime(b['acquisition_start_cdt'], errors='coerce')
    b['brooke_lifetime_ms'] = pd.to_numeric(b['brooke_lifetime_ms'], errors='coerce')
    b = b[b['assumed_time_CDT'].notna() & b['brooke_lifetime_ms'].gt(0)].copy()
    if b.empty:
        raise ValueError('Brooke CSV has no usable slide points')
    b = b.sort_values('assumed_time_CDT')
    if b['assumed_time_CDT'].duplicated().any():
        raise ValueError('Duplicate assumed Brooke timestamps; check input CSV.')
    return b


def add_series(ax, df, tag, name, marker, color, *, alpha=1.0, line=True, size=5.5):
    y = df[tag + '_tau_ms'].to_numpy(float)
    e = df[tag + '_tau_error_ms'].to_numpy(float)
    good = np.isfinite(y) & (y > 0)
    if not np.any(good):
        return 0
    x = df.loc[good, 'time_CDT']
    y = y[good]
    e = e[good]
    # Do not invent uncertainties: draw missing-error entries as markers only.
    err = np.where(np.isfinite(e) & (e > 0), e, 0.0)
    ax.errorbar(x, y, yerr=err, fmt=marker + ('-' if line else ''), color=color,
                elinewidth=1.25, linewidth=1.35, markersize=size, capsize=2.5,
                alpha=alpha, label=name, zorder=3)
    return len(x)


def add_brooke(ax, b):
    ax.plot(b['assumed_time_CDT'], b['brooke_lifetime_ms'],
            linestyle='--', marker='D', color='#333333', linewidth=1.6,
            markersize=6, markerfacecolor='white', markeredgewidth=1.7,
            label='Brooke 2024 combined detector (approx. slide digitization)', zorder=5)


def main():
    a = inputs()
    d = parse_v2(a.v2)
    b = parse_brooke(a.brooke)
    a.outdir.mkdir(parents=True, exist_ok=True)
    print(f'V2 acquisition rows: {len(d)}')
    for tag, *_ in SERIES:
        print(f'  valid {tag}: {int((d[tag+"_tau_ms"].notna() & (d[tag+"_tau_ms"] > 0)).sum())}')
    print(f'Brooke approximate slide points: {len(b)}')

    # Matched pairs for *visual context only*: CSV timestamps are assumed, not verified.
    matched = b.merge(d, left_on='assumed_time_CDT', right_on='time_CDT', how='left', validate='one_to_one')
    for tag, *_ in SERIES:
        matched[f'{tag}_minus_brooke_ms'] = matched[f'{tag}_tau_ms'] - matched['brooke_lifetime_ms']
    matchcols = ['acquisition_start_cdt', 'brooke_lifetime_ms', 'reference_status', 'packet_file']
    matchcols += [f'{tag}_{field}' for tag, *_ in SERIES for field in ('tau_ms', 'tau_error_ms')]
    matchcols += [f'{tag}_minus_brooke_ms' for tag, *_ in SERIES]
    matchcols = [c for c in matchcols if c in matched.columns]
    matchpath = a.outdir/'v2_vs_brooke_assumed_pairs.csv'
    matched[matchcols].to_csv(matchpath, index=False)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14.4, 9.0), sharex=True,
                                   gridspec_kw={'height_ratios': [1.45, 1.0]})
    add_series(ax1, d, *SERIES[0])
    add_series(ax1, d, *SERIES[1])
    add_brooke(ax1, b)
    ax1.set_ylim(0.35, a.focus_max_ms)
    ax1.set_ylabel('Effective attenuation time (ms)')
    ax1.set_title('DUNE 2x2 July 2024 | hit-based V2 vs Brooke Russell combined-detector result\n'
                  'Main comparison: n=2 and shared n=2+n=3', fontsize=14)
    ax1.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax1.grid(alpha=0.25)

    add_series(ax2, d, *SERIES[2], alpha=0.85)
    add_series(ax2, d, *SERIES[0], alpha=0.45, size=3.5)
    add_brooke(ax2, b)
    ax2.set_ylabel('Effective attenuation time (ms)')
    ax2.set_title('Diagnostic: n=3 alone (larger uncertainties and greater variability)', fontsize=11)
    ax2.grid(alpha=0.25)
    ax2.legend(loc='upper left', fontsize=8, framealpha=0.9)

    alltime = pd.concat([d['time_CDT'], b['assumed_time_CDT']])
    pad = pd.Timedelta(hours=3)
    ax1.set_xlim(alltime.min() - pad, alltime.max() + pad)
    ax2.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=6, maxticks=12))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
    ax2.set_xlabel('Acquisition start (CDT); Brooke horizontal positions are assumed')

    fig.text(0.07, 0.015,
             'Brooke reference: May 16, 2025 presentation, slide 4, "2x2 Demonstrator July 2024". '
             'Values digitized approximately; dates are assumed pairings.\n'
             'No Brooke error bars digitized. Methods differ in trigger grouping and selection; '
             'this is a qualitative time-trend comparison, NOT a validated run-matched agreement test.',
             fontsize=8.8, ha='left', va='bottom', color='#444444')
    fig.subplots_adjust(top=0.91, bottom=0.13, hspace=0.29, left=0.075, right=0.985)
    outplot = a.outdir/'v2_vs_brooke_combined.png'
    fig.savefig(outplot, dpi=190)
    plt.close(fig)

    # Separate focused plot (easier to read on slides)
    fig, ax = plt.subplots(figsize=(13.8, 5.8))
    add_series(ax, d, *SERIES[0]); add_series(ax, d, *SERIES[1]); add_brooke(ax, b)
    ax.set_ylim(0.35, a.focus_max_ms)
    ax.set_xlim(alltime.min() - pad, alltime.max() + pad)
    ax.set_title('DUNE 2x2: V2 (exact n=2 vs shared n=2+n=3) compared with Brooke (approx.)')
    ax.set_xlabel('Acquisition start (CDT); Brooke timestamps assumed')
    ax.set_ylabel('Electron lifetime / effective attenuation time (ms)')
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=6, maxticks=12))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
    ax.grid(alpha=0.25); ax.legend(loc='upper left', fontsize=8.5)
    fig.text(0.06, 0.01, 'Brooke slide 4: approximate digitization and ASSUMED timestamp pairings; no verified run matching or Brooke errors.',
             fontsize=9, color='#444444')
    fig.subplots_adjust(bottom=0.21, left=0.08, right=0.98, top=0.89)
    focusplot = a.outdir/'v2_vs_brooke_focus.png'
    fig.savefig(focusplot, dpi=190)
    plt.close(fig)

    nmatched = int(matched['packet_file'].notna().sum())
    print(f'Brooke assumed timestamps matching a V2 acquisition: {nmatched}/{len(b)}')
    print('NOTE: matches are based on previously assumed reference timestamps, not confirmed Brooke run IDs.')
    print('No Brooke uncertainties are available; do not interpret matched differences as significance.')
    print('Saved:', outplot)
    print('Saved:', focusplot)
    print('Saved:', matchpath)


if __name__ == '__main__':
    main()
