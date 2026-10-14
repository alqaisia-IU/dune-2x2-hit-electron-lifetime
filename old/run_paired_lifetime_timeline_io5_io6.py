#!/usr/bin/env python3
"""Compare frozen IO5 lifetimes with compact-group methods on IO5 AND IO6.

Requires the existing scripts in hit_based_LAST:
  run_paired_lifetime_timeline_all.py   (reuses the same grouping and fitting)
  packet_lifetime_FROZEN.py             (same master-clock unrolling as frozen)
  packet_lifetime_final.py
  packet_lifetime_v3_brooke.py

Frozen output files are READ ONLY.  New IO6 packet selections are built from the
original raw HDF5 (not from the IO5-selected frozen NPZ) and cached in outdir.
"""
import argparse
import json
import re
import zlib
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import run_paired_lifetime_timeline_all as core
from packet_lifetime_FROZEN import unroll_mclk, TICK_US


def args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--frozen-root', type=Path, action='append', default=None,
                    help='May be repeated; defaults to both frozen directories')
    ap.add_argument('--flow-root', type=Path, default=Path('/global/cfs/cdirs/dune/www/data/2x2/reflows/v9/flow/beam'))
    ap.add_argument('--outdir', type=Path, default=Path('paired_lifetime_IO5_IO6'))
    ap.add_argument('--only', type=str, default=None, help='Substring of filename for first-file validation')
    ap.add_argument('--bootstrap', type=int, default=30)
    ap.add_argument('--ext-io', type=int, choices=[5, 6], action='append', default=None,
                    help='New-method EXT choice (repeat to use both); defaults to IO5 and IO6')
    ap.add_argument('--ext-trigger-type', type=int, default=2)
    ap.add_argument('--chunk', type=int, default=2_000_000)
    ap.add_argument('--list-only', action='store_true')
    ap.add_argument('--force-io6-extract', action='store_true')
    ap.add_argument('--validate-io6-npz', type=Path, default=None,
                    help='Compare newly extracted IO6 selection to known original IO6 selection')
    return ap.parse_args()


def extract_npz_from_raw(packet_path, out_npz, ext_io, trig_type, chunk):
    """Reproduce frozen selection logic, changing only EXT IO to ext_io."""
    packet_path = Path(packet_path)
    if not packet_path.is_file():
        raise FileNotFoundError(f'Raw packet file missing: {packet_path}')
    with h5py.File(packet_path, 'r') as h:
        p = h['packets']
        tfields = ['packet_type', 'io_group', 'trigger_type', 'timestamp', 'receipt_timestamp']
        print(f'    reading timing for {len(p):,} raw packets ...', flush=True)
        tim = p.fields(tfields)[:]
        mclk = unroll_mclk(tim)
        ext_mask = ((tim['packet_type'] == 7) &
                    (tim['io_group'] == ext_io) &
                    (tim['trigger_type'] == trig_type))
        ext_t = np.sort(mclk[ext_mask])
        n_ext = len(ext_t)
        if n_ext == 0:
            raise RuntimeError(f'No EXT markers for IO{ext_io}, type {trig_type}')
        data_mask = tim['packet_type'] == 0
        data_mclk = mclk[data_mask]
        del mclk, tim

        ext_index = np.searchsorted(ext_t, data_mclk, side='right') - 1
        good = ext_index >= 0
        dt = np.full(len(data_mclk), -1, dtype=np.int64)
        dt[good] = data_mclk[good] - ext_t[ext_index[good]]
        dt_us = dt * TICK_US
        keep = good & (dt_us >= 0) & (dt_us <= 200.0)
        print(f'    IO{ext_io}: {n_ext:,} EXT markers, {len(data_mclk):,} charge packets, '
              f'{int(keep.sum()):,} selected', flush=True)

        parts = []
        cursor = 0
        sf = ['packet_type', 'io_group', 'io_channel', 'chip_id', 'channel_id', 'dataword']
        q_fields = ['io_group', 'io_channel', 'chip_id', 'channel_id', 'dataword']
        for start in range(0, len(p), chunk):
            z = p.fields(sf)[start:min(start + chunk, len(p))]
            data_only = z['packet_type'] == 0
            n_data = int(data_only.sum())
            use = keep[cursor:cursor+n_data]
            if use.any():
                parts.append(z[data_only][use][q_fields])
            cursor += n_data
        assert cursor == len(data_mclk), (cursor, len(data_mclk))
        if not parts:
            raise RuntimeError('No selected charge packets')
        selected = np.concatenate(parts)
        sel_dt = dt_us[keep]
        sel_ext = ext_index[keep]
        assert len(selected) == len(sel_dt) == len(sel_ext)
        out_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_npz, selected=selected, selected_dt_us=sel_dt,
                            selected_ext_index=sel_ext)
    return n_ext, len(selected)


def check_original_io6(candidate, known_npz):
    if known_npz is None:
        return
    if not known_npz.is_file():
        raise FileNotFoundError(f'Validation IO6 NPZ not found: {known_npz}')
    with np.load(candidate, allow_pickle=False) as got, np.load(known_npz, allow_pickle=False) as exp:
        for name in ['selected', 'selected_ext_index', 'selected_dt_us']:
            a = got[name]; b = exp[name]
            equal = (np.allclose(a, b, atol=1e-9, rtol=0) if name == 'selected_dt_us'
                     else np.array_equal(a, b))
            if not equal:
                raise RuntimeError(f'IO6 validation FAILED for {name}: '
                                   f'new shape {a.shape}, reference shape {b.shape}')
            print(f'    IO6 validation: {name} matches original ({len(a):,} rows)', flush=True)


def save_outputs(rows, mpv_rows, outdir):
    if rows:
        pd.DataFrame(rows).sort_values('time_CDT').to_csv(outdir / 'paired_lifetime_timeseries.csv', index=False)
    if mpv_rows:
        pd.DataFrame(mpv_rows).to_csv(outdir / 'new_method_mpvs.csv', index=False)


def plot_results(df, dest):
    df = df.sort_values('time_CDT')
    fig, ax = plt.subplots(figsize=(13, 6.5))
    styles = [
        ('frozen_unweighted_ms', 'Frozen IO5: unweighted 28-tick', 's--'),
        ('frozen_weighted_ms', 'Frozen IO5: weighted 28-tick', '^--'),
        ('new_io5_n2_tau_ms', 'New IO5: compact n=2', 'o-'),
        ('new_io5_nge2_tau_ms', 'New IO5: compact n>=2', 'd-'),
        ('new_io6_n2_tau_ms', 'New IO6: compact n=2', 'o-'),
        ('new_io6_nge2_tau_ms', 'New IO6: compact n>=2', 'd-'),
    ]
    for k, label, fmt in styles:
        if k not in df:
            continue
        y = pd.to_numeric(df[k], errors='coerce')
        mask = y.notna() & np.isfinite(y)
        if not mask.any():
            continue
        if k.startswith('new_'):
            e = pd.to_numeric(df[k.replace('_tau_ms','_tau_err_ms')], errors='coerce')
            ax.errorbar(df.loc[mask, 'time_CDT'], y.loc[mask], yerr=e.loc[mask],
                        fmt=fmt, capsize=2, label=label)
        else:
            ax.plot(df.loc[mask, 'time_CDT'], y.loc[mask], fmt, label=label)
    ax.set_xlabel('Acquisition start (CDT)')
    ax.set_ylabel('Effective attenuation time [ms]')
    ax.set_title('DUNE 2x2 July 2024: frozen IO5 vs compact-group IO5 and IO6')
    ax.grid(alpha=.25)
    ax.legend(ncol=2, fontsize=9)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(dest, dpi=200)
    plt.close(fig)


def main():
    a = args()
    if a.bootstrap < 2 and not a.list_only:
        raise SystemExit('--bootstrap must be at least 2')
    ext_ios = sorted(set(a.ext_io or [5,6]))
    roots = a.frozen_root or [Path('frozen_brooke10_bs20'), Path('frozen_hit_timeseries_bs20')]
    runs = core.collect_runs(roots, a.flow_root)
    if a.only:
        runs = [r for r in runs if a.only in Path(r['packet_file']).name]
    if not runs:
        raise SystemExit('No frozen runs found')
    print(f'Found {len(runs)} runs; new-method EXT selections: {ext_ios}', flush=True)
    for r in runs:
        name = Path(r['packet_file']).name
        print(f'  {name}: frozen EXT IO{r["ext_io"]}, '
              f'raw={Path(r["packet_file"]).is_file()}, '
              f'FLOW={bool(r["flow_path"] and r["flow_path"].is_file())}', flush=True)
    if a.list_only:
        return
    a.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    mpv_rows = []
    for idx, r in enumerate(runs, 1):
        name = Path(r['packet_file']).name
        print(f'\n[{idx}/{len(runs)}] {name}', flush=True)
        row = {'packet_file': name, 'time_CDT': r['timestamp'],
               'frozen_source': r['frozen_source'],
               'frozen_ext_io': r['ext_io'],
               'frozen_unweighted_ms': r['frozen_unweighted_ms'],
               'frozen_weighted_ms': r['frozen_weighted_ms']}
        if not r['flow_path'] or not r['flow_path'].is_file():
            row['problem'] = 'Missing FLOW file'; rows.append(row)
            save_outputs(rows, mpv_rows, a.outdir)
            continue
        for io in ext_ios:
            try:
                if io == r['ext_io']:
                    npz = r['npz']
                    print(f'  new-method IO{io}: reuse frozen selection {npz}', flush=True)
                else:
                    # All previously listed frozen runs use IO5; reconstruct IO6
                    # from raw packets rather than mislabeling an IO5 selection.
                    npz = a.outdir / 'selected_cache' / f'{Path(name).stem}.IO{io}.npz'
                    if a.force_io6_extract or not npz.is_file():
                        extract_npz_from_raw(r['packet_file'], npz, io,
                                             a.ext_trigger_type, a.chunk)
                    else:
                        print(f'  new-method IO{io}: reuse cached {npz}', flush=True)
                if io == 6 and a.validate_io6_npz is not None:
                    check_original_io6(npz, a.validate_io6_npz)
                groups, n_all, n_good, coverage = core.load_and_group(npz, r['flow_path'])
                row[f'new_io{io}_selected_packets'] = n_all
                row[f'new_io{io}_matched_packets'] = n_good
                row[f'new_io{io}_pedestal_coverage'] = coverage
                for suffix, cond in [ ('n2', lambda x: x.n == 2),
                                      ('nge2', lambda x: x.n >= 2) ]:
                    sample = f'new_io{io}_{suffix}'
                    try:
                        seed = zlib.crc32(f'{name}|IO{io}|{suffix}'.encode())
                        d, result = core.fit_one_sample(groups, cond,
                            np.random.default_rng(seed), a.bootstrap)
                        row[f'{sample}_tau_ms'] = result['tau_us'] / 1000
                        row[f'{sample}_tau_err_ms'] = result['tau_err_us'] / 1000
                        row[f'{sample}_chi2_ndf'] = result['chi2'] / result['ndf']
                        row[f'{sample}_n_groups'] = result['n_groups']
                        mpv_rows.extend([{'packet_file': name, 'time_CDT': r['timestamp'],
                                          'sample': sample, **record}
                                         for record in d.to_dict('records')])
                        print(f'    {sample}: tau={row[f"{sample}_tau_ms"]:.3f} '
                              f'+/- {row[f"{sample}_tau_err_ms"]:.3f} ms, '
                              f'chi2/ndf={row[f"{sample}_chi2_ndf"]:.2f}', flush=True)
                    except Exception as exc:
                        row[f'{sample}_failure'] = f'{type(exc).__name__}: {exc}'
                        print(f'    {sample}: NO FIT: {exc}', flush=True)
            except Exception as exc:
                row[f'new_io{io}_failure'] = f'{type(exc).__name__}: {exc}'
                print(f'  IO{io} ERROR: {type(exc).__name__}: {exc}', flush=True)
                if io == 6 and a.validate_io6_npz is not None:
                    # A failed explicit exact-match check must stop validation.
                    raise
        rows.append(row)
        save_outputs(rows, mpv_rows, a.outdir)
    df = pd.DataFrame(rows)
    plot_results(df, a.outdir / 'paired_lifetime_timeseries.png')
    print(f'\nSaved {len(df)} timeline rows to {a.outdir.resolve()}', flush=True)
    print('CSV: paired_lifetime_timeseries.csv')
    print('PNG: paired_lifetime_timeseries.png')
    print('MPVs: new_method_mpvs.csv (if any fits succeeded)')


if __name__ == '__main__':
    main()
