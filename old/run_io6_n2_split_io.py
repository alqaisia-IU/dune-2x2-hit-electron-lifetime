#!/usr/bin/env python3
"""Control test: fixed IO6 EXT, compact n=2, but separate charge-IO MPVs.

Reuses the *exact* previously tested cache and group construction from
run_paired_lifetime_timeline_all.py.  Only charge-IO pooling is changed:
  old: one Langau MPV per time slice across charge IO groups
  new: one Langau MPV per (charge IO group, time slice), then one shared
       attenuation slope alpha with separate IO normalization A_io.

Read-only with respect to the frozen methods and earlier paired timeline.
Run `python run_io6_n2_split_io.py --help` for instructions.
"""
from __future__ import annotations

import argparse
import json
import re
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

TIME_RANGES = [(20, 60), (60, 100), (100, 140), (140, 180)]
TIME_PIVOT_US = 40.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--cached-root', type=Path, default=Path('paired_lifetime_full_io5_io6/selected_cache'),
                   help='Previous raw-IO6-selection cache; files end in .IO6.npz')
    p.add_argument('--paired-csv', type=Path, default=Path('paired_lifetime_full_io5_io6/paired_lifetime_timeseries.csv'),
                   help='Existing results for matching pooled IO6 n=2 comparisons')
    p.add_argument('--outdir', type=Path, default=Path('io6_n2_split_charge_io_40us_bs30'))
    p.add_argument('--bootstrap', type=int, default=30,
                   help='Per-IO/time-bin charge-group resampling iterations')
    p.add_argument('--min-bin-entries', type=int, default=100,
                   help='Same minimum as the new pooled four-bin analysis')
    p.add_argument('--min-io-bins', type=int, default=4,
                   help='At least this many valid MPV time bins per charge IO group (4 keeps complete 40-us coverage)')
    p.add_argument('--only', type=str, default=None,
                   help='Run filename substring, e.g. 0050018-2024_07_10_09_36_12')
    p.add_argument('--list-only', action='store_true')
    p.add_argument('--rebuild-missing', action='store_true',
                   help='If missing a cached IO6 NPZ, reconstruct it from raw packets (slow)')
    p.add_argument('--brooke-csv', type=Path, default=None,
                   help='OPTIONAL approximate digitized Brooke points, matched to assumed acquisition starts')
    p.add_argument('--per-file-plots', action='store_true',
                   help='Save a per-file MPV plot for detailed checks')
    p.add_argument('--self-test', action='store_true',
                   help='Run mathematical/selection self-tests; does not require detector files/pylandau')
    p.add_argument('--flow-root', type=Path, default=Path('/global/cfs/cdirs/dune/www/data/2x2/reflows/v9/flow/beam'))
    p.add_argument('--frozen-root', type=Path, action='append', default=None,
                   help='Repeat to override discovery folders; default is both frozen timeline folders')
    return p.parse_args()


def selected_compact_n2(groups):
    """An EXACT selection copy: group all hits sharing EXT and channel first;
    then require exactly 2 total packets and a total first-to-last span <=10 us.
    No 10-us successive-gap splitting is performed.
    """
    return groups.loc[(groups['n'] == 2) & (groups['span_us'] <= 10.0) &
                      (groups['t_last'] >= 20.0) & (groups['t_last'] < 180.0)].copy()


def fit_common_alpha(mpvs, min_io_bins=3):
    """Weighted shared-alpha, separate charge-IO amplitudes at 40 us.

    Model: MPV_io(t) = A_io,40 * exp[-alpha_ms^-1*(t_us-40)/1000].
    MPV bootstrap standard deviations treated as KNOWN/independent;
    parameter errors conditional, not total physical lifetime errors.
    """
    d = mpvs.copy()
    d = d[np.isfinite(d['MPV_ADC']) & np.isfinite(d['MPV_error']) &
          (d['MPV_ADC'] > 0) & (d['MPV_error'] > 0)].copy()
    counts = d.groupby('io_group', sort=True).size()
    valid_ios = sorted(int(io) for io in counts.index if counts.loc[io] >= min_io_bins)
    d = d[d['io_group'].isin(valid_ios)].copy()
    if not valid_ios:
        raise RuntimeError('No charge IO groups have enough valid MPV/time bins')
    d = d.sort_values(['io_group', 't_center'])
    io_to_idx = {io: i for i, io in enumerate(valid_ios)}
    index = d['io_group'].map(io_to_idx).to_numpy(dtype=int)
    t = d['t_center'].to_numpy(float)
    q = d['MPV_ADC'].to_numpy(float)
    sig = d['MPV_error'].to_numpy(float)
    n_par = 1 + len(valid_ios)
    if len(q) <= n_par:
        raise RuntimeError(f'Not enough MPVs for slope and {len(valid_ios)} IO normalizations')

    amps0 = [float(d.loc[d.io_group == io, 'MPV_ADC'].median()) for io in valid_ios]

    def predict(p):
        return p[1:][index] * np.exp(-p[0] * (t - TIME_PIVOT_US) / 1000.0)

    def resid(p):
        return (q - predict(p)) / sig

    opt = least_squares(resid, x0=[0.8, *amps0],
                        bounds=([-5.0, *([0.0] * len(valid_ios))],
                                [5.0, *([200.0] * len(valid_ios))]),
                        max_nfev=20000, x_scale='jac')
    if not opt.success:
        raise RuntimeError(f'Common-slope fit did not converge: {opt.message}')
    jac = opt.jac
    jtj = jac.T @ jac
    if np.linalg.matrix_rank(jtj) != n_par:
        raise RuntimeError('Singular common-fit covariance matrix')
    covariance = np.linalg.inv(jtj)  # absolute_sigma=True analogue
    if not np.all(np.isfinite(covariance)):
        raise RuntimeError('Non-finite common-fit covariance')
    alpha = float(opt.x[0])
    alpha_err = float(np.sqrt(max(covariance[0, 0], 0)))
    chi2 = float(np.sum(resid(opt.x)**2))
    ndf = len(q) - n_par
    result = dict(
        alpha_per_ms=alpha, alpha_error_per_ms=alpha_err,
        tau_ms=(1.0/alpha if alpha > 0 else np.nan),
        tau_error_ms=(alpha_err/(alpha*alpha) if alpha > 0 else np.nan),
        chi2=chi2, ndf=ndf, chi2_ndf=chi2/ndf,
        n_io_groups=len(valid_ios), n_mpvs=len(d),
        io_groups=','.join(map(str, valid_ios)),
        n_accepted_spectral_fits=int(len(d)),
        n_poor_spectral_fits=int(np.sum(d['chi2_ndf'].to_numpy(float) > 3)),
    )
    amps = pd.DataFrame([dict(io_group=io, amplitude_at_40us_ADC=float(opt.x[i+1]),
                               amplitude_error_ADC=float(np.sqrt(max(covariance[i+1, i+1],0))),
                               n_time_bins=int(counts.loc[io]))
                         for i, io in enumerate(valid_ios)])
    accepted = d.copy()
    accepted['model_ADC'] = predict(opt.x)
    accepted['pull'] = (accepted['MPV_ADC'] - accepted['model_ADC']) / accepted['MPV_error']
    return result, amps, accepted


def fit_per_io_spectra(groups, filename, core, n_boot=30, min_bin=100):
    selected = selected_compact_n2(groups)
    if selected.empty:
        raise RuntimeError('No compact n=2 channel groups in 20-180 us')
    rows = []
    for io, d in selected.groupby('io_group', sort=True):
        for left, right in TIME_RANGES:
            q = d.loc[(d.t_last >= left) & (d.t_last < right), 'q_sum_ADC'].to_numpy(float)
            record = dict(io_group=int(io), t_start=left, t_end=right,
                          t_center=(left+right)/2, N=int(len(q)),
                          MPV_ADC=np.nan, MPV_error=np.nan, chi2_ndf=np.nan,
                          n_boot_success=0, fit_fraction=np.nan, failure='')
            if len(q) < min_bin:
                record['failure'] = f'Only {len(q)} groups; require {min_bin}'
                rows.append(record)
                continue
            record['fit_fraction'] = float(np.mean((q >= 0) & (q < 50)))
            try:
                central = core.fit_langau(q, bin_width=1., qmin=0., qmax=50.)
                record['MPV_ADC'] = float(central['mpv'])
                record['chi2_ndf'] = float(central['red_chi2'])
            except Exception as exc:
                record['failure'] = f'Central Langau: {type(exc).__name__}: {exc}'
                rows.append(record)
                continue
            seed = zlib.crc32(f'{filename}|charge-IO{int(io)}|{left}-{right}|n2'.encode())
            rng = np.random.default_rng(seed)
            boot = []
            for _ in range(n_boot):
                trial = q[rng.integers(0, len(q), size=len(q))]
                try:
                    mpv = float(core.fit_langau(trial, bin_width=1., qmin=0., qmax=50.)['mpv'])
                    if np.isfinite(mpv):
                        boot.append(mpv)
                except Exception:
                    pass
            record['n_boot_success'] = len(boot)
            if len(boot) >= max(2, (n_boot + 1)//2):
                record['MPV_error'] = float(np.std(boot, ddof=1))
            else:
                record['failure'] = f'Only {len(boot)}/{n_boot} bootstrap fits succeeded'
            if not np.isfinite(record['MPV_error']) or record['MPV_error'] <= 0:
                if not record['failure']:
                    record['failure'] = 'Bootstrap MPV uncertainty non-positive/non-finite'
            rows.append(record)
    return pd.DataFrame(rows), int(len(selected))


def save_incremental(timeline, mpvs, amps, outdir):
    if timeline:
        pd.DataFrame(timeline).sort_values('time_CDT').to_csv(outdir/'io6_n2_split_timeline.csv', index=False)
    if mpvs:
        pd.concat(mpvs, ignore_index=True).to_csv(outdir/'io6_n2_per_io_mpvs.csv', index=False)
    if amps:
        pd.concat(amps, ignore_index=True).to_csv(outdir/'io6_n2_io_normalizations.csv', index=False)


def plot_timeline(result, outfile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = result.sort_values('time_CDT').copy()
    d['time_CDT'] = pd.to_datetime(d['time_CDT'])
    fig, ax = plt.subplots(figsize=(12,6))
    for col, err, lbl, fmt in [
        ('pooled_io6_n2_tau_ms','pooled_io6_n2_tau_err_ms','Current pooled IO6 compact n=2', 'o-'),
        ('split_io6_n2_tau_ms','split_io6_n2_tau_err_ms','New: separate charge-IO MPVs, common slope', 's-')]:
        if col not in d:
            continue
        y = pd.to_numeric(d[col], errors='coerce')
        yy = pd.to_numeric(d[err], errors='coerce')
        okay = np.isfinite(y) & np.isfinite(yy) & (yy >= 0)
        if okay.any():
            ax.errorbar(d.loc[okay, 'time_CDT'], y[okay], yerr=yy[okay], fmt=fmt,
                        capsize=2.5, label=lbl)
    ax.set_title('DUNE 2x2: control test of pooling charge IO groups\n'
                 'IO6 EXT | compact n=2 | span <= 10 us | 4 x 40-us slices')
    ax.set_xlabel('Acquisition start (CDT)')
    ax.set_ylabel('Effective attenuation timescale [ms]')
    ax.legend()
    ax.grid(alpha=.25)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def plot_per_file(accepted, amps, result, outfile):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9,6))
    alpha = result['alpha_per_ms']
    xx = np.linspace(20,180,300)
    for a in amps.itertuples():
        d = accepted[accepted['io_group'] == a.io_group]
        line = ax.errorbar(d.t_center, d.MPV_ADC, yerr=d.MPV_error, fmt='o', capsize=2.5,
                           label=f'Charge IO{a.io_group} (N={len(d)})')
        color = line.lines[0].get_color()
        ax.plot(xx, a.amplitude_at_40us_ADC * np.exp(-alpha*(xx-TIME_PIVOT_US)/1000),
                color=color, alpha=.6)
    ax.set_title(f"IO6 compact n=2 | separate charge-IO fit\n"
                 f"tau={result['tau_ms']:.3f} +/- {result['tau_error_ms']:.3f} ms, "
                 f"chi2/ndf={result['chi2_ndf']:.2f}")
    ax.set_xlabel('Last packet time relative to EXT [us]')
    ax.set_ylabel('Landau-Gaussian fitted MPV [ADC]')
    ax.legend(ncol=2, fontsize=8)
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def compare_brooke(result, brooke_path, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b = pd.read_csv(brooke_path)
    need = {'acquisition_start_cdt','brooke_lifetime_ms'}
    if not need.issubset(b.columns):
        raise ValueError(f'Brooke CSV missing columns: {sorted(need - set(b.columns))}')
    b['time_CDT'] = pd.to_datetime(b['acquisition_start_cdt'])
    d = result.copy()
    d['time_CDT'] = pd.to_datetime(d['time_CDT'])
    m = b.merge(d, on='time_CDT', how='left', validate='one_to_one')
    for col in ['pooled_io6_n2_tau_ms', 'split_io6_n2_tau_ms']:
        m[f'{col}_minus_brooke_ms'] = pd.to_numeric(m[col], errors='coerce') - m['brooke_lifetime_ms']
    m.to_csv(outdir/'brooke_approx_paired_residuals.csv',index=False)
    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(m.time_CDT, m.brooke_lifetime_ms, 'o-', label='Brooke slide: APPROX. digitized/ASSUMED times')
    for col, err, lbl, fmt in [
        ('pooled_io6_n2_tau_ms', 'pooled_io6_n2_tau_err_ms', 'Pooled IO6 n=2','s-'),
        ('split_io6_n2_tau_ms', 'split_io6_n2_tau_err_ms', 'Separate charge-IO, common slope','d-')]:
        y = pd.to_numeric(m[col], errors='coerce')
        e = pd.to_numeric(m[err], errors='coerce')
        okay = np.isfinite(y) & np.isfinite(e)
        if okay.any():
            ax.errorbar(m.loc[okay,'time_CDT'],y[okay],yerr=e[okay],fmt=fmt,capsize=2,label=lbl)
    ax.set_title('Exploratory Brooke comparison: reference timestamps assumed, NOT verified')
    ax.set_xlabel('Assumed acquisition time (CDT)')
    ax.set_ylabel('Effective attenuation timescale [ms]')
    ax.grid(alpha=.25)
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(outdir/'brooke_approx_pooled_vs_split.png', dpi=180)
    plt.close(fig)
    for col in ['pooled_io6_n2_tau_ms','split_io6_n2_tau_ms']:
        dd = m[[col, 'brooke_lifetime_ms']].dropna()
        if len(dd):
            print(f"Brooke APPROX. comparison {col}: N={len(dd)}, "
                  f"MAE={float(np.abs(dd[col]-dd['brooke_lifetime_ms']).mean()):.3f} ms")


def self_test():
    true_alpha = 0.8
    t = np.tile(np.array([40.,80.,120.,160.]),2)
    io = np.repeat([1,2],4)
    a = np.repeat([20., 34.],4)
    q = a*np.exp(-true_alpha*(t-TIME_PIVOT_US)/1000)
    toy = pd.DataFrame(dict(io_group=io, t_center=t, MPV_ADC=q,
                            MPV_error=np.repeat(.1,len(q)), chi2_ndf=np.ones(len(q))))
    fit, amps, accepted = fit_common_alpha(toy)
    assert np.isclose(fit['alpha_per_ms'],true_alpha,atol=1e-6),fit
    assert np.isclose(fit['tau_ms'],1.25,atol=1e-6),fit
    assert fit['ndf']==5 and fit['n_io_groups']==2
    toy_groups = pd.DataFrame(dict(io_group=[1,1,2,2],n=[2,2,3,2],span_us=[2.8,12,3,8],
                                   t_last=[50,60,80,190]))
    assert list(selected_compact_n2(toy_groups).index)==[0]
    print('SELF-TEST PASSED: common-slope math and compact group selection')


def main():
    a = parse_args()
    if a.self_test:
        self_test()
        return
    if a.bootstrap < 2:
        raise SystemExit('--bootstrap must be >=2')
    # Imported only after --self-test so math can be tested without detector installations.
    import run_paired_lifetime_timeline_all as core
    import run_paired_lifetime_timeline_io5_io6 as io_runner

    roots = a.frozen_root or [Path('frozen_brooke10_bs20'), Path('frozen_hit_timeseries_bs20')]
    runs = core.collect_runs(roots, a.flow_root)
    if a.only:
        runs = [r for r in runs if a.only in Path(r['packet_file']).name]
    if not runs:
        raise SystemExit('No source runs found: check --only and frozen directories')
    baseline = pd.read_csv(a.paired_csv)
    if 'packet_file' not in baseline or 'new_io6_n2_tau_ms' not in baseline:
        raise SystemExit('Existing paired CSV lacks packet_file or new_io6_n2_tau_ms')
    if baseline['packet_file'].duplicated().any():
        raise SystemExit('Duplicate packet_file rows in paired baseline CSV')
    old = baseline.set_index('packet_file')
    print(f'Found {len(runs)} frozen file references and {len(old)} old pooled comparison rows')
    print('Only new split-by-charge-IO analysis will be recomputed. EXT IO6 is fixed.')

    for r in runs:
        name = Path(r['packet_file']).name
        cache = a.cached_root / (Path(name).stem + '.IO6.npz')
        print(f"  {name}: IO6 cache={'OK' if cache.is_file() else 'MISSING'}, "
              f"FLOW={'OK' if r['flow_path'] and r['flow_path'].is_file() else 'MISSING'}, "
              f"pooled={old.loc[name, 'new_io6_n2_tau_ms'] if name in old.index else 'NO ROW'}")
    if a.list_only:
        return
    a.outdir.mkdir(parents=True, exist_ok=True)
    timeline, mpv_parts, amp_parts = [], [], []
    for idx,r in enumerate(runs,1):
        name=Path(r['packet_file']).name
        print(f'\n[{idx}/{len(runs)}] {name}',flush=True)
        row=dict(packet_file=name,time_CDT=r['timestamp'],
                 ext_reference='IO6', selection='compact n=2, total span <=10us, t_last 20-180us',
                 n_raw_selected=np.nan,n_ped_matched=np.nan,ped_coverage=np.nan,
                 n_selected_groups=np.nan,status='FAILED', failure='')
        if name in old.index:
            x=old.loc[name]
            for prior, now in [('new_io6_n2_tau_ms','pooled_io6_n2_tau_ms'),
                               ('new_io6_n2_tau_err_ms','pooled_io6_n2_tau_err_ms'),
                               ('new_io6_n2_chi2_ndf','pooled_io6_n2_chi2_ndf')]:
                row[now]=x.get(prior, np.nan)
        cache=a.cached_root/(Path(name).stem+'.IO6.npz')
        try:
            if not r['flow_path'] or not r['flow_path'].is_file():
                raise FileNotFoundError('Matching reflow v9 FLOW file missing')
            if not cache.is_file():
                if not a.rebuild_missing:
                    raise FileNotFoundError(f'IO6 cache missing: {cache}; rerun earlier timeline or pass --rebuild-missing')
                print('  IO6 cache missing. Rebuilding from original raw HDF5...',flush=True)
                io_runner.extract_npz_from_raw(r['packet_file'],cache,6,2,2_000_000)
            groups, n_all, n_good, coverage=core.load_and_group(cache,r['flow_path'])
            row.update(n_raw_selected=n_all,n_ped_matched=n_good,ped_coverage=coverage)
            df, n_selected=fit_per_io_spectra(groups,name,core,a.bootstrap,a.min_bin_entries)
            row['n_selected_groups']=n_selected
            df['packet_file']=name
            df['time_CDT']=r['timestamp']
            mpv_parts.append(df)
            print(f'  n=2 compact selected groups: {n_selected}; IO/time bin fits: '
                  f'{len(df)}, failed: {int((df.failure != "").sum())}',flush=True)
            result, amps, accepted=fit_common_alpha(df, min_io_bins=a.min_io_bins)
            row.update(split_io6_n2_tau_ms=result['tau_ms'],
                       split_io6_n2_tau_err_ms=result['tau_error_ms'],
                       split_io6_n2_alpha_per_ms=result['alpha_per_ms'],
                       split_io6_n2_alpha_error_per_ms=result['alpha_error_per_ms'],
                       split_io6_n2_chi2_ndf=result['chi2_ndf'],
                       split_io6_n2_ndf=result['ndf'],
                       split_io6_n2_n_io_groups=result['n_io_groups'],
                       split_io6_n2_n_mpvs=result['n_mpvs'],
                       split_io6_n2_n_bad_langau=result['n_poor_spectral_fits'],
                       split_io6_n2_used_io_groups=result['io_groups'],status='OK')
            amps['packet_file']=name
            amps['time_CDT']=r['timestamp']
            amp_parts.append(amps)
            print(f"  tau={result['tau_ms']:.3f} +/- {result['tau_error_ms']:.3f} ms; "
                  f"chi2/ndf={result['chi2_ndf']:.2f}; charge IO groups="
                  f"{result['io_groups']}; high-Langau-chi2={result['n_poor_spectral_fits']}/{result['n_mpvs']}",flush=True)
            if a.per_file_plots:
                safe=Path(name).stem
                plot_per_file(accepted,amps,result,a.outdir/f'{safe}_charge_io_MPVs.png')
        except Exception as exc:
            row['failure']=f'{type(exc).__name__}: {exc}'
            print(f'  FAILED: {row["failure"]}',flush=True)
        timeline.append(row)
        save_incremental(timeline,mpv_parts,amp_parts,a.outdir)
    result=pd.DataFrame(timeline).sort_values('time_CDT')
    plot_timeline(result,a.outdir/'io6_n2_pooled_vs_split.png')
    if a.brooke_csv is not None:
        compare_brooke(result,a.brooke_csv,a.outdir)
    print(f'\nSaved {int((result.status=="OK").sum())}/{len(result)} successful files to {a.outdir.resolve()}')
    print('Primary plot: io6_n2_pooled_vs_split.png')
    print('Primary CSV: io6_n2_split_timeline.csv')
    print('Per-IO MPVs: io6_n2_per_io_mpvs.csv')


if __name__ == '__main__':
    main()
