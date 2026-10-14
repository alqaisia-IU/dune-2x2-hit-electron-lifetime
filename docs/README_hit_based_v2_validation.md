# Hit-based electron-lifetime V2 diagnostic (DUNE 2x2)

This new script **does not change** `packet_lifetime_final.py` or the frozen scripts. It was developed from the supplied `run_paired_lifetime_timeline_all.py` and `run_io6_n2_split_io.py` (and reuses their existing functions rather than reimplementing the packet selection and charge fitter).

## Install

Copy `run_hit_based_v2_validation.py` to the existing NERSC directory:

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST
ls run_paired_lifetime_timeline_all.py run_io6_n2_split_io.py
python run_hit_based_v2_validation.py --self-test
```

## First, one-file smoke test

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST
python -u run_hit_based_v2_validation.py \
  --only 2024_07_11_12_09_02 \
  --bootstrap 10 \
  --daq-audit \
  --per-file-plots \
  --outdir hit_based_v2_smoke
```

Default `--reuse-n2-csv` points to `io6_n2_split_full/io6_n2_per_io_mpvs.csv`. If present, the previous n=2 30-bootstrap MPV estimates are reused and **only n=3** MPVs are newly fitted. The smoke-test flag `--bootstrap 10` applies only to newly fit spectra. Do not interpret the smoke test as a finished statistical comparison because the two bootstrap sizes differ.

## Full 18-file diagnostic (30 bootstrap resamples)

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST
python -u run_hit_based_v2_validation.py \
  --bootstrap 30 \
  --daq-audit \
  --per-file-plots \
  --outdir hit_based_v2_full \
  2>&1 | tee hit_based_v2_full.log
```

Run from an appropriate NERSC compute allocation; this may take a while. Files are written incrementally. Running again in the same output directory will recompute and overwrite V2 outputs; it does not modify the original analysis outputs.

To recompute both n=2 and n=3 instead of reusing existing n=2 MPVs: add `--refit-n2`.

For an alternative **within-bin** EXT-cluster bootstrap: use `--refit-n2 --bootstrap-unit ext`. Note that this still does not propagate cross-IO or cross-time-bin EXT correlations into the slope-fit covariance; it is a sensitivity check, not a complete event-level bootstrap.

## Exact physics choices kept unchanged

- Reference time: **external trigger from IO group 6** (from existing IO6 selected-packet caches); charge groups are IO1–IO8.
- Electron drift-time proxy: `t_last` relative to the IO6 EXT, restricted to 20–180 µs.
- Full `(EXT index, IO group, IO channel, ASIC chip, channel)` grouping, **without** successive-gap splitting.
- Exactly **two** or exactly **three** packets per group, separately, with **total first-to-last span ≤10 µs**.
- Static FLOW-v9 channel pedestal subtraction, summed signed charge, four drift bins centered on 40/80/120/160 µs.
- Same `fit_langau(..., bin_width=1.0, qmin=0, qmax=50)` fitter, and same `run_io6_n2_split_io.fit_common_alpha` weighted shared-slope fit.

## Three fits per acquisition

1. `n2`: common `alpha2` for all IO groups, separate amplitudes for charge IO1–IO8 (same previous split-IO baseline).
2. `n3`: common `alpha3` for all IO groups, separate amplitudes for charge IO1–IO8.
3. `shared`: one `alpha` for both n=2 and n=3, but **16 separate normalization parameters** `A_(IO,n)`. No pooling of their charge spectra.

Model: `MPV_(IO,n)(t) = A_(IO,n) * exp[-alpha * (t-40µs)/1000]`, where `alpha` is in inverse ms and a positive alpha corresponds to `tau_ms = 1/alpha`.

The script calculates `delta_chi2 = chi2(shared) - [chi2(n2) + chi2(n3)]` and an **approximate** 1-degree-of-freedom chi-square tail probability only if the same MPVs and amplitude groups enter both nested models. The p-value assumes fixed/independent MPV uncertainties and does **not** account for shared-EXT correlations or detector systematics. Interpret the tau and p-values as diagnostic, not yet publication-grade.

## Files produced

- `v2_lifetime_timeline.csv` — for each acquisition: selected group counts, n2/n3/shared slopes and lifetimes, errors, chi2/ndf, and nested compatibility statistics.
- `v2_per_io_multiplicity_mpvs.csv` — charge-spectrum MPVs/errors, fit quality, fit failures, and provenance for all `(IO,n,drift-bin)`.
- `v2_n2_n3_shared_lifetimes.png` — time-series comparison for positive lifetimes.
- `packet-..._normalized_mpvs.png` — optional per-acquisition illustrative MPV attenuation curves.
- `v2_daq_channel_audit.csv` — optional per-file number of accepted groups whose channels are configured/enabled; **audit only; no filtering**.

## Important limitations

- This is a **comparison of grouped-packet populations**, not a confirmed physical classification of distinct charge deposits. Whether compact same-channel n=3 should be summed at all still needs CRS expert input.
- No n=1 groups are introduced. Simulation threshold maps are **not** applied as cuts or ADC corrections; calibration provenance and threshold/ADC relation remain unresolved.
- Errors come from finite per-spectrum bootstraps and a weighted least-squares covariance, not a full correlated systematic or full-run EXT-level bootstrap.
- For the fastest first run, reuse the existing n=2 CSV; checking that the resulting n2 curve matches the old split-IO baseline is an important cross-check.
