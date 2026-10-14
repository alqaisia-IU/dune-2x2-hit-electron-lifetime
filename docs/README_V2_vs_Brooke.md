# DUNE 2x2 Hit-based V2 versus Brooke Russell combined result

The plotting script overlays the three already-computed V2 lifetime time series with Brooke Russell's **2x2 Demonstrator July 2024** lifetime plot from **slide 4**, *ND Prototypes Electron Lifetime*, May 16, 2025.

## Provenance and warning

The included `brooke_digitized_approx.csv` contains **nine approximate slide digitizations**, with **assumed acquisition timestamps** inherited from a previous exploratory overlay. Brooke's underlying timestamps, run/file identifiers, and point uncertainties were not obtained. These are **not independently verified run matches**. A close-looking overlay can establish suggestive time-trend agreement, **not** point-by-point statistical agreement. The digitized points are shown as open diamonds, with NO invented Brooke error bars. Brooke's hit selection and grouping differ from the V2 exact n=2/n=3 population definitions.

## Install at NERSC

Copy `plot_v2_vs_brooke_combined.py` and `brooke_digitized_approx.csv` into:

`/pscratch/sd/a/alqaisia/hit_based_LAST/`

Then, in your existing `e-lifetime` environment:

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST
python plot_v2_vs_brooke_combined.py \
  --v2 hit_based_v2_full/v2_lifetime_timeline.csv \
  --outdir hit_based_v2_full/brooke_overlay
```

If the results are in a different directory (e.g. the one-file `hit_based_v2_smoke`), change `--v2`. For a meaningful full comparison, prefer the completed `hit_based_v2_full` (18 files, n3 bootstrap 30), not the smoke test with mixed bootstrap counts.

## Output

- `v2_vs_brooke_combined.png`: 2-panel comparison, focused n2/shared top, n3 diagnostic bottom.
- `v2_vs_brooke_focus.png`: n2/shared/Brooke overlay suitable for slides.
- `v2_vs_brooke_assumed_pairs.csv`: all nine Brooke reference rows aligned to V2 values **only at the assumed timestamps**, with illustrative signed differences (not significance tests).

No packet selection, charge fit, or lifetime fit is rerun. `packet_lifetime_final.py` and `run_hit_based_v2_validation.py` remain untouched.
