# DUNE 2×2 electron lifetime — current hit-based V2 analysis

**Presentation and reproducibility README — 13 October 2026**  
**Current V2 driver:** `run_hit_based_v2_validation.py`  
**Primary NERSC working directory:** `/pscratch/sd/a/alqaisia/hit_based_LAST`  
**Status:** exploratory, reproducible *effective charge-attenuation-timescale* analysis; not yet a validated detector-purity or Brooke-method reproduction.

> **Short description for a talk.** Starting from raw DUNE 2×2 charge packets previously associated with an IO-group-6 external-trigger (EXT) stream, we subtract a static per-electronics-channel pedestal derived from matching FLOW-v9 data, sum short groups of packets on the same electronic channel, and use the last packet time relative to EXT as the drift-time proxy. We study **exactly n=2** and **exactly n=3** packets separately. A Landau–Gaussian fit gives the charge MPV in four time bins; an exponential fit shares a slope across IO groups while allowing their charge normalizations to differ. V2 additionally tests whether n=2 and n=3 are compatible with one common slope.

## 1. Exactly which scripts to keep locally

**Source files: preserve the exact versions currently on NERSC**, not just renamed or older exploratory notebooks.

| File | Role | Essential? |
|---|---|---|
| `run_hit_based_v2_validation.py` | **CURRENT V2 DRIVER.** Generates n2, n3 and shared-slope results, optional DAQ audit and plots | **Yes** |
| `run_paired_lifetime_timeline_all.py` | V2's input discovery, loading cached IO6 selections, FLOW pedestal subtraction, exact EXT/channel grouping, Langau fitter import | **Yes** |
| `run_io6_n2_split_io.py` | V2's weighted common-exponential-slope fitter with independent amplitudes; earlier n2 split-IO control | **Yes** |
| `packet_lifetime_final.py` | Provides `flow_ped_map`, channel `key`, `KE_PER_ADC` to the helper; older method itself is **not** V2 | **Yes (dependency)** |
| `packet_lifetime_v3_brooke.py` | Provides the `fit_langau` Landau–Gaussian peak fitter used by the helper | **Yes (dependency)** |
| `plot_v2_vs_brooke_combined.py` | Produces the qualitative overlay plots *from saved V2 CSV*, no refitting | Yes for presentations |
| `brooke_digitized_approx.csv` | Nine approximate Brooke slide points with **assumed** dates; not validated run matches | Yes for reproducing current overlay |
| `packet_lifetime_FROZEN.py` | Original MCLK / cached-IO6-selection supporting code | Keep for provenance / cache regeneration |
| `run_paired_lifetime_timeline_io5_io6.py` | Previous IO5/IO6 timeline runner that generated IO6 selections | Keep for provenance / cache regeneration |
| `README_HIT_BASED_V2_CURRENT_PRESENTATION.md` | This README | **Yes** |

**Important packaging note:** The ChatGPT source snapshot supplied with this README includes the V2, helper, split-IO and plotting scripts and the Brooke reference CSV. **It does not contain the actual NERSC `packet_lifetime_final.py`, `packet_lifetime_v3_brooke.py`, `packet_lifetime_FROZEN.py` or `run_paired_lifetime_timeline_io5_io6.py`.** Copy those files *from NERSC*, where your completed runs already have their actual imported versions. Do not substitute `packet_lifetime(1).py` or a similarly named historical file for `packet_lifetime_final.py`.

### Code snapshot commands to run ON NERSC

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST

# Inspect and ensure the actual imported code is present.
ls -l \
  run_hit_based_v2_validation.py \
  run_paired_lifetime_timeline_all.py \
  run_io6_n2_split_io.py \
  packet_lifetime_final.py \
  packet_lifetime_v3_brooke.py \
  packet_lifetime_FROZEN.py \
  run_paired_lifetime_timeline_io5_io6.py \
  plot_v2_vs_brooke_combined.py \
  brooke_digitized_approx.csv

# Once you have also copied this README into this directory:
mkdir -p CURRENT_V2_CODE_SNAPSHOT
cp -p \
  run_hit_based_v2_validation.py \
  run_paired_lifetime_timeline_all.py \
  run_io6_n2_split_io.py \
  packet_lifetime_final.py \
  packet_lifetime_v3_brooke.py \
  packet_lifetime_FROZEN.py \
  run_paired_lifetime_timeline_io5_io6.py \
  plot_v2_vs_brooke_combined.py \
  brooke_digitized_approx.csv \
  README_HIT_BASED_V2_CURRENT_PRESENTATION.md \
  CURRENT_V2_CODE_SNAPSHOT/

(cd CURRENT_V2_CODE_SNAPSHOT && sha256sum * > SHA256SUMS.txt)
tar -czf CURRENT_V2_CODE_SNAPSHOT_2026-10-13.tar.gz CURRENT_V2_CODE_SNAPSHOT/
ls -lh CURRENT_V2_CODE_SNAPSHOT_2026-10-13.tar.gz
```

This creates a small, portable **source-code snapshot**, not a raw-data backup. Download the `.tar.gz` to your computer using your usual NERSC file-transfer method (working DTN/Globus/scp) and retain the checksums. For a *bit-for-bit reproducible numerical rerun*, preserve the correct input caches/metadata and source data described below as well.

### Files to preserve from the completed analysis results

```text
hit_based_v2_full/v2_lifetime_timeline.csv
hit_based_v2_full/v2_per_io_multiplicity_mpvs.csv
hit_based_v2_full/v2_daq_channel_audit.csv              # if --daq-audit was used
hit_based_v2_full/v2_n2_n3_shared_lifetimes.png
hit_based_v2_full/packet-*_normalized_mpvs.png           # if --per-file-plots was used
hit_based_v2_full/brooke_overlay/v2_vs_brooke_combined.png
hit_based_v2_full/brooke_overlay/v2_vs_brooke_focus.png
hit_based_v2_full/brooke_overlay/v2_vs_brooke_assumed_pairs.csv
io6_n2_split_full/io6_n2_per_io_mpvs.csv                 # IMPORTANT: reused n2 MPVs
hit_based_v2_full.log                                    # if retained
```

**Inputs needed to rerun at NERSC** (not included in the small source ZIP):

- `paired_lifetime_full_io5_io6/selected_cache/packet-....IO6.npz` — preselected raw-charge packets, packet-relative times and EXT indices.
- `frozen_brooke10_bs20/` and `frozen_hit_timeseries_bs20/` — at least the `selected_packets.npz`/`summary.json` file discovery pairs; V2 uses the summaries to find acquisitions. The IO6 selection itself is loaded from the separate IO6 cache, **not** by treating the old IO5 selection as IO6.
- Matching FLOW-v9 charge files under `/global/cfs/cdirs/dune/www/data/2x2/reflows/v9/flow/beam/...` — pedestal recovery.
- Raw `.h5` files under `/global/cfs/cdirs/dune/www/data/2x2/nearline/packet/beam/...` — needed for regenerating caches or running the optional DAQ audit.

A local code snapshot can be read, presented and reviewed without copying these large HDF5 inputs. Running V2 on a laptop will require equivalent local paths/options and dependencies.

## 2. Input data and the role of the external trigger

**Raw data.** DUNE 2×2 LArPix CRS HDF5 charge packets. A charge packet (`packet_type == 0`) carries, among other information, `io_group`, `io_channel`, `chip_id`, `channel_id`, `dataword` (raw ADC), and electronics timing. An external-trigger marker is `packet_type == 7`, `trigger_type == 2` in the tested upstream packet-selection code.

**IO6 EXT reference.** The acquisition's **IO-group-6 external trigger markers** define the reference times for charge packets from **all charge IO groups 1–8**. This is the most important terminology distinction: the **EXT reference source = IO6**, while **charge IO6** is only one of eight charge-electronics populations being fitted. It is not eight separate time-zero sources.

**Clock and packet timing.** The upstream selection unwraps/aligns electronics master-clock timestamps (MCLK; 0.1 µs per tick) and associates a charge packet with the **most recent preceding IO6 EXT marker**. The cached array stores its assigned `ext_index` and time difference, essentially

\[
\Delta t_j = 0.1\;\mu\mathrm{s}\;[\mathrm{MCLK}_{j}-\mathrm{MCLK}_{\mathrm{EXT}(j)}].
\]

The prior IO6 workflow selected charge packets in approximately `0 <= Δt < 200 µs` and stored them in `.IO6.npz`. The current V2 script **loads this cache**; it does not perform a fresh raw-trigger search on every run. V2 checks the cache exists and fails rather than silently substituting another EXT stream.

**Is this the actual physical drift time?** **Not established.** It is a *drift-time proxy* conditional on correct EXT-to-charge association and on the EXT representing the relevant ionization start time (`t0`). The approximately 200 µs span is of the order of the maximum 2×2 TPC drift time (about 30.4 cm maximum drift divided by the drift velocity), but the real trigger origin and association require CRS confirmation. There is no independent reconstructed event `t0` correction in V2.

## 3. Pedestal subtraction and the meaning of charge

**Static pedestal source: matching REFLOW/FLOW v9.** The helper imports `flow_ped_map` from `packet_lifetime_final.py`, recovering the static channel pedestal implicit in FLOW-v9 raw charge calibration (`Q_raw`, with the imported fixed `KE_PER_ADC` convention). This is **not** a pedestal averaged from this packet-file's charge histogram, and **not** the earlier panel-pedestal JSON selection. The pedestal is looked up by full electronics channel `(io_group, io_channel, chip_id, channel_id)`.

For each retained packet,

\[
q_j^{\mathrm{ADC}} = \mathrm{dataword}_j-P_{\mathrm{FLOW},c}.
\]

The **signed** charges are preserved; no `abs(ADC-P)` operation is used. If the group has n packets from channel c, its summed charge is

\[
Q_{\mathrm{group}}^{\mathrm{ADC}}=\sum_{j=1}^{n}(\mathrm{ADC}_j-P_{\mathrm{FLOW},c})
=\sum_j\mathrm{ADC}_j-nP_{\mathrm{FLOW},c}.
\]

For example, with ADC readings 35 and 31 and a pedestal of 28, n=2 gives `(35−28)+(31−28)=10` pedestal-subtracted ADC, **not** `35+31−28=38`.

Missing pedestal entries, non-finite times, and packet times outside `[0,200)` µs are removed by the existing helper; its matched-packet coverage must be ≥95% or the file fails. `pedestal_coverage` is saved in the V2 timeline. **The MPVs are in summed, pedestal-subtracted ADC counts**, not automatically electrons (e⁻), energy deposition, or calibrated `dQ/dx`. The analysis does not apply electron-attachment correction, recombination correction, charge equalization, or track-path-length normalization before these lifetime fits.

## 4. How one provisional hit is constructed

A *group* is all cached packets sharing the same

```text
(ext_index, io_group, io_channel, chip_id, channel_id)
```

The software groups **all packets in that EXT × channel**, not by consecutive intervals separated by a maximum time gap. Only **after grouping** it computes:

```text
n          = number of packets in the full group
t_first    = earliest packet Δt within the group
t_last     = latest packet Δt within the group
span_us    = t_last − t_first
q_sum_ADC  = sum of signed per-packet (ADC − FLOW pedestal)
```

The current selection is:

- `n == 2` **or** `n == 3`, studied as **separate** samples (not n≥2).
- `span_us <= 10` µs, the **total first-to-last span**.
- `20 <= t_last < 180` µs.

For instance, a same-channel three-packet group at 70, 73 and 78 µs is n=3, total span 8 µs, assigned time 78 µs, and accepted in the n=3 population. A four-packet group is not used, however compact. A singleton is not used, however large its measured ADC. A group of three packets with span 11 µs is rejected even if some consecutive pairs are within 10 µs.

**Physics interpretation warning:** These are provisional **same-channel, same-EXT charge groups**, **not necessarily one isolated particle hit**. Repeated self-triggering, integration/retrigger electronics response, threshold bias, pileup, and accidental combinations can influence both n and charge. The 10 µs span and exact-n classification are **our exploratory criteria**, not documented Brooke Russell settings. In particular V2 **does not impose the older empirically observed exact 28-clock-tick spacing (2.8 µs)**.

**Time assigned:** `t_last`, not `t_first`, packet-average time, or reconstructed t0. This is a deliberate part of the provisional method; when discussing physical drift time, say “last packet time relative to IO6 EXT.”

## 5. Four charge spectra and Landau–Gaussian MPVs

For each acquisition, charge IO group **and** packet multiplicity, the code fills four independent summed-charge histograms:

| Relative-time selection | Plotting / fit time center |
|---|---|
| [20, 60) µs | 40 µs |
| [60, 100) µs | 80 µs |
| [100, 140) µs | 120 µs |
| [140, 180) µs | 160 µs |

These are **bin centers**, not individually reconstructed drift times. Up to 8 × 2 × 4 = **64 charge-spectrum MPVs per acquisition** (fewer if a fit fails). V2 requires ≥100 accepted groups for a *newly fitted* spectrum; an IO/multiplicity family needs four valid MPVs by default to enter the exponential fit. If the existing n=2 CSV is reused, the preexisting fitted MPVs/errors are read instead of refitting.

The imported `fit_langau(q, bin_width=1.0, qmin=0, qmax=50)` fits a Landau convolved with a Gaussian (“Langau”):

\[
f(Q)=\mathrm{Landau}(Q)\otimes\mathcal N(Q;\sigma),
\]

with a 1-ADC histogram across **0 ≤ Q < 50 ADC**, retaining only the central histogram region around the smoothed peak (a full-width-at-half-maximum/FWHM-style window). It fits the central part with `curve_fit`, and reports the location of the **maximum of the fitted convolution** as the MPV. It is **not** a fit to the entire high-charge tail and is **not** simply the sample median. The legacy fitter performs unweighted histogram fitting in its core region and calculates a Poisson-like χ² diagnostic afterward; `chi2_ndf > 3` is **counted/flagged**, not an automatic cut in V2.

### MPV uncertainty

Default **30 bootstrap resamples** for newly fitted spectra, resampling charge **groups** with replacement *within each IO × multiplicity × time bin*, refitting Langau and taking the standard deviation of successful bootstrap MPVs. With the default `--reuse-n2-csv`, earlier n=2 results (typically 30 bootstraps) are loaded from `io6_n2_split_full/io6_n2_per_io_mpvs.csv`, while n=3 MPVs are newly bootstrapped. For a sensitivity test, `--refit-n2 --bootstrap-unit ext` resamples EXT clusters within each separate spectrum bin. This is **not** a full acquisition-wide correlated bootstrap.

## 6. The three V2 lifetime fits

V2 performs three independent *fit hypotheses*, with **different amplitudes for charge IO populations**:

1. **`n2`**: eight charge IO groups, each with its own amplitude \(A_{g,2}\), sharing slope \(\alpha_2\).
2. **`n3`**: eight charge IO groups, each with amplitude \(A_{g,3}\), sharing slope \(\alpha_3\).
3. **`shared`**: up to **16 separate** amplitudes \(A_{g,n}\) (one per charge IO × n), all sharing **one** \(\alpha\). This fit combines **MPV points**, **not** the n=2 and n=3 charge histograms.

The fitted charge-MPV model is:

\[
\mathrm{MPV}_{g,n}(t)= A_{g,n}(40\,\mu\mathrm{s})\;\exp\!\left[-\alpha\frac{t-40\,\mu\mathrm{s}}{1000\,\mu\mathrm{s/ms}}\right],\qquad\tau_{\mathrm{ms}}=\frac{1}{\alpha_{\mathrm{ms}^{-1}}}\quad(\alpha>0).
\]

A positive decay slope corresponds to electron attenuation with later drift time; the fitted \(A_{g,n}\) represents the predicted MPV at 40 µs, *not* an absolute unattenuated charge calibration at t=0. The fit minimizes

\[
\chi^2=\sum_{g,n,t}\left[\frac{\mathrm{MPV}_{g,n,t}-\mathrm{model}_{g,n,t}}{\sigma_{\mathrm{MPV},g,n,t}}\right]^2.
\]

`run_io6_n2_split_io.fit_common_alpha()` performs the weighted least-squares fit. It uses its inverse Jacobian normal matrix to estimate parameter uncertainties, and propagates the fitted slope error approximately via \(\sigma_\tau=\sigma_\alpha/\alpha^2\). The reported lifetime is an **effective attenuation time constant**; an actual liquid-argon electron lifetime interpretation additionally requires the signal selection, t0 and response to be physically appropriate.

**Why separate amplitudes?** Charge IO groups can have different pedestal/gain/threshold response, and n=3 summed charges occupy a different scale from n=2. A common slope allows their charge scales to differ without directly pooling unlike charge spectra. **It does not correct for IO-dependent drift-time distortions**.

### Common-slope compatibility check

V2 also reports:

\[
\Delta\chi^2 = \chi^2_{\mathrm{shared}}-(\chi^2_{n2}+\chi^2_{n3}),
\]

and an **approximate χ² tail probability with 1 degree of freedom**, for comparing one shared α versus separate α2 and α3, provided the model's retained sample counts are compatible. It is **diagnostic only**: independent fixed MPV errors, ideal nested-model assumptions, no shared EXT/systematic covariance. Do not describe this p-value as proof that n=2 and n=3 physically measure the same phenomenon.

## 7. What V2 does with thresholds, calibrations and geometry — and what it does NOT

### Present in V2 **today**

- **Pedestal**: yes — the channel-dependent static FLOW-v9 baseline is subtracted **from every individual packet** before summing.
- **Trigger threshold**: **no charge threshold-based correction or acceptance rule.** Selection uses multiplicity, compact span, time, and fitting-range constraints; the 0–50 ADC fit range is *not* the hardware self-trigger threshold.
- **Actual DAQ configuration**: optional `--daq-audit` decompresses `daq_configs` (XZ-compressed TAR of ASIC JSON), checks whether selected n2/n3 channels are configured and whether `channel_mask==0` and `csa_enable==1`, and writes an audit **without rejecting or changing any selected group**.
- **Simulation threshold maps**: **not applied to hits or lifetime fits**. The files `simulation/larndsim_data/2x2_CRS_thresholds/thresholds_module{0,1,2,3}.npz` were studied separately as candidate maps likely in electrons. Their date applicability and ADC-electron/threshold response remain unconfirmed.
- **Per-pixel geometry**: **not required** for lifetime extraction. An exploratory separate study geometrically matched CRS electronics IDs to FLOW pixel (z,y) coordinates, with Fortran-order LUT indexing plus an offset of 1. This is **not** currently used to alter V2 hit charge, timing, or selection.
- **Drift-distance or track correction**: no track reconstruction, track angle, path length, dQ/dx, drift-position reconstruction, or space-charge correction is used in V2.
- **Purity monitor / gas analyzers**: no comparison or constraints enter this fit.

### What we learned but have *not yet adopted*

Actual July 10/11 ASIC configurations were recovered. For the tested example `(IO8,1,22,61)` both periods had `threshold_global=34`, `pixel_trim_dac=2`, unmasked, CSA enabled; for tested IO3/IO4 ASIC configs the July 10→July 11 pixel trims changed on 162 of 101,568 channels, while global thresholds did not differ in the comparison. For July 11 IO3/IO4, direct acquisition settings matched 38,421 of 38,424 tested singletons, resolving the earlier incomplete CSV match. A candidate geometry-to-simulation-module1 threshold mapping was produced for 77,763 pixels and showed a positive ASIC-level association with global DAC settings, but the individual pixel correspondence/provenance is **not a validated operational threshold cut**.

**Crucial interpretation:** `threshold_global=34` is an ASIC DAC value, **not** a pedestal-subtracted packet charge of 34 ADC. Likewise 5,000 electrons in a simulation threshold map is **not** 5,000 ADC nor automatically a valid July 2024 measured cutoff. We cannot yet classify any singleton as “above threshold” solely from these numbers.

## 8. Current observed results and how to present them

For the completed July 2024 comparison (approximately 18 acquisitions in the prior analyses), the exact-n2 lifetime and the shared-n2+n3 lifetime trajectories are **close** on the plots, while n3 alone is noisier and can have much larger fitted errors. This is evidence that including the separate n3 MPV population does not drastically change the shared fitted slope in many acquisitions; **it does not, on its own, establish improved precision or an unbiased physical lifetime**.

The `plot_v2_vs_brooke_combined.py` script reads the already-produced `v2_lifetime_timeline.csv` and overlays Brooke Russell's combined-detector July 2024 slide values. **Important:** Brooke's nine reference points were *approximately digitized*, and the paired timestamps are **assumed**; original measurement timestamps/run IDs and errors were not available. Agreement is **qualitative time-trend agreement only**. The script produces a two-panel diagnostic plot (`v2_vs_brooke_combined.png`) and a simpler talk plot (`v2_vs_brooke_focus.png`); they use the same n2/shared/Brooke numerical series. The two-panel version additionally shows n3-alone.

**Safe talk wording:** “Our exploratory IO6-EXT compact-packet method produces an effective attenuation timescale of order 1 ms over these acquisitions, and the n2 and shared n2+n3 fits show similar time evolution. An approximately digitized Brooke reference exhibits a broadly comparable timescale, but event/trigger definitions and run matching are not yet validated.”

**Avoid claiming:** “We have reproduced Brooke's lifetime exactly,” “the threshold map proves which singleton signals are valid,” “IO6 EXT is definitely the physical interaction t0,” or “adding n3 has been proven statistically superior.”

## 9. Re-running the current code at NERSC

Use the `e-lifetime` conda environment with `numpy`, `scipy`, `pandas`, `h5py`, `matplotlib`, `pylandau` and standard-library Python packages. The two imported custom physics modules must exist in the same directory (see §1).

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST
conda activate e-lifetime

# Verify the V2 selection and weighted-model implementation:
python run_hit_based_v2_validation.py --self-test

# Check input discovery; read-only:
python run_hit_based_v2_validation.py --list-only

# One acquisition:
python -u run_hit_based_v2_validation.py \
  --only 2024_07_11_12_09_02 \
  --bootstrap 30 --daq-audit --per-file-plots \
  --outdir hit_based_v2_single

# All discovered acquisitions (typically 18 in our July 2024 work):
python -u run_hit_based_v2_validation.py \
  --bootstrap 30 --daq-audit --per-file-plots \
  --outdir hit_based_v2_full \
  2>&1 | tee hit_based_v2_full.log

# Reproduce the Brooke *visual* overlay, without refitting the physics:
python plot_v2_vs_brooke_combined.py \
  --v2 hit_based_v2_full/v2_lifetime_timeline.csv \
  --outdir hit_based_v2_full/brooke_overlay
```

**Reused n2 MPVs:** unless `--refit-n2` is supplied, V2 reuses `io6_n2_split_full/io6_n2_per_io_mpvs.csv` if it exists, so the V2 n2 curve is meant to reproduce the earlier split-IO control. To regenerate n2 using the present selection and fitting routine, add `--refit-n2`. `--bootstrap-unit ext` only causes new fitting when paired with `--refit-n2`; the code deliberately does not reuse group-bootstrap MPVs for EXT-bootstrap mode.

**Outputs:** `v2_lifetime_timeline.csv` (three fitted taus, errors, χ² and slope comparison by acquisition), `v2_per_io_multiplicity_mpvs.csv` (per spectrum charge MPV, error, fit quality), optional `v2_daq_channel_audit.csv`, `v2_n2_n3_shared_lifetimes.png`, and optional `packet-*_normalized_mpvs.png`.

**Restore warning:** Source files by themselves will not reproduce the fit without the cached IO6 selected packets, frozen run metadata, matching FLOW-v9 files, and/or the original raw files. NERSC remains the data-backed execution environment; a local snapshot is for preservation and inspection unless those inputs are also copied and paths configured.

## 10. Presentation-ready sequence and open physics questions

A five-slide methods/validation talk can follow this structure:

1. **Aim:** electron attachment in liquid argon, signal attenuation \(Q(t)\propto e^{-t/\tau}\); distinguish electron lifetime from particle lifetime.
2. **Timing and hit definition:** charge type-0 packets; IO6 type-7, trigger-type-2 EXT; MCLK unrolling; \(t_j-t_\mathrm{EXT}\); EXT × channel groups; last-trigger time; compact exact n=2 versus exact n=3.
3. **Charge and statistics:** FLOW-v9 channel pedestal subtraction **per packet**, signed sum; four 40-µs bins; FWHM Langau MPV; bootstrap error.
4. **Lifetime model and V2 checks:** per-IO and per-multiplicity independent amplitudes, shared \(\alpha\), three fit hypotheses, \(\Delta\chi^2\), explain what is/not shared.
5. **July 2024 comparison and caveats:** n2 and shared trends, noisy n3, Brooke slide overlay as qualitative only; unknown t0/EXT interpretation, successive-trigger grouping, singleton hardware-threshold treatment, and incomplete systematics.

**Open CRS/physics questions:** Which EXT stream/time corresponds to the relevant physical `t0`? What exactly constitutes *successive* same-channel retriggers (including n≥3 and deadtime)? Should single-trigger groups above a calibrated discriminator threshold be retained, and how is that threshold related to ADC? How reliable is the static pedestal over time? Can charge spectra change shape or selection efficiency with drift? Would full EXT-level cluster bootstrapping or cross-method PM/track measurements expose bias?

---

### Bottom line

**Present code:** V2 uses **IO6 EXT**, **FLOW-v9 static pedestals**, **signed per-packet ADC subtraction**, **compact exact n2/n3 same-channel groups**, **last packet Δt**, **four time bins**, **Landau–Gaussian MPVs**, and **IO-/multiplicity-separated amplitudes with shared-exponential-slope diagnostics**. **It does not use actual discriminator threshold values to admit singletons**, nor does it apply simulation threshold maps or reconstructed track information. This is a transparent **method-development lifetime measurement**, not yet the final calibrated 2×2 electron lifetime.
