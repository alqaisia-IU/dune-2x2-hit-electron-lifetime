# DUNE 2x2 — Hit-Based Electron Lifetime

Exploratory electron lifetime measurement using charge packets from
the DUNE 2x2 liquid argon time projection chamber (LArTPC).

## Current analysis

**Current version:** Hit-Based V2 (October 2026)

Main script:

`src/run_hit_based_v2_validation.py`

The analysis:

- Uses IO6 external-trigger (EXT) timing as a reference.
- Subtracts static, channel-dependent FLOW-v9 pedestals.
- Groups packets by EXT and electronics channel.
- Selects exactly n=2 or n=3 packets within a 10-us total span.
- Uses the last packet timestamp as the drift-time proxy.
- Fits Landau-Gaussian charge distributions in four drift-time bins.
- Fits MPV attenuation with independent IO-group and multiplicity amplitudes.
- Compares n=2, n=3, and shared n=2+n=3 attenuation lifetimes.

## Results

- `results/v2/` — current numerical lifetime and MPV results.
- `results/brooke/` — comparison plots with Brooke Russell.
- `results/legacy/` — earlier n=2 split-IO analysis.

The Brooke reference is approximately digitized, with assumed
timestamps. It is not a verified run-matched comparison.

## Documentation

See [Full physics and presentation README](docs/README_HIT_BASED_V2_CURRENT_PRESENTATION.md).

Additional documentation:

- `docs/README_hit_based_v2_validation.md`
- `docs/README_V2_vs_Brooke.md`

## Environment

Python dependencies are listed in `requirements.txt`.

The numerical analysis was performed at NERSC using the
`e-lifetime` conda environment.

## Data availability

Raw packet HDF5 files, FLOW files, and selected-packet caches
are not distributed in this repository.

A complete rerun requires the relevant DUNE 2x2 datasets and
cached selection metadata on NERSC.

## Current limitations

The physical interpretation of IO6 EXT as interaction t0,
the precise same-channel retrigger grouping, and the treatment
of above-threshold singletons are still under investigation.

No experimental or simulation discriminator-threshold correction
is currently applied to the lifetime measurement.

This repository preserves a preliminary method-development
analysis, not a finalized calibrated electron lifetime result.

## Repository source organization

### Current analysis — `src/`

- `run_hit_based_v2_validation.py` — Current V2 electron lifetime analysis.

### Previous and supporting scripts — `old/`

Contains the earlier lifetime implementations, split-IO studies,
packet-processing dependencies, and Brooke comparison plotting code.

These scripts are retained for reproducibility.

The current V2 driver automatically imports its required historical
analysis dependencies from `old/`.

To run the V2 mathematical self-test from the repository root:

```bash
python src/run_hit_based_v2_validation.py --self-test
Full detector-data processing requires the appropriate NERSC environment,
FLOW data, and selected-packet caches.


## Brooke Russell — Electron Lifetime Reference

**Presentation:** [ND Prototypes — Electron Lifetime,
Brooke Russell (May 1, 2025)](
https://docs.dunescience.org/cgi-bin/sso/RetrieveFile?docid=33816&filename=NDPrototypes_electron_lifetime___Russell___1May2025.pdf&version=1)

**Comparison data:** [`references/brooke_digitized_approx.csv`](references/brooke_digitized_approx.csv)

The CSV contains approximately digitized electron lifetime
values used for comparison with the current hit-based analysis.

The digitization and assumed timestamps are approximate.
The reference presentation is provided for context; the exact
source of each digitized point has not been independently verified.
