# DUNE 2x2 — Hit-Based Electron Lifetime

Exploratory electron lifetime measurement using charge packets from the DUNE 2x2 liquid argon time projection chamber (LArTPC).

## Current analysis

## Command to run the script

At NERSC, using the existing working directory and `e-lifetime` environment:

```bash
cd /pscratch/sd/a/alqaisia/hit_based_LAST

conda activate e-lifetime

python -u run_hit_based_v2_validation.py \
  --bootstrap 30 \
  --daq-audit \
  --per-file-plots \
  --outdir hit_based_v2_full
Results
Environment
Python dependencies are listed in requirements.txt.
The numerical analysis was performed at NERSC using the e-lifetime conda environment.
Data availability
Raw packet HDF5 files, FLOW files, and selected-packet caches are not distributed in this repository.
A complete rerun requires the relevant DUNE 2x2 datasets and cached selection metadata on NERSC.
Current limitations
The physical interpretation of IO6 EXT as interaction t0, the precise same-channel retrigger grouping, and the treatment of above-threshold singletons are still under investigation.
No experimental or simulation discriminator-threshold correction is currently applied to the lifetime measurement.
This repository preserves a preliminary method-development analysis, not a finalized calibrated electron lifetime result.
