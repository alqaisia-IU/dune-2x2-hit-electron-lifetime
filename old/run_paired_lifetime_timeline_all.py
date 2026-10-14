#!/usr/bin/env python3
"""Extend the 2024 July 2x2 compact-group analysis over frozen-run packet selections.

This NEVER edits the frozen script or existing frozen products. For each existing
selected_packets.npz + summary.json, it reuses the exact saved EXT/time selection,
recovers the FLOW static pedestal, and independently fits two compact-group samples.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Import the exact tested pedestal mapper and existing FWHM Langau fitter.
from packet_lifetime_final import flow_ped_map, key as flow_key, KE_PER_ADC
from packet_lifetime_v3_brooke import fit_langau

CHANNEL_KEYS = ["io_group", "io_channel", "chip_id", "channel_id"]
GROUP_KEYS = ["ext_index"] + CHANNEL_KEYS
TIME_RANGES = [(20, 60), (60, 100), (100, 140), (140, 180)]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--frozen-root", type=Path, action="append", default=None,
        help="Frozen output directory; may be repeated. Default: both frozen_brooke10_bs20 and frozen_hit_timeseries_bs20."
    )
    p.add_argument("--flow-root", type=Path, default=Path(
        "/global/cfs/cdirs/dune/www/data/2x2/reflows/v9/flow/beam"))
    p.add_argument("--outdir", type=Path, default=Path("paired_lifetime_timeline_v9"))
    p.add_argument("--bootstrap", type=int, default=30)
    p.add_argument("--list-only", action="store_true")
    p.add_argument("--only", default=None, help="Limit to packet filenames containing this substring")
    p.add_argument("--allow-non-io6", action="store_true",
                   help="Allow files for which the frozen selection used an EXT IO other than 6")
    return p.parse_args()


def get_flow_path(packet_file, flow_root):
    """Match exact filename and acquisition subdirectory, not a broad search."""
    raw = str(packet_file)
    match = re.search(r"[/\\]packet[/\\]beam[/\\](.+)$", raw)
    if not match:
        return None
    relative = Path(match.group(1))
    if not relative.name.endswith(".h5"):
        return None
    return flow_root / relative.parent / (relative.name[:-3] + ".FLOW.hdf5")


def extract_cdt_timestamp(packet_file):
    m = re.search(r"(20\d\d)_(\d\d)_(\d\d)_(\d\d)_(\d\d)_(\d\d)_CDT", str(packet_file))
    if m is None:
        return None
    yyyy, mm, dd, hh, minute, ss = map(int, m.groups())
    return pd.Timestamp(yyyy, mm, dd, hh, minute, ss)


def collect_runs(roots, flow_root):
    rows = []
    for root in roots:
        if not root.exists():
            raise RuntimeError(f"Frozen root does not exist: {root}")
        for npz in sorted(root.rglob("selected_packets.npz")):
            summary_path = npz.with_name("summary.json")
            if not summary_path.exists():
                print(f"SKIP missing summary.json: {npz}", flush=True)
                continue
            with summary_path.open() as f:
                summary = json.load(f)
            packet_file = summary.get("packet_file")
            if not packet_file:
                print(f"SKIP missing packet_file: {summary_path}", flush=True)
                continue
            flow_path = get_flow_path(packet_file, flow_root)
            rows.append({
                "npz": npz,
                "frozen_source": str(root),
                "summary": summary,
                "packet_file": packet_file,
                "flow_path": flow_path,
                "timestamp": extract_cdt_timestamp(packet_file),
                "ext_io": int(summary.get("ext_io_group", -1)),
                "frozen_unweighted_ms": summary.get("unweighted_tau_ms"),
                "frozen_weighted_ms": summary.get("weighted_tau_ms"),
            })
    # Guard against inadvertently mixing multiple runs for the same raw packet file.
    seen = set()
    for row in rows:
        if row["packet_file"] in seen:
            raise RuntimeError(f"Duplicate frozen source file: {row['packet_file']}")
        seen.add(row["packet_file"])
    return sorted(rows, key=lambda x: (str(x["timestamp"]), x["packet_file"]))


def load_and_group(npz_path, flow_path):
    with np.load(npz_path, allow_pickle=False) as z:
        selected = z["selected"]
        dt = np.asarray(z["selected_dt_us"], dtype=float)
        ext = np.asarray(z["selected_ext_index"], dtype=np.int64)
    if not len(selected) == len(dt) == len(ext):
        raise RuntimeError("Saved selected/DT/EXT arrays are misaligned in length")

    ped = flow_ped_map(flow_path, KE_PER_ADC, chunk=1_000_000)
    electronic_id = flow_key(*(selected[k].astype(np.int64) for k in CHANNEL_KEYS))
    baseline = np.fromiter((ped.get(int(k), np.nan) for k in electronic_id),
                           dtype=float, count=len(selected))
    valid = np.isfinite(baseline) & np.isfinite(dt) & (dt >= 0) & (dt < 200)
    coverage = float(np.mean(valid))
    if coverage < .95:
        raise RuntimeError(f"Unacceptably low matched-packet coverage: {coverage:.3%}")

    s = pd.DataFrame({k: selected[k][valid].astype(np.int64) for k in CHANNEL_KEYS})
    s["ext_index"] = ext[valid]
    s["t_us"] = dt[valid]
    s["q_ADC"] = selected["dataword"][valid].astype(float) - baseline[valid]

    # EXACTLY as the exploratory comparison: one provisional group per EXT/channel.
    # The span <=10 us requirement is applied AFTER grouping, not as a gap split.
    groups = s.groupby(GROUP_KEYS, sort=False).agg(
        n=("q_ADC", "size"), q_sum_ADC=("q_ADC", "sum"),
        t_first=("t_us", "min"), t_last=("t_us", "max")
    ).reset_index()
    groups["span_us"] = groups["t_last"] - groups["t_first"]
    return groups, len(selected), int(valid.sum()), coverage


def attenuation(t, q40, tau):
    return q40 * np.exp(-(np.asarray(t) - 40.0) / tau)


def fit_one_sample(groups, selection, rng, n_boot):
    selected = groups[(groups["span_us"] <= 10) & selection(groups)]
    mpvs = []
    for left, right in TIME_RANGES:
        q = selected.loc[(selected.t_last >= left) & (selected.t_last < right),
                         "q_sum_ADC"].to_numpy()
        if len(q) < 100:
            raise RuntimeError(f"Only {len(q)} groups in [{left}, {right}) us")
        central = fit_langau(q, bin_width=1.0, qmin=0.0, qmax=50.0)
        boot = []
        for _ in range(n_boot):
            qb = q[rng.integers(0, len(q), size=len(q))]
            try:
                val = fit_langau(qb, bin_width=1.0, qmin=0.0, qmax=50.0)["mpv"]
                if np.isfinite(val):
                    boot.append(val)
            except Exception:
                pass
        err = np.std(boot, ddof=1) if len(boot) >= max(2, n_boot // 2) else np.nan
        mpvs.append({
            "t_start": left, "t_end": right, "t_center": (left + right) / 2,
            "N": len(q), "MPV_ADC": float(central["mpv"]),
            "MPV_error": float(err), "chi2_ndf": float(central["red_chi2"]),
            "n_boot_success": len(boot)
        })
    d = pd.DataFrame(mpvs)
    t = d.t_center.to_numpy(dtype=float)
    q = d.MPV_ADC.to_numpy(dtype=float)
    sigma = d.MPV_error.to_numpy(dtype=float)
    if not np.all(np.isfinite(sigma) & (sigma > 0)):
        raise RuntimeError("Some bootstrap MPV uncertainties are missing or zero")
    pars, cov = curve_fit(attenuation, t, q, sigma=sigma, absolute_sigma=True,
                          p0=[q[0], 1300], bounds=([0, 1], [np.inf, np.inf]),
                          maxfev=20000)
    pred = attenuation(t, *pars)
    chi2 = float(np.sum(((q - pred) / sigma) ** 2))
    return d, {"tau_us": float(pars[1]), "tau_err_us": float(np.sqrt(cov[1, 1])),
               "chi2": chi2, "ndf": len(t) - 2, "n_groups": int(len(selected))}


def main():
    a = parse_args()
    roots = a.frozen_root or [Path("frozen_brooke10_bs20"), Path("frozen_hit_timeseries_bs20")]
    print("Frozen source directories:", ", ".join(str(p) for p in roots), flush=True)
    runs = collect_runs(roots, a.flow_root)
    if a.only:
        runs = [r for r in runs if a.only in Path(r["packet_file"]).name]
    if not runs:
        raise SystemExit("No selected_packets.npz + summary.json pairs found. Inspect your frozen output folder.")
    print(f"Found {len(runs)} candidate frozen runs", flush=True)
    for r in runs:
        label = Path(r["packet_file"]).name
        is_present = bool(r["flow_path"] and r["flow_path"].is_file())
        print(f"[{r['frozen_source']}] {label}: EXT IO {r['ext_io']}, FLOW {'OK' if is_present else 'MISSING'}, "
              f"frozen w={r['frozen_weighted_ms']} ms", flush=True)
    eligible = sum(
        bool(r["flow_path"] and r["flow_path"].is_file())
        and (r["ext_io"] == 6 or a.allow_non_io6)
        for r in runs
    )
    print(f"Eligible for processing: {eligible}/{len(runs)} ",
          "(matching FLOW exists; EXT selection accepted)", flush=True)
    if a.list_only:
        return

    a.outdir.mkdir(parents=True, exist_ok=True)
    result_rows = []
    all_mpvs = []
    for index, r in enumerate(runs, 1):
        rawname = Path(r["packet_file"]).name
        if r["ext_io"] != 6 and not a.allow_non_io6:
            print(f"SKIP {rawname}: EXT IO {r['ext_io']} (not IO6)", flush=True)
            continue
        if r["flow_path"] is None or not r["flow_path"].is_file():
            print(f"SKIP {rawname}: matching FLOW file not found under {a.flow_root}", flush=True)
            continue
        print(f"\n[{index}/{len(runs)}] {rawname}", flush=True)
        try:
            groups, npackets, nmatched, coverage = load_and_group(r["npz"], r["flow_path"])
            row = {
                "time_CDT": r["timestamp"], "packet_file": rawname,
                "frozen_source": r["frozen_source"],
                "flow_file": str(r["flow_path"]),
                "ext_io_group": r["ext_io"], "n_packets": npackets,
                "matched_packets": nmatched, "pedestal_coverage": coverage,
                "frozen_unweighted_ms": r["frozen_unweighted_ms"],
                "frozen_weighted_ms": r["frozen_weighted_ms"],
            }
            for name, selection in [
                ("compact_n2", lambda x: x.n == 2),
                ("compact_nge2", lambda x: x.n >= 2),
            ]:
                # Deterministic per-file and per-sample bootstrap.
                run_match = re.search(r"packet-(\d+)", rawname)
                file_seed = int(run_match.group(1)) if run_match else index
                seed = 20261010 + file_seed * 3 + (0 if name == "compact_n2" else 1)
                d, fit = fit_one_sample(groups, selection, np.random.default_rng(seed), a.bootstrap)
                row[f"{name}_tau_ms"] = fit["tau_us"] / 1000
                row[f"{name}_tau_err_ms"] = fit["tau_err_us"] / 1000
                row[f"{name}_chi2_ndf"] = fit["chi2"] / fit["ndf"]
                row[f"{name}_n_groups"] = fit["n_groups"]
                for item in d.to_dict("records"):
                    all_mpvs.append({"packet_file": rawname, "time_CDT": r["timestamp"],
                                     "sample": name, **item})
                print(f"  {name}: tau = {row[f'{name}_tau_ms']:.3f} +/- "
                      f"{row[f'{name}_tau_err_ms']:.3f} ms, "
                      f"chi2/ndf={row[f'{name}_chi2_ndf']:.2f}", flush=True)
            result_rows.append(row)
            # Save partial results after every file, in case a long series is interrupted.
            pd.DataFrame(result_rows).to_csv(a.outdir / "paired_lifetime_timeseries.csv", index=False)
            pd.DataFrame(all_mpvs).to_csv(a.outdir / "new_method_mpvs.csv", index=False)
        except Exception as exc:
            print(f"FAILED {rawname}: {type(exc).__name__}: {exc}", flush=True)

    if not result_rows:
        raise SystemExit("No files completed. Check FLOW paths, EXT IO, and inputs above.")

    result = pd.DataFrame(result_rows).sort_values("time_CDT")
    result.to_csv(a.outdir / "paired_lifetime_timeseries.csv", index=False)
    plt.figure(figsize=(12, 6))
    for key, label, style in [
        ("frozen_unweighted_ms", "Frozen unweighted (28 ticks)", "s-"),
        ("frozen_weighted_ms", "Frozen weighted (28 ticks)", "^-"),
        ("compact_n2_tau_ms", "New compact n=2", "o-"),
        ("compact_nge2_tau_ms", "New compact n>=2", "d-"),
    ]:
        x = result["time_CDT"]
        y = pd.to_numeric(result[key], errors="coerce")
        if key.startswith("compact"):
            err = pd.to_numeric(result[key.replace("_tau_ms", "_tau_err_ms")], errors="coerce")
            plt.errorbar(x, y, yerr=err, fmt=style, capsize=3, label=label)
        else:
            plt.plot(x, y, style, label=label)
    plt.xlabel("Acquisition start time (CDT)")
    plt.ylabel("Effective lifetime [ms]")
    plt.title("DUNE 2x2: paired frozen vs compact-group packet lifetime")
    plt.grid(alpha=.25)
    plt.legend()
    plt.gcf().autofmt_xdate()
    plt.tight_layout()
    plt.savefig(a.outdir / "paired_lifetime_timeseries.png", dpi=200)
    plt.close()
    print(f"\nFinished {len(result)} runs. Saved comparison CSV, MPVs, and plot to {a.outdir}", flush=True)


if __name__ == "__main__":
    main()
