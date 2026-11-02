"""
postprocess_GPR_DELv1.py

Compares FLORIDyn's own predicted DEL (CLC_history_<port>.xlsx) against the REAL,
MLife-measured DEL computed separately on each of the 5 contiguous 600s time windows
DEL1..DEL5 carved out of [1200, 4200]s (FAST.Farm .out cut at 1200 s, TMax 4200).
FLORIDyn's prediction is a time-weighted mean over the same span (whole run after
1200 s, and each 600 s window; CLC_history `time` column, rows before 1200 s ignored).
Besides the per-turbine comparison, the farm sum (T1+T2+T3) per window is compared:
summary_DEL*_farm_sum.csv and figures_DEL*/DEL_farm_sum_<load>_layout_<A|B>.png.
Measured DEL = max across the 6 MLife channel columns (MEASURED_CHANNEL; "last" =
Tower_6 / Blade_6 only). Where postprocess.py can only compare FLORIDyn's
average against one whole-run MLife number, this script shows whether that single
average is actually representative of the window-to-window spread of the real
short-term loads, for every port (all 12 "cases" in this family: the 10 active
lambda-sweep ports plus the 2 no-yaw baselines).

This is the SeaState_v1 family (ports 15101-15112, MLife workbooks
"FLORIDyn_zeroMQ_DEL<chunk>_v1_Short-term.xlsx"). See postprocess_GPR_DELv1.py /
postprocess_GPR_DELv1.py for the original/base / SeaState_v2 families
(ports 15001-15012 / 15201-15212).

Inputs:
    <RESULTS_DIR>/CLC_history_<port>.xlsx
        FLORIDyn's own predicted DEL history (see postprocess.py's docstring) --
        only the T{1,2,3}_DEL_tower/blade columns are used here, averaged per
        turbine over the run minus its first N_TRANSIENT_ROWS (=4) rows,
        which are the start-up transient. Baseline ports have no CLC_history (they never
        run FLORIDyn) and so contribute no predicted value.
    <MLIFE_DIR>/FLORIDyn_zeroMQ_DEL<chunk>_Short-term.xlsx   (chunk = 1..6)
        One MLife 'Short-term DELs' workbook per 600s window (see
        Mlife/CertTest/FLORIDyn_zeroMQ_DEL<chunk>.mlif for the corresponding MLife
        input-file list, which points at splitter_DEL.py's
        <port>DEL<chunk>.T<turbine>.out chunk files). Read the same way as
        postprocess.py's load_measured_DEL, except the 'Filename' column here reads
        like '15101DEL3.T1' (port + DEL<chunk> + turbine) instead of '15101.T1'.

Output (written under postprocess/):
    summary_DEL.csv    one row per (port, turbine, load type): FLORIDyn's predicted
                        post-transient average, each of the 6 measured chunk values, their mean
                        and std, and the predicted-vs-measured-mean %% difference.
    figures_DEL/DEL_<load_type>_layout_<A|B>.png   (4 total: blade/tower x A/B)
                        FLORIDyn's single predicted post-transient average (diamond)
                        plotted alongside the 6 measured per-window values (small dots) and
                        their mean +/- std (errorbar), per turbine per run.

Usage:
    python3 postprocess_GPR_DELv1.py
    python3 postprocess_GPR_DELv1.py --results-dir /path --mlife-dir /path --ports 15101 15102

GPR VARIANT: ports 15113-15122 (tower residual corrector = GPR instead of FFNN; same sims as
15101-15111). MLife workbooks FLORIDyn_zeroMQ_DEL<k>_GPR_v1_Short-term.xlsx. Baselines 15106/15112
have no measured data in these workbooks (not in the GPR .mlif) -> only runs are drawn.
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Large fonts everywhere (user request 2026-10-12)
plt.rcParams.update({
    "font.size": 24,
    "axes.labelsize": 28,
    "axes.titlesize": 26,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 20,
    "figure.titlesize": 28,
})
import numpy as np
import pandas as pd

# Windows: the FAST.Farm .out are cut to t >= T_CUT_S and split into N_CHUNKS
# windows of CHUNK_S seconds (<port>DEL<k>.T<n>.out, DEL_k = [T_CUT_S + CHUNK_S*(k-1),
# T_CUT_S + CHUNK_S*k)); TMax 4200 -> DEL1..DEL5.
N_CHUNKS = 5
T_CUT_S = 1200.0
CHUNK_S = 600.0

# Same lambda-sweep / layout structure as postprocess.py (see that file's docstring
# for what these mean); kept in sync by construction since both families ran the
# same sweep, just under a different sea state.
LAMBDA_WEIGHTS = {
    15113: {"lambda_blade": 0, "lambda_tower": 0, "lambda_P": 1, "layout": "A"},
    15114: {"lambda_blade": 1, "lambda_tower": 0, "lambda_P": 1, "layout": "A"},
    15115: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 1, "layout": "A"},
    15116: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 0, "layout": "A"},
    15117: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 5, "layout": "A"},
    15118: {"lambda_blade": 0, "lambda_tower": 0, "lambda_P": 1, "layout": "B"},
    15119: {"lambda_blade": 1, "lambda_tower": 0, "lambda_P": 1, "layout": "B"},
    15120: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 1, "layout": "B"},
    15121: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 0, "layout": "B"},
    15122: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 5, "layout": "B"},
}
BASELINE_PORTS = {"A": 15106, "B": 15112}
PORT_LAYOUT = {**{p: v["layout"] for p, v in LAMBDA_WEIGHTS.items()},
               **{port: layout for layout, port in BASELINE_PORTS.items()}}
ALL_PORTS = sorted(LAMBDA_WEIGHTS) + sorted(BASELINE_PORTS.values())

# Optional FLORIDyn farm-sum DEL (T1+T2+T3) of the baselines PER WINDOW, e.g.
# {"A": {"tower": [DEL1, ..., DEL5], "blade": [...]}, "B": {...}}. The baselines
# never run FLORIDyn, so without values here their bars show the measured side only.
# FLORIDyn predicted baseline DEL per turbine WITH the GPR corrector (user-supplied, whole-run value;
# baselines never run FLORIDyn). {port: {turbine: {load: value}}}
BASELINE_PRED_DEL = {
    15106: {1: {"blade": 18550.0, "tower": 34517.0}, 2: {"blade": 16068.0, "tower": 33343.0}, 3: {"blade": 16161.0, "tower": 31667.0}},
    15112: {1: {"blade": 18550.0, "tower": 34517.0}, 2: {"blade": 15999.0, "tower": 30468.0}, 3: {"blade": 16050.0, "tower": 29422.0}},
}
# Per-window farm-sum bars for the baselines: the whole-run value repeated for every window (no per-window FLORIDyn value).
BASELINE_PRED_DEL_SUM_CHUNKS = {
    "A": {"tower": [99527.0] * N_CHUNKS, "blade": [50779.0] * N_CHUNKS},
    "B": {"tower": [94407.0] * N_CHUNKS, "blade": [50599.0] * N_CHUNKS},
}

DEFAULT_RESULTS_DIR = (
    "/mnt/d/zero_model/FLORIDyn_Matlab-main_3_floating/FLORIDyn_Matlab-main_3/"
    "Simulations/2026_3T_DEL/3T_IEA15_floating_00/Results"
)
DEFAULT_MLIFE_DIR = "/mnt/d/zero_model/fastfarm_simulation/Mlife/CertTest"
MLIFE_ROOT = "FLORIDyn_zeroMQ"
MLIFE_FAMILY_SUFFIX = "_GPR_v1"   # "" here (base family); "_v1"/"_v2" in the other two scripts.
                           # Actual workbook names are "<root>_DEL<chunk><suffix>_Short-term.xlsx"
                           # -- the suffix comes AFTER "_DEL<chunk>", not before.

N_TURBINES = 3
LOAD_TYPES = ("tower", "blade")

# How the 6 MLife channel columns of a row become one measured DEL: "max" (max
# across Tower_1..6 / Blade_1..6, matches FLORIDyn's magnitude on the 15xxx
# workbooks) or "last" (Tower_6 / Blade_6 only).
MEASURED_CHANNEL = "max"


def reduce_channels(values: list):
    """One measured DEL from a row's channel values (see MEASURED_CHANNEL)."""
    if MEASURED_CHANNEL == "last":
        return float(values[-1]) if pd.notna(values[-1]) else None
    vals = [v for v in values if pd.notna(v)]
    return float(max(vals)) if vals else None


def window_bounds(chunk: int) -> tuple:
    """[t0, t1] FAST.Farm time of DEL window `chunk` (1-based)."""
    t0 = T_CUT_S + CHUNK_S * (chunk - 1)
    return t0, t0 + CHUNK_S


def predicted_mean(clc: pd.DataFrame, col: str, t0: float, t1: float,
                   label: str = "") -> float:
    """Time-weighted mean of a CLC_history column over [t0, t1] FAST.Farm time
    (`time` column = time.current + 60, see FLORIDyn_step_*.m).

    Each row is a re-optimization: its predicted DEL holds from its `time` until the
    next row (zero-order hold; the last row holds until t1). Weighting by the time
    each row is in force, instead of a plain mean of the rows inside [t0, t1], keeps
    the result meaningful when a 600 s window has only one or two re-optimizations
    (or none), and it counts the row already in force at t0.
    Falls back to a plain mean of all rows if the file has no `time` column."""
    if col not in clc.columns or clc.empty:
        return np.nan
    if "time" not in clc.columns:
        print(f"  [warn] {label or 'CLC_history'}: no 'time' column -- plain mean of all rows")
        return float(clc[col].mean())
    d = clc[["time", col]].dropna().sort_values("time")
    start = d["time"].to_numpy(dtype=float)
    end = np.append(start[1:], np.inf)
    w = np.clip(np.minimum(end, t1) - np.maximum(start, t0), 0.0, None)
    if w.sum() <= 0:
        return np.nan
    return float(np.sum(w * d[col].to_numpy(dtype=float)) / w.sum())


def load_predicted_DEL(results_dir: Path, port: int, n_chunks: int = N_CHUNKS) -> dict:
    """FLORIDyn's own predicted DEL for one port, per turbine:
    {turbine: {"tower": {"all": mean, "chunks": [w_1..w_n]}, "blade": {...}}}
    "all" = time-weighted mean over the whole run after T_CUT_S (= union of the
    windows), "chunks" = the same over each 600 s window (see predicted_mean).
    Empty dict if this port has no CLC_history (baseline ports never run FLORIDyn)."""
    path = results_dir / f"CLC_history_{port}.xlsx"
    if not path.exists():
        return {}
    clc = pd.read_excel(path)
    if clc.empty:
        return {}
    t_end = window_bounds(n_chunks)[1]
    out = {}
    for t in range(1, N_TURBINES + 1):
        vals = {}
        for load_type in LOAD_TYPES:
            col = f"T{t}_DEL_{load_type}"
            if col in clc.columns:
                vals[load_type] = {
                    "all": predicted_mean(clc, col, T_CUT_S, t_end, label=path.name),
                    "chunks": [predicted_mean(clc, col, *window_bounds(k), label=path.name)
                               for k in range(1, n_chunks + 1)],
                }
        # GPR campaign: load-model uncertainty of the tower DEL logged by FLORIDyn_step_<port>.m
        # (sigma of the FFNN predictor, of the GPR corrector, and sqrt(s1^2+s2^2)); time-weighted
        # mean over the same span as "all". Blade has no sigma (FFNN predictor only, not logged).
        if "tower" in vals:
            for kind in ("FFNN", "GPR", "comb"):
                c = f"T{t}_sigma_{kind}_tower"
                if c in clc.columns:
                    vals["tower"][f"sigma_{kind}"] = predicted_mean(clc, c, T_CUT_S, t_end, label=path.name)
        if vals:
            out[t] = vals
    return out


def load_measured_DEL_chunks(mlife_dir: Path, mlife_root: str, family_suffix: str = "",
                              n_chunks: int = N_CHUNKS) -> dict:
    """Real per-window DEL from the n_chunks MLife 'Short-term DELs' workbooks
    (<mlife_root>_DEL<chunk><family_suffix>_Short-term.xlsx -- family_suffix, e.g.
    "_v1", comes AFTER "_DEL<chunk>", not before). Mirrors postprocess.py's
    load_measured_DEL block-parsing, except 'Filename' here is '<port>DEL<chunk>.T
    <turbine>' -- the DEL<chunk> infix is parsed out (and cross-checked against the
    chunk number of the workbook it came from) rather than kept as part of the port.

    Returns {port: {turbine: {"tower": [v_1..v_n_chunks], "blade": [...]}}}, with
    a chunk's value left as None if that chunk's workbook is missing or doesn't
    contain that port/turbine.
    """
    fname_re = re.compile(r"^(\d+)DEL(\d+)\.T(\d+)$")
    result: dict = {}

    def _parse_block(raw: pd.DataFrame, block_title: str, channel_prefix: str,
                      key: str, chunk: int):
        title_rows = raw.index[raw[0].astype(str).str.startswith(block_title, na=False)]
        if len(title_rows) == 0:
            raise ValueError(f"no '{block_title}' block")
        header_row = title_rows[0] + 1
        headers = raw.iloc[header_row].tolist()
        filename_col = headers.index("Filename")
        channel_cols = [i for i, h in enumerate(headers)
                        if isinstance(h, str) and h.startswith(channel_prefix)]

        row = header_row + 1
        while row < len(raw) and not (isinstance(raw.iat[row, filename_col], str)
                                       and fname_re.match(raw.iat[row, filename_col])):
            row += 1  # skip the L_Ult / L_MF summary rows
        while row < len(raw):
            fname = raw.iat[row, filename_col]
            if not (isinstance(fname, str) and fname_re.match(fname)):
                break
            port, fchunk, turbine = (int(g) for g in fname_re.match(fname).groups())
            if fchunk != chunk:
                print(f"  [warn] {fname}: chunk in filename ({fchunk}) != "
                      f"workbook chunk ({chunk}), using the filename's own value slot")
            value = reduce_channels([raw.iat[row, c] for c in channel_cols])
            if value is not None and 1 <= fchunk <= n_chunks:
                slot = result.setdefault(port, {}).setdefault(turbine, {}) \
                             .setdefault(key, [None] * n_chunks)
                slot[fchunk - 1] = value
            row += 1

    for chunk in range(1, n_chunks + 1):
        xlsx_path = mlife_dir / f"{mlife_root}_DEL{chunk}{family_suffix}_Short-term.xlsx"
        if not xlsx_path.exists():
            print(f"  [missing] {xlsx_path.name}")
            continue
        raw = pd.read_excel(xlsx_path, sheet_name="Short-term DELs", header=None)
        try:
            _parse_block(raw, "Tower_Base", "Tower_", "tower", chunk)
        except ValueError as e:
            print(f"  [warn] {xlsx_path.name}: {e}")
        try:
            _parse_block(raw, "Blade_Root", "Blade_", "blade", chunk)
        except ValueError as e:
            print(f"  [warn] {xlsx_path.name}: {e}")

    return result


def plot_DEL_chunks_vs_predicted(predicted: dict, measured: dict, run_number: dict,
                                  layout: str, load_type: str, fig_dir: Path):
    """One figure per (layout, load_type): FLORIDyn's single predicted average
    (diamond) plotted against the real per-600s-window measured DEL (small dots for
    each of the up-to-6 chunks, plus a mean +/- std errorbar), one x slot per run,
    turbines offset side by side within a run and color-coded. Ports with no
    predicted value (the baseline) are appended after the runs, labeled 'BL', with
    only the measured side shown."""
    ports = [p for p in ALL_PORTS if PORT_LAYOUT.get(p) == layout and p in measured]
    if not ports:
        return
    ports = sorted(ports, key=lambda p: (p not in run_number, run_number.get(p, p)))

    turbine_colors = {1: "C0", 2: "C1", 3: "C2"}
    turbine_dx = {1: -0.22, 2: 0.0, 3: 0.22}
    fig, ax = plt.subplots(figsize=(max(15, 1.3 * len(ports) + 6), 7))

    xticklabels = []
    legend_done = set()
    for i, port in enumerate(ports):
        xticklabels.append(str(run_number[port]) if port in run_number else "BL")
        pred_here = predicted.get(port, {})
        for t in range(1, N_TURBINES + 1):
            chunks = measured.get(port, {}).get(t, {}).get(load_type)
            if not chunks:
                continue
            x = i + turbine_dx[t]
            vals = np.array([v for v in chunks if v is not None], dtype=float)
            if vals.size:
                mean, std = vals.mean(), vals.std()
                ax.scatter([x] * vals.size, vals, color=turbine_colors[t], alpha=0.35,
                           s=22, zorder=2, linewidths=0)
                ax.errorbar([x], [mean], yerr=[std], fmt="o", color=turbine_colors[t],
                            ecolor=turbine_colors[t], elinewidth=1.5, capsize=4,
                            markersize=7, zorder=3,
                            label=f"T{t} measured (mean\u00b1std)"
                            if ("m", t) not in legend_done else None)
                legend_done.add(("m", t))

            pv = pred_here.get(t, {}).get(load_type, {}).get("all")
            if pv is not None and np.isfinite(pv):
                ax.scatter([x], [pv], color=turbine_colors[t], marker="D", s=90,
                           edgecolors="k", linewidths=0.8, zorder=4,
                           label=f"T{t} FLORIDyn"
                           if ("p", t) not in legend_done else None)
                legend_done.add(("p", t))

    ax.set_xticks(range(len(ports)))
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel("Run")
    ax.set_ylabel(f"{load_type.capitalize()} DEL [kNm]")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=20, loc="upper left", bbox_to_anchor=(1.01, 1.0))

    fig.tight_layout()
    fig.savefig(fig_dir / f"DEL_{load_type}_layout_{layout}.png", dpi=150)
    plt.close(fig)


def write_summary(predicted_by_port: dict, measured: dict, run_number: dict,
                   out_csv: Path):
    rows = []
    for port in ALL_PORTS:
        layout = PORT_LAYOUT.get(port)
        for t in range(1, N_TURBINES + 1):
            for load_type in LOAD_TYPES:
                chunks = measured.get(port, {}).get(t, {}).get(load_type)
                pred = predicted_by_port.get(port, {}).get(t, {}).get(load_type, {})
                pv = pred.get("all")
                if chunks is None and pv is None:
                    continue
                chunks = chunks or [None] * N_CHUNKS
                vals = np.array([v for v in chunks if v is not None], dtype=float)
                mean = float(vals.mean()) if vals.size else np.nan
                std = float(vals.std()) if vals.size else np.nan
                pct_diff = (mean - pv) / pv * 100 if (pv not in (None, 0) and not np.isnan(mean)) else np.nan
                row = {
                    "port": port, "run": run_number.get(port, ""), "layout": layout,
                    "turbine": t, "load_type": load_type,
                    "predicted_mean": pv, "measured_mean": mean, "measured_std": std,
                    "pct_diff_measured_vs_predicted": pct_diff,
                }
                for c, v in enumerate(chunks, start=1):
                    row[f"measured_chunk{c}"] = v
                for c, v in enumerate(pred.get("chunks", [None] * N_CHUNKS), start=1):
                    row[f"predicted_chunk{c}"] = v
                rows.append(row)

    pd.DataFrame(rows).to_csv(out_csv, index=False)


def farm_sums(predicted: dict, measured: dict, port: int, load_type: str,
              n_chunks: int = N_CHUNKS) -> tuple:
    """Per window, the farm sum (T1+T2+T3) of the MLife DEL and of FLORIDyn's
    time-weighted mean over that window: (meas[n_chunks], pred[n_chunks]).
    NaN wherever one of the 3 turbines is missing."""
    meas = np.full(n_chunks, np.nan)
    pred = np.full(n_chunks, np.nan)
    for k in range(n_chunks):
        m = [(measured.get(port, {}).get(t, {}).get(load_type) or [None] * n_chunks)[k]
             for t in range(1, N_TURBINES + 1)]
        if all(v is not None for v in m):
            meas[k] = sum(m)
        p = [(predicted.get(port, {}).get(t, {}).get(load_type, {}).get("chunks")
              or [np.nan] * n_chunks)[k] for t in range(1, N_TURBINES + 1)]
        p = np.array([np.nan if v is None else v for v in p], dtype=float)
        if np.all(np.isfinite(p)):
            pred[k] = p.sum()
    return meas, pred


def write_farm_sum_summary(predicted: dict, measured: dict, run_number: dict,
                           ports: list, out_csv: Path, n_chunks: int = N_CHUNKS):
    """One row per (port, load type, window): farm sum of the measured DEL of that
    600 s window vs. farm sum of FLORIDyn's means over the same window."""
    rows = []
    for port in ports:
        for load_type in LOAD_TYPES:
            meas, pred = farm_sums(predicted, measured, port, load_type, n_chunks)
            if np.all(np.isnan(meas)) and np.all(np.isnan(pred)):
                continue
            for k in range(n_chunks):
                t0, t1 = window_bounds(k + 1)
                rows.append({
                    "port": port, "run": run_number.get(port, "BL"),
                    "layout": PORT_LAYOUT.get(port), "load_type": load_type,
                    "window": f"DEL{k + 1}", "t_start": t0, "t_end": t1,
                    "sum_measured": meas[k], "sum_predicted": pred[k],
                    "pct_pred_vs_meas": (pred[k] - meas[k]) / meas[k] * 100,
                })
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def plot_farm_sum_windows(predicted: dict, measured: dict, run_number: dict,
                          layout: str, load_type: str, fig_dir: Path,
                          n_chunks: int = N_CHUNKS):
    """One subplot per run (+ baseline): farm sum of the DEL in each 600 s window,
    measured (MLife) vs. FLORIDyn's time-weighted mean over that window."""
    ports = [p for p in ALL_PORTS if PORT_LAYOUT.get(p) == layout
             and (p in measured or p in predicted)]
    if not ports:
        return
    ports = sorted(ports, key=lambda p: (p not in run_number, run_number.get(p, p)))

    ncols = 3
    nrows = -(-len(ports) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(17, 4.6 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    x = np.arange(1, n_chunks + 1)
    for i, port in enumerate(ports):
        ax = axes[i // ncols][i % ncols]
        meas, pred = farm_sums(predicted, measured, port, load_type, n_chunks)
        ax.plot(x, meas, "o-", color="C1", label="FAST.Farm/MLife")
        err_txt = ""
        if np.any(np.isfinite(pred)):
            ax.plot(x, pred, "D--", color="C0", label="FLORIDyn")
            if np.all(np.isfinite(meas)) and np.all(np.isfinite(pred)):
                err = np.mean((pred - meas) / meas * 100)
                err_txt = f" (mean err {err:+.1f}%)"
        ax.set_title((f"Run {run_number[port]}" if port in run_number
                      else f"BL-{layout}") + err_txt, fontsize=22)
        ax.grid(alpha=0.3)
    for i in range(len(ports), nrows * ncols):
        axes[i // ncols][i % ncols].axis("off")
    for c in range(ncols):
        axes[-1][c].set_xticks(x)
        axes[-1][c].set_xticklabels([f"DEL{k}" for k in x])

    fig.supylabel("Farm DEL sum [kNm]")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, fontsize=22, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_dir / f"DEL_farm_sum_{load_type}_layout_{layout}.png", dpi=150)
    plt.close(fig)


def plot_farm_sum_bars_per_window(predicted: dict, measured: dict, run_number: dict,
                                  load_type: str, fig_dir: Path,
                                  n_chunks: int = N_CHUNKS):
    """Same bar chart as postprocess.py's DEL_farm_sum_<load>.png, one figure per
    600 s window: farm sum (T1+T2+T3) measured vs FLORIDyn for every run, ordered
    1-5, BL-A, 6-10, BL-B, with the (pred-meas)/meas %% panel below. Same y-axis
    across the windows so they can be compared side by side. Baselines get a
    FLORIDyn bar only if given in BASELINE_PRED_DEL_SUM_CHUNKS."""
    columns = []   # (label, port)
    for layout, bl_port in BASELINE_PORTS.items():
        columns += [(str(run_number[p]), p) for p in sorted(LAMBDA_WEIGHTS, key=run_number.get)
                    if PORT_LAYOUT[p] == layout]
        columns.append((f"BL-{layout}", bl_port))
    labels = [c[0] for c in columns]

    meas = np.full((n_chunks, len(columns)), np.nan)
    pred = np.full((n_chunks, len(columns)), np.nan)
    for c, (label, port) in enumerate(columns):
        meas[:, c], pred[:, c] = farm_sums(predicted, measured, port, load_type, n_chunks)
        if label.startswith("BL-"):
            bl = BASELINE_PRED_DEL_SUM_CHUNKS.get(label[3:], {}).get(load_type)
            if bl is not None:
                pred[:, c] = bl[:n_chunks]
    finite = np.concatenate([meas[np.isfinite(meas)], pred[np.isfinite(pred)]])
    if finite.size == 0:
        return
    lo, hi = finite.min(), finite.max()
    pad = 0.15 * (hi - lo) if hi > lo else 0.05 * hi
    pct_all = (pred - meas) / meas * 100
    pct_lim = np.nanmax(np.abs(pct_all)) * 1.15 if np.any(np.isfinite(pct_all)) else 1.0

    x = np.arange(len(labels))
    for k in range(n_chunks):
        t0, t1 = window_bounds(k + 1)
        fig, (ax, ax_pct) = plt.subplots(2, 1, figsize=(max(11, 1.1 * len(labels)), 9),
                                         sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        ax.bar(x - 0.2, meas[k], width=0.4, color="C1", label="FAST.Farm/MLife (sum T1-T3)")
        ax.bar(x + 0.2, pred[k], width=0.4, color="C0",
               label="FLORIDyn predicted (sum of T1-T3 window means)")
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel(f"Farm {load_type} DEL sum [kNm]")
        ax.set_title(f"Farm {load_type} DEL, DEL{k + 1} ({t0:.0f}-{t1:.0f} s): "
                     f"measured vs FLORIDyn")
        ax.legend(fontsize=20)
        ax.grid(axis="y", alpha=0.3)

        pct = np.nan_to_num(pct_all[k])
        ax_pct.bar(x, pct, color=["C2" if p >= 0 else "C3" for p in pct])
        ax_pct.axhline(0, color="k", linewidth=0.8)
        ax_pct.set_ylim(-pct_lim, pct_lim)
        ax_pct.set_ylabel("Error [%]")
        ax_pct.set_xticks(x)
        ax_pct.set_xticklabels(labels)
        ax_pct.set_xlabel("Run")
        ax_pct.grid(axis="y", alpha=0.3)

        fig.tight_layout()
        fig.savefig(fig_dir / f"DEL_farm_sum_{load_type}_DEL{k + 1}.png", dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                         help="Folder containing CLC_history_<port>.xlsx")
    parser.add_argument("--mlife-dir", default=DEFAULT_MLIFE_DIR,
                         help="Folder containing <root>_DEL<chunk>_Short-term.xlsx")
    parser.add_argument("--ports", nargs="+", type=int, default=ALL_PORTS,
                         help="Ports to include (default: all 12 in this family)")
    parser.add_argument("--n-chunks", type=int, default=N_CHUNKS,
                         help=f"Number of 600s chunks (default: {N_CHUNKS})")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        raise SystemExit(
            f"--results-dir does not exist or is not a directory: {results_dir}\n"
            f"(hint: WSL-style '/mnt/d/...' paths only resolve under WSL -- pass a path "
            f"this Python interpreter can actually see, e.g. 'D:\\...' on native Windows)"
        )
    mlife_dir = Path(args.mlife_dir)
    if not mlife_dir.is_dir():
        raise SystemExit(
            f"--mlife-dir does not exist or is not a directory: {mlife_dir}\n"
            f"(hint: WSL-style '/mnt/d/...' paths only resolve under WSL -- pass a path "
            f"this Python interpreter can actually see, e.g. 'D:\\...' on native Windows)"
        )

    out_dir = Path(__file__).parent
    fig_dir = out_dir / "figures_GPR_DELv1"
    fig_dir.mkdir(exist_ok=True)

    run_number = {port: i + 1 for i, port in enumerate(sorted(LAMBDA_WEIGHTS))}

    print(f"Reading FLORIDyn predictions from: {results_dir}")
    print(f"  (time-weighted means, CLC_history before {T_CUT_S:.0f} s discarded)")
    predicted = {port: load_predicted_DEL(results_dir, port, args.n_chunks)
                 for port in args.ports}

    print(f"\nReading MLife per-window measurements from: {mlife_dir}")
    measured = load_measured_DEL_chunks(mlife_dir, MLIFE_ROOT, MLIFE_FAMILY_SUFFIX,
                                         n_chunks=args.n_chunks)
    measured = {p: v for p, v in measured.items() if p in args.ports}

    # Baselines: measured per-window DEL from the original (FFNN-campaign) workbooks (same FAST.Farm runs);
    # predicted = user-supplied GPR-corrected FLORIDyn values (BASELINE_PRED_DEL).
    base_meas = load_measured_DEL_chunks(mlife_dir, MLIFE_ROOT, MLIFE_FAMILY_SUFFIX.replace("_GPR", ""),
                                          n_chunks=args.n_chunks)
    args.ports = list(args.ports) + [p for p in BASELINE_PORTS.values() if p not in args.ports]
    for bp in BASELINE_PORTS.values():
        if bp in base_meas:
            measured[bp] = base_meas[bp]
        predicted[bp] = {t: {lt: {"all": v} for lt, v in d.items()} for t, d in BASELINE_PRED_DEL[bp].items()}

    summary_csv = out_dir / "summary_GPR_DELv1.csv"
    write_summary(predicted, measured, run_number, summary_csv)
    print(f"\nWrote summary: {summary_csv}")

    farm_csv = summary_csv.with_name(summary_csv.stem + "_farm_sum.csv")
    write_farm_sum_summary(predicted, measured, run_number, args.ports, farm_csv,
                           args.n_chunks)
    print(f"Wrote farm-sum summary: {farm_csv}")

    for layout in ("A", "B"):
        for load_type in LOAD_TYPES:
            plot_DEL_chunks_vs_predicted(predicted, measured, run_number,
                                          layout, load_type, fig_dir)
            plot_farm_sum_windows(predicted, measured, run_number,
                                  layout, load_type, fig_dir, args.n_chunks)
    for load_type in LOAD_TYPES:
        plot_farm_sum_bars_per_window(predicted, measured, run_number,
                                      load_type, fig_dir, args.n_chunks)

    print(f"Wrote figures to: {fig_dir}")


if __name__ == "__main__":
    main()
