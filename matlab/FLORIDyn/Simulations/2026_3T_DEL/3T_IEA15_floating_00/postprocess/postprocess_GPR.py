"""
postprocess_GPR.py

Post-processes the FLORIDyn closed-loop control (CLC) results for the
3T_IEA15_floating_00 lambda-weight sensitivity sweep (ports 15013-15017).

Each port ran the same 3-turbine floating layout with a different cost-function
weighting (lambda_blade, lambda_tower, lambda_P) between blade DEL, tower DEL,
and power in the yaw optimizer. For each port this script reads the three
result files FLORIDyn writes during a run:

    CLC_history_<port>.xlsx      one row per applied yaw re-optimization event
                                  (T{1,2,3}_yaw_deg, T{1,2,3}_DEL_tower, T{1,2,3}_DEL_blade)
    enkf_log_<port>.xlsx         EnKF wind state estimate, one row per EnKF cycle
                                  (time_s, ws_hat, ws_std, wd_hat, wd_std, wd_drift_deg,
                                   wd_stable, pred_checked, TI_use_T{1,2,3})
    reliability_log_<port>.xlsx  one row per reliability check (predicted vs the
                                  FAST.Farm-measured outcome of an applied yaw)
                                  (t_applied, t_check, controller_rate, gam_T{1,2,3},
                                   ws_pred, wd_pred, ws_check, wd_check, P_pred_W,
                                   P_bl_pred_W, dP_pred_pct, P_FASTfarm_W, P_bl_check_W,
                                   dP_check_pct, rel_power, err_dP_pct,
                                   rel_DEL_tower_T{1,2,3}, rel_DEL_blade_T{1,2,3})

Also compares against a no-yaw-steering baseline -- but there are two layouts (A:
15013-15017, B: 15018-15022, mirrored), and each layout has its own baseline: port
15006 (see controller_floris/FLORIDyn_15006.py) for layout A, port 15012 (same
controller, mirrored layout) for layout B. Neither baseline ever runs FLORIDyn/the
EnKF, so there's no CLC_history/enkf_log/reliability_log for either -- only their raw
FAST.Farm .out files exist, which get read directly (see load_baseline_power) for a
power-only comparison, each port compared against its own layout's baseline.

Output (written under postprocess/):
    summary_by_port_GPR.csv   one row per port: lambda weights + key aggregate stats,
                           plus power_ratio_vs_baseline (computed against that row's
                           own layout's baseline) if that baseline loaded
    figures_GPR/*.png         comparison plots across all 10 ports, incl. power_vs_baseline_layout_<A|B>.png
                           and, split one figure per layout per load type (4 total:
                           DEL_blade_layout_A.png, DEL_tower_layout_A.png,
                           DEL_blade_layout_B.png, DEL_tower_layout_B.png), each with
                           FLORIDyn-predicted DEL on the left and FAST.Farm/MLife
                           -measured DEL on the right (see plot_DEL_layout_comparison),
                           and DEL_farm_sum_tower.png / DEL_farm_sum_blade.png: farm sum
                           (T1+T2+T3) of the MLife DEL vs the sum of FLORIDyn's
                           time-weighted means, both over the whole run after T_CUT_S
                           (see plot_DEL_farm_sum; numbers also in summary_by_port_GPR.csv:
                           mean_DEL_*_sum, meas_DEL_*_sum, pct_DEL_*_sum_pred_vs_meas)

Usage:
    python3 postprocess_GPR.py
    python3 postprocess_GPR.py --results-dir /path/to/Results --ports 15013 15014 15016
    python3 postprocess_GPR.py --baseline-port-a 15006 --baseline-port-b 15012
                                                # (defaults shown; pass 0 to skip either)

GPR VARIANT: ports 15013-15022 = same sims as 15001-15011 (06/12 baselines reused) but the
tower residual corrector is the GPR instead of the FFNN. Baselines stay 15006 (A) / 15012 (B).
Extra summary columns mean_sigma_{FFNN,GPR,comb}_tower_farm: load-model uncertainty of the tower
DEL, time-weighted mean per turbine, quadrature sum over T1..T3. MLife workbook expected:
FLORIDyn_zeroMQ_GPR_Short-term.xlsx (must include the baseline ports 15006/15012
for the measured baseline lines).
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write PNGs without needing a display
import matplotlib.pyplot as plt

# Large fonts everywhere (user request 2026-10-12)
plt.rcParams.update({
    "font.size": 24,
    "axes.labelsize": 28,
    "axes.titlesize": 26,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 20,
    "legend.framealpha": 1.0,  # opaque legend box
    "figure.titlesize": 28,
})
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

# Cost-function weights per port, read directly from each port's
# CLC_init_floating_fun_ZMQ_<port>.m (state.lambda_blade / lambda_tower / lambda_P).
# 15018-15022 are the same 5 lambda combos re-run on the mirrored turbine layout
# (see the .fstf WT_X/WT_Y comparison against 15013-15017) -- "layout" flags which.
LAMBDA_WEIGHTS = {
    15013: {"lambda_blade": 0, "lambda_tower": 0, "lambda_P": 1, "layout": "A"},
    15014: {"lambda_blade": 1, "lambda_tower": 0, "lambda_P": 1, "layout": "A"},
    15015: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 1, "layout": "A"},
    15016: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 0, "layout": "A"},
    15017: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 5, "layout": "A"},
    15018: {"lambda_blade": 0, "lambda_tower": 0, "lambda_P": 1, "layout": "B"},
    15019: {"lambda_blade": 1, "lambda_tower": 0, "lambda_P": 1, "layout": "B"},
    15020: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 1, "layout": "B"},
    15021: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 0, "layout": "B"},
    15022: {"lambda_blade": 1, "lambda_tower": 1, "lambda_P": 5, "layout": "B"},
}

DEFAULT_RESULTS_DIR = (
    "/mnt/d/zero_model/FLORIDyn_Matlab-main_3_floating/FLORIDyn_Matlab-main_3/"
    "Simulations/2026_3T_DEL/3T_IEA15_floating_00/Results"
)

DEFAULT_FASTFARM_DIR = "/mnt/d/zero_model/fastfarm_simulation"

DEFAULT_MLIFE_XLSX = (
    "/mnt/d/zero_model/fastfarm_simulation/Mlife/CertTest/FLORIDyn_zeroMQ_GPR_Short-term.xlsx"
)

# Baseline (no yaw steering, greedy control) runs for comparison -- one per layout,
# see controller_floris/FLORIDyn_15006.py (layout A) and FLORIDyn_15012.py (layout B,
# same controller, mirrored layout). Neither ever runs FLORIDyn/the EnKF, so neither
# has a CLC_history/enkf_log/reliability_log; only their raw FAST.Farm .out files exist.
BASELINE_PORTS = {"A": 15006, "B": 15012}
N_TURBINES = 3

# Bar color per layout in the power and yaw plots (A blue, B purple).
LAYOUT_COLORS = {"A": "C0", "B": "C4"}

# FLORIDyn's predicted farm-sum DEL (T1+T2+T3) for the baselines (yaw 0), computed
# separately (the baseline ports never run FLORIDyn, so there's no CLC_history).
# Shown next to the measured baseline bars in DEL_farm_sum_<load>.png.
BASELINE_PRED_DEL_SUM = {   # FLORIDyn predicted baseline farm sum WITH the GPR corrector (user-supplied per-turbine values)
    "A": {"tower": 143278.0, "blade": 50907.0},
    "B": {"tower": 142926.0, "blade": 51060.0},
}

# Every known port (lambda-sweep or baseline) mapped to its layout -- used to split
# plots by layout without repeating the LAMBDA_WEIGHTS/BASELINE_PORTS lookup logic.
PORT_LAYOUT = {**{p: v["layout"] for p, v in LAMBDA_WEIGHTS.items()},
               **{port: layout for layout, port in BASELINE_PORTS.items()}}

# Signal-processing constants copied from the live controller_floris/FLORIDyn_<port>.py
# bridge (_average_FF_signals), so the "real" wind reconstructed here matches what the
# EnKF actually received during the run.
FF_FS = 10.0          # assumed sample rate of the .out file (10 Hz, i.e. dt=0.1s)
FF_CUTOFF = 0.02       # low-pass cutoff used in the live pipeline
FF_ORDER = 2
FF_WINDOW_S = 50.0     # trailing window (N_eff=500 samples at 10 Hz) averaged/RMS'd
                       # at each query time, matching the live "last N samples" logic

# FLORIDyn's predicted DEL is averaged over the same span MLife sees: the FAST.Farm
# .out files are cut to t >= T_CUT_S before MLife (start-up transient), so the
# CLC_history rows are cut the same way, by their `time` column (FAST.Farm time,
# = time.current + 60, see FLORIDyn_step_*.m). T_END_S = TMax of the campaign.
T_CUT_S = 1200.0
T_END_S = 4200.0


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


def predicted_mean(clc: pd.DataFrame, col: str, t0: float = T_CUT_S,
                   t1: float = T_END_S, label: str = "") -> float:
    """Time-weighted mean of a CLC_history column over [t0, t1] FAST.Farm time.

    Each row is a re-optimization: its predicted DEL holds from its `time` until the
    next row (zero-order hold; the last row holds until t1). Weighting by the time
    each row is in force, instead of a plain mean of the rows inside [t0, t1], keeps
    the result meaningful when re-optimizations are irregular or a window has only
    one or two of them, and it counts the row already in force at t0.
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


# Older runs showed FLORIDyn's own clock (enkf_log/reliability_log) drifting away from
# FAST.Farm's raw clock (the .out files) by ~1000s over a run -- traced to
# controller_floris/FLORIDyn_<port>.py calling into FLORIDyn on essentially every ROSCO
# tick instead of every update_rate=12s: its `t_update_rate` was stored but never
# actually used to gate the call, so the two clocks weren't kept in step. That's fixed
# now (the gating condition uses `self.t_update_rate`), so both clocks should agree and
# nothing needs to be cut or shifted here anymore -- everything below reads straight
# through the full [0, 4600]s of both FLORIDyn's own logs and the raw .out files.


def load_port_data(results_dir: Path, port: int):
    """Read the three result files for one port. Missing files are skipped (None)."""
    files = {
        "clc_history": results_dir / f"CLC_history_{port}.xlsx",
        "enkf_log": results_dir / f"enkf_log_{port}.xlsx",
        "reliability_log": results_dir / f"reliability_log_{port}.xlsx",
    }
    data = {}
    for key, path in files.items():
        if path.exists():
            data[key] = pd.read_excel(path)
        else:
            print(f"  [{port}] missing {path.name}, skipping")
            data[key] = None

    return data


def summarize_port(port: int, data: dict) -> dict:
    """Build one summary row for a port from whatever files were available."""
    row = {"port": port, **LAMBDA_WEIGHTS.get(port, {})}

    clc = data.get("clc_history")
    if clc is not None and not clc.empty:
        row["n_reopt_events"] = len(clc)
        row["mean_yaw_T1_deg"] = clc["T1_yaw_deg"].mean()
        row["mean_yaw_T2_deg"] = clc["T2_yaw_deg"].mean()
        # DELs are averaged over [T_CUT_S, T_END_S] only (same span as MLife).
        if "time" in clc.columns:
            row["n_reopt_events_settled"] = int((clc["time"] >= T_CUT_S).sum())
        for load_type in ("tower", "blade"):
            row[f"mean_DEL_{load_type}_sum"] = sum(
                predicted_mean(clc, f"T{t}_DEL_{load_type}", label=f"CLC_history_{port}")
                for t in range(1, N_TURBINES + 1)
            )
        # GPR campaign: load-model uncertainty of the tower DEL (log-only columns written by
        # FLORIDyn_step_<port>.m). Time-weighted mean per turbine, farm value = quadrature sum
        # over T1..T3 (independence assumed, same assumption as sigma_comb).
        for kind in ("FFNN", "GPR", "comb"):
            cols = [f"T{t}_sigma_{kind}_tower" for t in range(1, N_TURBINES + 1)]
            if all(c in clc.columns for c in cols):
                m = [predicted_mean(clc, c, label=f"CLC_history_{port}") for c in cols]
                row[f"mean_sigma_{kind}_tower_farm"] = float(np.sqrt(np.nansum(np.square(m))))

    rel = data.get("reliability_log")
    if rel is not None and not rel.empty:
        row["n_reliability_checks"] = len(rel)
        row["mean_P_FASTfarm_W"] = rel["P_FASTfarm_W"].mean()
        row["mean_err_dP_pct"] = rel["err_dP_pct"].mean()
        row["mean_abs_err_dP_pct"] = rel["err_dP_pct"].abs().mean()
        row["mean_rel_power"] = rel["rel_power"].mean()

    enkf = data.get("enkf_log")
    if enkf is not None and not enkf.empty:
        row["mean_ws_hat"] = enkf["ws_hat"].mean()
        row["mean_wd_hat"] = enkf["wd_hat"].mean()
        row["frac_wd_stable"] = enkf["wd_stable"].mean()

    return row


def plot_yaw(summary_df: pd.DataFrame, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.bar(summary_df["run"].astype(str), summary_df["mean_yaw_T1_deg"],
           color=[LAYOUT_COLORS.get(l, "C0") for l in summary_df["layout"]])
    ax.set_ylabel("Mean applied yaw, T1 [deg]")
    ax.set_xlabel("Run")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "yaw_by_port.png", dpi=150)
    plt.close(fig)


def plot_DEL_layout_comparison(all_data: dict, measured_del: dict, run_number: dict,
                                layout: str, load_type: str, fig_dir: Path):
    """One figure for one (layout, load_type) combo -- e.g. layout A / blade -- with
    the FLORIDyn-predicted DEL per turbine (from CLC_history) on the left and the
    real MLife/FAST.Farm-measured DEL per turbine (see load_measured_DEL) on the
    right, restricted to the ports belonging to `layout` (via PORT_LAYOUT). All
    turbines for the same run share one x position (color distinguishes the
    turbine), mirroring the previous combined plots' style.

    The predicted side has no baseline (the baseline port never runs FLORIDyn, so
    it has no CLC_history); the measured side includes this layout's baseline port
    if it's present in measured_del, placed after the runs and labeled 'BL-<layout>'.
    """
    turbine_colors = {1: "C0", 2: "C1", 3: "C2"}

    pred_ports = sorted(
        (p for p, d in all_data.items()
         if PORT_LAYOUT.get(p) == layout
         and d.get("clc_history") is not None and not d["clc_history"].empty),
        key=lambda p: run_number[p],
    )
    meas_ports = sorted(
        (p for p in measured_del if PORT_LAYOUT.get(p) == layout),
        key=lambda p: (p not in run_number, run_number.get(p, p)),
    )
    if not pred_ports and not meas_ports:
        return

    baseline_label = f"BL-{layout}"
    fig, (ax_pred, ax_meas) = plt.subplots(1, 2, figsize=(15, 6.5))
    plotted_values = []  # collected across both panels -> shared y-axis at the end

    # --- Left: FLORIDyn-predicted ---
    turbine_ids = sorted({
        int(col.split("_")[0][1:])
        for port in pred_ports
        for col in all_data[port]["clc_history"].columns
        if col.startswith("T") and col.endswith(f"_DEL_{load_type}")
    })
    for i, port in enumerate(pred_ports):
        clc = all_data[port]["clc_history"]
        for t in turbine_ids:
            col = f"T{t}_DEL_{load_type}"
            if col not in clc.columns:
                continue
            val = predicted_mean(clc, col, label=f"CLC_history_{port}")
            plotted_values.append(val)
            ax_pred.scatter(i, val, color=turbine_colors.get(t, f"C{t}"),
                             label=f"T{t}" if i == 0 else None, s=60, zorder=3,
                             edgecolors="k", linewidths=0.5)
    ax_pred.set_xticks(range(len(pred_ports)))
    ax_pred.set_xticklabels([str(run_number[p]) for p in pred_ports])
    ax_pred.set_xlabel("Run")
    ax_pred.set_ylabel(f"Mean {load_type} DEL [kNm]")
    ax_pred.set_title(f"FLORIDyn-predicted {load_type} DEL "
                      f"(mean {T_CUT_S:.0f}-{T_END_S:.0f} s)")
    ax_pred.legend(title="Turbine")
    ax_pred.grid(axis="y", alpha=0.3)

    # --- Right: MLife/FAST.Farm-measured ---
    xticklabels = []
    for i, port in enumerate(meas_ports):
        for t, vals in sorted(measured_del[port].items()):
            if load_type not in vals:
                continue
            plotted_values.append(vals[load_type])
            ax_meas.scatter(i, vals[load_type], color=turbine_colors.get(t, f"C{t}"),
                             label=f"T{t}" if i == 0 else None, s=60, zorder=3,
                             edgecolors="k", linewidths=0.5)
        xticklabels.append(str(run_number[port]) if port in run_number else baseline_label)
    ax_meas.set_xticks(range(len(meas_ports)))
    ax_meas.set_xticklabels(xticklabels)
    ax_meas.set_xlabel("Run")
    ax_meas.set_ylabel(f"Measured {load_type} DEL [kNm]")
    ax_meas.set_title(f"FAST.Farm-measured {load_type} DEL")
    ax_meas.legend(title="Turbine")
    ax_meas.grid(axis="y", alpha=0.3)

    # Shared y-axis across both panels, so predicted and measured DEL magnitudes
    # are directly comparable at a glance instead of each panel auto-scaling on
    # its own (which can visually exaggerate or hide a gap between the two).
    if plotted_values:
        lo, hi = min(plotted_values), max(plotted_values)
        pad = 0.05 * (hi - lo) if hi > lo else max(abs(hi), 1.0) * 0.05
        ylim = (lo - pad, hi + pad)
        ax_pred.set_ylim(ylim)
        ax_meas.set_ylim(ylim)

    fig.suptitle(f"{load_type.capitalize()} DEL, layout {layout}: predicted vs. measured")
    fig.tight_layout()
    fig.savefig(fig_dir / f"DEL_{load_type}_layout_{layout}.png", dpi=150)
    plt.close(fig)


def measured_farm_sum(measured_del: dict, port: int, load_type: str) -> float:
    """Sum over the 3 turbines of the MLife DEL of one port; NaN if a turbine is missing."""
    vals = [measured_del.get(port, {}).get(t, {}).get(load_type)
            for t in range(1, N_TURBINES + 1)]
    return float(sum(vals)) if all(v is not None for v in vals) else np.nan


def plot_DEL_farm_sum(summary_df: pd.DataFrame, measured_del: dict, run_number: dict,
                      baseline_ports: dict, fig_dir: Path):
    """Farm sum (T1+T2+T3) of the DEL over the whole run after T_CUT_S: MLife on the
    cut FAST.Farm .out vs. the sum of FLORIDyn's time-weighted means over the same
    span. One figure per (load type, layout) -- DEL_farm_sum_<load>_layout_<A|B>.png --
    with that layout's runs side by side. Each layout's no-yaw baseline is drawn as
    two dashed horizontal lines in the bar colors (orange = FAST.Farm/MLife measured,
    blue = FLORIDyn predicted, the latter only if given in BASELINE_PRED_DEL_SUM).
    Y axis starts at 0 and is shared between layouts A and B; unit is kNm."""
    for load_type in ("tower", "blade"):
        pred_col = f"mean_DEL_{load_type}_sum"
        meas_col = f"meas_DEL_{load_type}_sum"
        if pred_col not in summary_df.columns or meas_col not in summary_df.columns:
            continue

        # Baseline values per layout (measured from MLife, predicted from the constant).
        bl = {}
        for layout, port in baseline_ports.items():
            bl[layout] = {
                "meas": measured_farm_sum(measured_del, port, load_type) if port else np.nan,
                "pred": BASELINE_PRED_DEL_SUM.get(layout, {}).get(load_type, np.nan),
            }

        # Shared upper y limit across both layouts so A and B are directly comparable.
        all_vals = np.concatenate([
            summary_df[pred_col].to_numpy(dtype=float),
            summary_df[meas_col].to_numpy(dtype=float),
            np.array([v for d in bl.values() for v in d.values()], dtype=float),
        ])
        all_vals = all_vals[np.isfinite(all_vals)]
        if all_vals.size == 0:
            continue
        ymax = 1.12 * all_vals.max()

        for layout in baseline_ports:
            rows = summary_df[summary_df["layout"] == layout]
            if rows.empty:
                continue
            labels = [str(r) for r in rows["run"]]
            pred = rows[pred_col].to_numpy(dtype=float)
            meas = rows[meas_col].to_numpy(dtype=float)
            x = np.arange(len(labels))

            fig, (ax, ax_pct) = plt.subplots(2, 1, figsize=(max(8, 1.3 * len(labels) + 3), 9),
                                             sharex=True, gridspec_kw={"height_ratios": [3, 1]})
            ax.bar(x - 0.2, meas, width=0.4, color="C1", label="FAST.Farm/MLife (sum T1-T3)")
            ax.bar(x + 0.2, pred, width=0.4, color="C0",
                   label="FLORIDyn predicted (sum of T1-T3 means)")
            if np.isfinite(bl[layout]["meas"]):
                ax.axhline(bl[layout]["meas"], color="C1", linestyle="--", linewidth=2,
                           label=f"Baseline {layout} measured")
            if np.isfinite(bl[layout]["pred"]):
                ax.axhline(bl[layout]["pred"], color="C0", linestyle="--", linewidth=2,
                           label=f"Baseline {layout} predicted")
            ax.set_ylim(0, ymax)
            ax.set_ylabel(f"Farm {load_type} DEL sum [kNm]")
            ax.set_title(f"Farm {load_type} DEL sum, layout {layout}")
            ax.legend(fontsize=20, loc="lower right")
            ax.grid(axis="y", alpha=0.3)

            with np.errstate(divide="ignore", invalid="ignore"):
                pct = np.nan_to_num((pred - meas) / meas * 100)
            ax_pct.bar(x, pct, color=["C2" if p >= 0 else "C3" for p in pct])
            ax_pct.axhline(0, color="k", linewidth=0.8)
            ax_pct.set_ylabel("Error [%]")
            ax_pct.set_xticks(x)
            ax_pct.set_xticklabels(labels)
            ax_pct.set_xlabel("Run")
            ax_pct.grid(axis="y", alpha=0.3)

            fig.tight_layout()
            fig.savefig(fig_dir / f"DEL_farm_sum_{load_type}_layout_{layout}.png", dpi=150)
            plt.close(fig)


def plot_DEL_farm_sum_combined(summary_df: pd.DataFrame, measured_del: dict, baseline_ports: dict,
                               fig_dir: Path):
    """DEL_farm_sum_<tower|blade>.png: one bar chart with both layouts side by side,
    order 1-5, BL-A, 6-10, BL-B (FAST.Farm/MLife vs FLORIDyn predicted farm sum T1-T3),
    with the (pred-meas)/meas error panel below. No title (the caption is in the text);
    the vertical axis of the upper panel does not start from zero."""
    for load_type in ("tower", "blade"):
        pred_col = f"mean_DEL_{load_type}_sum"
        meas_col = f"meas_DEL_{load_type}_sum"
        if pred_col not in summary_df.columns or meas_col not in summary_df.columns:
            continue
        labels, pred, meas = [], [], []
        for layout, bl_port in baseline_ports.items():
            rows = summary_df[summary_df["layout"] == layout]
            for _, r in rows.iterrows():
                labels.append(str(r["run"]))
                pred.append(float(r[pred_col]))
                meas.append(float(r[meas_col]))
            labels.append(f"BL-{layout}")
            pred.append(BASELINE_PRED_DEL_SUM.get(layout, {}).get(load_type, np.nan))
            meas.append(measured_farm_sum(measured_del, bl_port, load_type) if bl_port else np.nan)
        pred = np.array(pred, dtype=float)
        meas = np.array(meas, dtype=float)
        finite = np.concatenate([pred[np.isfinite(pred)], meas[np.isfinite(meas)]])
        if finite.size == 0:
            continue
        lo, hi = finite.min(), finite.max()
        pad = 0.05 * (hi - lo) if hi > lo else 0.05 * hi
        x = np.arange(len(labels))
        fig, (ax, ax_pct) = plt.subplots(2, 1, figsize=(max(12, 1.3 * len(labels) + 3), 10),
                                         sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        ax.bar(x - 0.2, meas, width=0.4, color="C1", label="FAST.Farm/MLife (sum T1-T3)")
        ax.bar(x + 0.2, pred, width=0.4, color="C0",
               label="FLORIDyn predicted (sum of T1-T3 means)")
        ax.set_ylim(lo - pad, hi + 7 * pad)
        ax.set_ylabel(f"Farm {load_type} DEL sum [kNm]")
        ax.legend(fontsize=20, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.nan_to_num((pred - meas) / meas * 100)
        ax_pct.bar(x, pct, color=["C2" if p >= 0 else "C3" for p in pct])
        ax_pct.axhline(0, color="k", linewidth=0.8)
        ax_pct.set_ylabel("Error [%]")
        ax_pct.set_xticks(x)
        ax_pct.set_xticklabels(labels)
        ax_pct.set_xlabel("Run")
        ax_pct.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / f"DEL_farm_sum_{load_type}.png", dpi=150)
        plt.close(fig)


def plot_power(summary_df: pd.DataFrame, fig_dir: Path, baseline_power_w: dict = None):
    """Mean FAST.Farm-measured farm power per run, one figure per layout
    (power_by_port_layout_<A|B>.png), with that layout's no-yaw baseline as a
    dashed horizontal line. Y axis starts at 0, shared between layouts."""
    baseline_power_w = baseline_power_w or {}
    ymax = summary_df["mean_P_FASTfarm_W"].max() / 1e6
    if baseline_power_w:
        ymax = max(ymax, max(baseline_power_w.values()) / 1e6)
    ymax *= 1.12
    for layout in ("A", "B"):
        rows = summary_df[summary_df["layout"] == layout]
        if rows.empty:
            continue
        fig, ax = plt.subplots(figsize=(4.2, 6))
        ax.bar(rows["run"].astype(str), rows["mean_P_FASTfarm_W"] / 1e6,
               color=LAYOUT_COLORS[layout])
        if layout in baseline_power_w:
            ax.axhline(baseline_power_w[layout] / 1e6, color="C1", linestyle="--",
                       linewidth=2, label=f"Baseline {layout}")
            ax.legend(loc="lower right")
        ax.set_ylim(0, ymax)
        ax.set_title(f"Layout {layout}")
        ax.set_ylabel("Mean farm power [MW]")
        ax.set_xlabel("Run")
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / f"power_by_port_layout_{layout}.png", dpi=150)
        plt.close(fig)


def plot_power_vs_baseline(summary_df: pd.DataFrame, baseline_power_w: dict, fig_dir: Path):
    """Each port's mean FAST.Farm power against its own layout's no-yaw-steering
    baseline (BASELINE_PORTS), one figure per layout
    (power_vs_baseline_layout_<A|B>.png): absolute MW with the layout's dashed
    baseline line, and the same numbers normalized to that baseline
    (1.0 = baseline, >1.0 = gain).

    `power_ratio_vs_baseline` must already be in summary_df (computed per-row against
    that row's own layout in main())."""
    if "power_ratio_vs_baseline" not in summary_df.columns:
        return

    ymax = summary_df["mean_P_FASTfarm_W"].max() / 1e6
    if baseline_power_w:
        ymax = max(ymax, max(baseline_power_w.values()) / 1e6)
    ymax *= 1.12

    for layout in ("A", "B"):
        rows = summary_df[summary_df["layout"] == layout]
        if rows.empty:
            continue
        runs = rows["run"].astype(str).tolist()
        powers = rows["mean_P_FASTfarm_W"].to_numpy()
        ratios = rows["power_ratio_vs_baseline"].to_numpy()

        fig, axes = plt.subplots(1, 2, figsize=(7.5, 6.5))

        ax = axes[0]
        ax.bar(runs, powers / 1e6, color=LAYOUT_COLORS[layout])
        if layout in baseline_power_w:
            ax.axhline(baseline_power_w[layout] / 1e6, color="C1", linestyle="--",
                       linewidth=2, label=f"Baseline layout {layout}")
            ax.legend(fontsize=20, loc="lower right")
        ax.set_ylim(0, ymax)
        ax.set_ylabel("Mean farm power [MW]", fontsize=28)
        ax.set_xlabel("Run", fontsize=28)
        ax.set_title(f"Layout {layout}")
        ax.grid(axis="y", alpha=0.3)
        ax.tick_params(axis="both", labelsize=18)

        ax = axes[1]
        ratio_colors = ["C2" if r >= 1.0 else "C3" for r in ratios]
        ax.bar(runs, ratios, color=ratio_colors)
        ax.axhline(1.0, color="C1", linestyle="--", linewidth=2, label=f"Baseline {layout}")
        rmax = np.nanmax(ratios) if np.isfinite(ratios).any() else 1.0
        ax.set_ylim(0, max(1.0, rmax) * 1.12)
        ax.set_ylabel("Power / own baseline [-]", fontsize=28)
        ax.set_xlabel("Run", fontsize=28)
        ax.set_title(f"Layout {layout}")
        ax.legend(loc="lower right")
        ax.grid(axis="y", alpha=0.3)
        ax.tick_params(axis="both", labelsize=18)

        fig.tight_layout()
        fig.savefig(fig_dir / f"power_vs_baseline_layout_{layout}.png", dpi=150)
        plt.close(fig)


def _read_fastfarm_out_columns(out_path: Path, columns: list) -> pd.DataFrame:
    """Read selected columns from a FAST.Farm .out text file.

    These files have a few description lines, then a tab-separated header row
    (starts with "Time"), then a units row, then tab-separated data. Mirrors the
    header-detection used by controller_floris/FLORIDyn_<port>.py so column
    indices line up the same way.
    """
    with open(out_path, "r", encoding="utf-8", errors="ignore") as f:
        head_lines = [f.readline() for _ in range(12)]

    header_idx = next(
        (i for i, ln in enumerate(head_lines) if ln.strip().startswith("Time")), None
    )
    if header_idx is None:
        raise ValueError(f"Could not find header row in {out_path}")

    raw_hdr = head_lines[header_idx].strip()
    header_cols = raw_hdr.split("\t") if "\t" in raw_hdr else re.split(r"\s+", raw_hdr)
    data_start = header_idx + 2  # skip header row + units row

    # Select by integer position, not by name: FAST.Farm's own header has genuine
    # duplicate names (e.g. "YawBrFxp" appears twice), which breaks pandas' `names=`.
    missing = [c for c in columns if c not in header_cols]
    if missing:
        raise ValueError(f"{out_path.name} is missing columns {missing}")
    col_positions = [header_cols.index(c) for c in columns]

    def _read(sep):
        d = pd.read_csv(
            out_path, sep=sep, skiprows=data_start, header=None,
            usecols=col_positions, engine="c",
        )
        return d[col_positions].set_axis(columns, axis=1)  # restore requested order + names

    df = _read("\t")
    if df.shape[1] != len(columns):  # fall back if the file wasn't actually tab-delimited
        df = _read(r"\s+")
    return df


def load_fastfarm_wind(fastfarm_dir: Path, port: int, turbine: int = 1) -> pd.DataFrame:
    """Reconstruct the wind speed/direction time series FLORIDyn's EnKF actually saw,
    from the raw <port>.T<turbine>.out file, using the same filter + trailing-window
    logic as controller_floris/FLORIDyn_<port>.py's _average_FF_signals():

      speed [m/s] = RMS of the disk-resultant velocity (WindDiskVelX/Y/Z) over a
                    trailing 50s window
      dir   [deg] = -mean(WindHubAngXY over the trailing 50s window) + 270.0

    Difference from the live pipeline: the live code re-filters a fresh 500-sample
    window on every call; here the Butterworth filter is applied once to the whole
    series and then a rolling window is taken, which is numerically very close but
    avoids repeated edge artifacts and is far cheaper to compute after the fact.

    Returns the full [0, 4600]s series -- nothing is cut or shifted (see the
    update_rate fix note near the top of this file for why that's no longer needed).
    """
    out_path = fastfarm_dir / f"{port}.T{turbine}.out"
    if not out_path.exists():
        raise FileNotFoundError(out_path)

    cols = ["Time", "WindHubAngXY", "WindDiskVelX", "WindDiskVelY", "WindDiskVelZ"]
    df = _read_fastfarm_out_columns(out_path, cols).apply(pd.to_numeric, errors="coerce").dropna()
    df = df.sort_values("Time").reset_index(drop=True)

    b, a = butter(FF_ORDER, FF_CUTOFF / (0.5 * FF_FS), btype="low")
    padlen = 3 * (max(len(a), len(b)) - 1)
    def lowpass(x):
        return filtfilt(b, a, x) if len(x) > padlen else x

    disk_res = np.sqrt(df["WindDiskVelX"] ** 2 + df["WindDiskVelY"] ** 2 + df["WindDiskVelZ"] ** 2)
    disk_res_f = lowpass(disk_res.to_numpy())
    wind_hub_ang_f = lowpass(df["WindHubAngXY"].to_numpy())

    window_n = max(1, int(round(FF_WINDOW_S * FF_FS)))
    speed_rms = pd.Series(disk_res_f ** 2).rolling(window_n, min_periods=1).mean() ** 0.5
    dir_mean = pd.Series(wind_hub_ang_f).rolling(window_n, min_periods=1).mean()

    return pd.DataFrame({
        "Time": df["Time"].to_numpy(),
        "ws_ff_ms": speed_rms.to_numpy(),
        "wd_ff_deg": -dir_mean.to_numpy() + 270.0,
    })


def load_measured_DEL(xlsx_path: Path) -> dict:
    """Per-(port, turbine) tower/blade DEL from an MLife-generated 'Short-term DELs'
    sheet (e.g. fastfarm_simulation/Mlife/CertTest/FLORIDyn_zeroMQ_GPR_Short-term.xlsx):
    for each of the two load-channel blocks in that sheet (Tower_Base -> tower,
    Blade_Root -> blade), reduces that row's 6 channel columns (Tower_1..6 /
    Blade_1..6) according to MEASURED_CHANNEL ("max" = max across them, "last" =
    Tower_6 / Blade_6). Block boundaries are located by their title text, not hardcoded row
    numbers, so this tolerates the sheet growing as more runs are added.

    Returns {port: {turbine: {"tower": val, "blade": val}}}.
    """
    raw = pd.read_excel(xlsx_path, sheet_name="Short-term DELs", header=None)
    fname_re = re.compile(r"^(\d+)\.T(\d+)$")
    result: dict = {}

    def _parse_block(block_title: str, channel_prefix: str, key: str):
        title_rows = raw.index[raw[0].astype(str).str.startswith(block_title, na=False)]
        if len(title_rows) == 0:
            raise ValueError(f"Could not find a '{block_title}' block in {xlsx_path.name}")
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
            port, turbine = (int(g) for g in fname_re.match(fname).groups())
            value = reduce_channels([raw.iat[row, c] for c in channel_cols])
            if value is not None:
                result.setdefault(port, {}).setdefault(turbine, {})[key] = value
            row += 1

    _parse_block("Tower_Base", "Tower_", "tower")
    _parse_block("Blade_Root", "Blade_", "blade")
    return result


def load_baseline_power(fastfarm_dir: Path, baseline_port: int,
                         n_turbines: int = N_TURBINES) -> float:
    """Total farm power [W] for the baseline (greedy, no yaw steering) run, read
    straight from the raw <baseline_port>.T<n>.out files' GenPwr channel [kW].
    Averaged over the full run (nothing discarded) and summed across turbines."""
    total_w = 0.0
    for tid in range(1, n_turbines + 1):
        out_path = fastfarm_dir / f"{baseline_port}.T{tid}.out"
        if not out_path.exists():
            raise FileNotFoundError(out_path)
        df = _read_fastfarm_out_columns(out_path, ["Time", "GenPwr"])
        df = df.apply(pd.to_numeric, errors="coerce").dropna()
        if df.empty:
            raise ValueError(f"{out_path.name} has no usable GenPwr data")
        total_w += df["GenPwr"].mean() * 1e3  # kW -> W
    return total_w


def _farm_power_from_out(fastfarm_dir: Path, port: int) -> float:
    try:
        return load_baseline_power(fastfarm_dir, port)
    except (FileNotFoundError, ValueError) as e:
        print(f"  [{port}] could not load farm power from .out: {e}")
        return float("nan")


def _plot_wind_grid(ports: list, all_data: dict, fastfarm_wind: dict, run_number: dict,
                     enkf_col: str, ff_col: str, ylabel: str, title: str,
                     out_path: Path):
    """One metric (ws or wd), one small subplot per run in a 2-column grid, instead
    of every run stacked as a full-width row -- with 10 runs a single-column stack
    gets very tall and hard to scan; a grid reads in about a third the height."""
    ncols = 2
    nrows = -(-len(ports) // ncols)  # ceil
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.8 * nrows), sharex=True)
    axes = np.atleast_2d(axes)

    for i, port in enumerate(ports):
        ax = axes[i // ncols][i % ncols]
        enkf = all_data[port]["enkf_log"]
        ax.plot(enkf["time_s"], enkf[enkf_col], color="C0", label="EnKF")

        ff = fastfarm_wind.get(port)
        if ff is not None:
            ax.plot(ff["Time"], ff[ff_col], color="C1", alpha=0.8, label="FAST.Farm (T1)")
        else:
            ax.text(0.5, 0.5, "no .out data", transform=ax.transAxes, ha="center", color="gray")

        ax.set_title(f"Run {run_number[port]}", fontsize=24)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=20)

    for i in range(len(ports), nrows * ncols):  # hide unused trailing axes
        axes[i // ncols][i % ncols].axis("off")

    for c in range(ncols):
        axes[-1][c].set_xlabel("Time [s]")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_wind_conditions(all_data: dict, fastfarm_wind: dict, run_number: dict, fig_dir: Path):
    """EnKF wind estimate vs. the real wind reconstructed from FAST.Farm's own T1
    output (see load_fastfarm_wind), as two separate grids (speed, direction) rather
    than one combined figure -- keeps each grid a reasonable size regardless of how
    many runs there are.

    Both curves are already on the same clock: FLORIDyn's own logs (enkf_log) never
    contain the pre-1000s transient in the first place, and the FAST.Farm curve has
    been cut-and-shifted onto that same clock in load_fastfarm_wind. No further
    trimming needed here."""
    ports = [p for p, d in all_data.items() if d.get("enkf_log") is not None]
    if not ports:
        return
    ports = sorted(ports, key=lambda p: run_number[p])

    _plot_wind_grid(ports, all_data, fastfarm_wind, run_number,
                     "ws_hat", "ws_ff_ms", "ws [m/s]",
                     "EnKF estimate vs. real FAST.Farm wind speed (from T1 .out), per run",
                     fig_dir / "wind_speed_by_port.png")
    _plot_wind_grid(ports, all_data, fastfarm_wind, run_number,
                     "wd_hat", "wd_ff_deg", "wd [deg]",
                     "EnKF estimate vs. real FAST.Farm wind direction (from T1 .out), per run",
                     fig_dir / "wind_direction_by_port.png")


def plot_prediction_accuracy(all_data: dict, run_number: dict, fig_dir: Path):
    """Predicted vs FAST.Farm-measured power at each reliability check, one small
    subplot per run instead of every run overlaid on a single axes -- with 10 runs
    the overlay just turns into a wall of overlapping colors, so a small-multiples
    grid (all sharing the same axis limits, for an at-a-glance comparison) reads
    much better."""
    ports = [p for p, d in all_data.items()
             if d.get("reliability_log") is not None and not d["reliability_log"].empty]
    if not ports:
        return
    ports = sorted(ports, key=lambda p: run_number[p])

    # Shared axis limits across all subplots, so the perfect-prediction diagonal
    # lands in the same place everywhere and run-to-run spread is comparable.
    all_pred  = pd.concat([all_data[p]["reliability_log"]["P_pred_W"] for p in ports]) / 1e6
    all_meas  = pd.concat([all_data[p]["reliability_log"]["P_FASTfarm_W"] for p in ports]) / 1e6
    lo = min(all_pred.min(), all_meas.min())
    hi = max(all_pred.max(), all_meas.max())
    pad = 0.05 * (hi - lo) if hi > lo else 1.0
    lims = (lo - pad, hi + pad)

    ncols = 2
    nrows = -(-len(ports) // ncols)  # ceil
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 5.5 * nrows), squeeze=False)

    for i, port in enumerate(ports):
        ax = axes[i // ncols][i % ncols]
        rel = all_data[port]["reliability_log"]
        ax.scatter(rel["P_pred_W"] / 1e6, rel["P_FASTfarm_W"] / 1e6, alpha=0.7, s=25)
        ax.plot(lims, lims, "k--", linewidth=1, label="perfect prediction")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_title(f"Run {run_number[port]}")
        ax.set_xlabel("Predicted power [MW]")
        ax.set_ylabel("FAST.Farm-measured power [MW]")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=20)

    for i in range(len(ports), nrows * ncols):  # hide unused trailing axes
        axes[i // ncols][i % ncols].axis("off")

    fig.suptitle("Predicted vs measured power at each reliability check, per run")
    fig.tight_layout()
    fig.savefig(fig_dir / "prediction_accuracy.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                         help="Folder containing CLC_history_*, enkf_log_*, reliability_log_* files")
    parser.add_argument("--fastfarm-dir", default=DEFAULT_FASTFARM_DIR,
                         help="Folder containing <port>.T1.out etc. (for the real wind overlay)")
    parser.add_argument("--ports", nargs="+", type=int, default=sorted(LAMBDA_WEIGHTS),
                         help="Ports to include (default: 15013-15017)")
    parser.add_argument("--baseline-port-a", type=int, default=BASELINE_PORTS["A"],
                         help="No-yaw-steering baseline port for layout A (15013-15017) "
                              "power normalization (default: 15006). Pass 0 to skip.")
    parser.add_argument("--baseline-port-b", type=int, default=BASELINE_PORTS["B"],
                         help="No-yaw-steering baseline port for layout B (15018-15022, "
                              "mirrored) power normalization (default: 15012). Pass 0 to skip.")
    parser.add_argument("--mlife-xlsx", default=DEFAULT_MLIFE_XLSX,
                         help="MLife 'Short-term DELs' workbook for the measured-DEL "
                              "plot. Pass '' to skip it.")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        raise SystemExit(
            f"--results-dir does not exist or is not a directory: {results_dir}\n"
            f"(hint: WSL-style '/mnt/d/...' paths only resolve under WSL -- pass a path "
            f"this Python interpreter can actually see, e.g. 'D:\\...' on native Windows)"
        )

    fastfarm_dir = Path(args.fastfarm_dir)
    if not fastfarm_dir.is_dir():
        print(f"[warn] --fastfarm-dir does not exist or is not a directory: {fastfarm_dir}\n"
              f"       wind overlay plots will be empty for every port")

    out_dir = Path(__file__).parent
    fig_dir = out_dir / "figures_GPR"
    fig_dir.mkdir(exist_ok=True)

    # Progressive run number in place of the raw port number, in the order ports
    # were requested (default: ascending port number) -- used on every plot/legend.
    run_number = {port: i + 1 for i, port in enumerate(args.ports)}

    print(f"Reading results from: {results_dir}")
    all_data = {}
    fastfarm_wind = {}
    summary_rows = []
    for port in args.ports:
        print(f"Port {port}:")
        data = load_port_data(results_dir, port)
        all_data[port] = data
        summary_rows.append(summarize_port(port, data))
        try:
            fastfarm_wind[port] = load_fastfarm_wind(fastfarm_dir, port)
        except (FileNotFoundError, ValueError) as e:
            print(f"  [{port}] could not load FAST.Farm wind data: {e}")

    measured_del = None
    if args.mlife_xlsx:
        mlife_path = Path(args.mlife_xlsx)
        try:
            measured_del = load_measured_DEL(mlife_path)
            print(f"\nLoaded measured DEL from: {mlife_path} ({len(measured_del)} ports)")
            # The GPR workbooks have no baseline ports: take the measured baselines (same FAST.Farm runs)
            # from the original FFNN-campaign workbook.
            bl_path = Path(args.mlife_xlsx.replace("_GPR", ""))
            try:
                bl_meas = load_measured_DEL(bl_path)
                for bp in BASELINE_PORTS.values():
                    if bp in bl_meas:
                        measured_del[bp] = bl_meas[bp]
                print(f"Baseline measured DEL from: {bl_path}")
            except (FileNotFoundError, ValueError) as e:
                print(f"[warn] could not load baseline measured DEL from {bl_path}: {e}")
        except (FileNotFoundError, ValueError) as e:
            print(f"\n[warn] could not load measured DEL from {mlife_path}: {e}")

    baseline_ports = {"A": args.baseline_port_a, "B": args.baseline_port_b}
    baseline_power_w = {}
    for layout, baseline_port in baseline_ports.items():
        if not baseline_port:
            continue
        try:
            baseline_power_w[layout] = load_baseline_power(fastfarm_dir, baseline_port)
            print(f"\nBaseline layout {layout} (port {baseline_port}, no yaw steering) "
                  f"farm power: {baseline_power_w[layout] / 1e6:.4f} MW")
        except (FileNotFoundError, ValueError) as e:
            print(f"\n[warn] could not load baseline power for layout {layout} "
                  f"from port {baseline_port}: {e}")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.insert(0, "run", summary_df["port"].map(run_number))
    # Farm power of every port from the FAST.Farm .out files (GenPwr, T1..T3 summed, same window and
    # same function as the baselines). The reliability_log value is kept only for reference.
    if "mean_P_FASTfarm_W" in summary_df.columns:
        summary_df["mean_P_FASTfarm_W_rellog"] = summary_df["mean_P_FASTfarm_W"]
    summary_df["mean_P_FASTfarm_W"] = summary_df["port"].map(
        lambda p: _farm_power_from_out(fastfarm_dir, p))
    if baseline_power_w and "mean_P_FASTfarm_W" in summary_df.columns and "layout" in summary_df.columns:
        summary_df["power_ratio_vs_baseline"] = summary_df.apply(
            lambda row: row["mean_P_FASTfarm_W"] / baseline_power_w[row["layout"]]
            if row["layout"] in baseline_power_w else np.nan,
            axis=1,
        )
    if measured_del:
        for load_type in ("tower", "blade"):
            meas_col = f"meas_DEL_{load_type}_sum"
            summary_df[meas_col] = summary_df["port"].map(
                lambda p: measured_farm_sum(measured_del, p, load_type))
            pred_col = f"mean_DEL_{load_type}_sum"
            if pred_col in summary_df.columns:
                summary_df[f"pct_DEL_{load_type}_sum_pred_vs_meas"] = (
                    (summary_df[pred_col] - summary_df[meas_col]) / summary_df[meas_col] * 100)

    summary_csv = out_dir / "summary_by_port_GPR.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\nWrote summary: {summary_csv}")
    print(summary_df.to_string(index=False))

    if "mean_yaw_T1_deg" in summary_df.columns:
        plot_yaw(summary_df, fig_dir)
    for layout in ("A", "B"):
        for load_type in ("blade", "tower"):
            plot_DEL_layout_comparison(all_data, measured_del or {}, run_number,
                                        layout, load_type, fig_dir)
    if measured_del:
        plot_DEL_farm_sum(summary_df, measured_del, run_number, baseline_ports, fig_dir)
        plot_DEL_farm_sum_combined(summary_df, measured_del, baseline_ports, fig_dir)
    if "mean_P_FASTfarm_W" in summary_df.columns:
        plot_power(summary_df, fig_dir, baseline_power_w)
        if baseline_power_w:
            plot_power_vs_baseline(summary_df, baseline_power_w, fig_dir)
    plot_wind_conditions(all_data, fastfarm_wind, run_number, fig_dir)
    plot_prediction_accuracy(all_data, run_number, fig_dir)

    print(f"Wrote figures to: {fig_dir}")


if __name__ == "__main__":
    main()
