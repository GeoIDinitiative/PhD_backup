"""
run_pipeline.py  —  single entry point for the whole updated Etna tilt workflow
==============================================================================
Runs every stage end-to-end and emits the standard outputs (the modern equivalents
of the old swcc_edit.py + comprehensive_signal_analysis_complete.py deliverables).

Stages (each lives in its own validated module; this script orchestrates them):
  1. build_clean_bandpassed        raw(P-wave-removed) → despike+detrend+per-seg bandpass
  2. swcc_comprehensive            vectorised segment-aware SWCC, both components
  3. flag_significant_peaks        phase-randomised null floor → significance
  4. characterize_significant_peaks  synchrony / volcanic coincidence / timeline
  5. credibility_checks            filter response, de-ringing, null test
Then this script adds the old-style figures:
  6. per (station,component) SWCC overview + cumulative (|r| vs time, threshold + null floor,
     volcanic overlay for INGV)
  7. INGV-vs-experiment dataset comparison (|r|, SNR)
  8. candidate panels: waveform + spectrogram + Morlet wavelet for each significant detection

Deliberately DROPPED (obsolete under the new design — documented in MANIFEST):
  · P-wave contamination-probability / impact plots and contaminated-vs-clean peak lists
    — P-waves are excised upfront, so there is no contaminated/clean split; the
      phase-randomised null floor is the modern significance criterion instead.

Usage:
  python3 run_pipeline.py                 # full run (recompute + all figures)
  python3 run_pipeline.py --plots-only    # skip recompute, regenerate figures from existing data
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import spectrogram

import build_clean_bandpassed
import swcc_comprehensive
import flag_significant_peaks
import swcc_template4
import characterize_significant_peaks
import credibility_checks
import swcc_oldstyle_plots

try:
    import pywt
    HAVE_PYWT = True
except Exception:
    HAVE_PYWT = False

BASE  = Path("/home/owen/etna_signals_phd")
CLEAN = BASE / "clean_bandpassed"
SWCC  = BASE / "SWCC_comprehensive"
VOLC  = BASE / "etna_volcanic_events_cleaned.csv"
THRESHOLD = swcc_comprehensive.THRESHOLD
PERIODS = {"ingv": ("2022-11-14", "2023-03-01"),
           "experiment": ("2023-07-23", "2023-08-03")}


# ── stages 1–5 ────────────────────────────────────────────────────────────────
def run_compute():
    print("\n" + "="*70 + "\n[1/5] clean bandpass\n" + "="*70)
    build_clean_bandpassed.main()
    print("\n" + "="*70 + "\n[2/5] comprehensive SWCC\n" + "="*70)
    swcc_comprehensive.main()
    print("\n" + "="*70 + "\n[3/6] significance (null floor, T1-3)\n" + "="*70)
    flag_significant_peaks.main()
    print("\n" + "="*70 + "\n[4/6] template4 long-segment SWCC (separate, then merged)\n" + "="*70)
    swcc_template4.main()
    print("\n" + "="*70 + "\n[5/6] characterization\n" + "="*70)
    characterize_significant_peaks.main()
    print("\n" + "="*70 + "\n[6/6] credibility checks\n" + "="*70)
    credibility_checks.main()


# ── helpers ───────────────────────────────────────────────────────────────────
def load_volc():
    v = pd.read_csv(VOLC)
    v["start"] = pd.to_datetime(v["Date"] + " " + v["Time"], errors="coerce")
    return v.dropna(subset=["start"])


def load_clean(ds, st, comp):
    f = CLEAN / ds / f"{st}_{comp}_0p001-0p01Hz_clean_bp.feather"
    if not f.exists():
        return None
    d = pd.read_feather(f); d["datetime"] = pd.to_datetime(d["datetime"])
    return d


# ── 6. per (station,component) SWCC overview + cumulative ─────────────────────
def swcc_overview_plots(peaks, volc):
    for (ds, st, comp), g in peaks.groupby(["dataset", "station", "component"]):
        g = g.sort_values("peak_time")
        floor = g["null_floor"].iloc[0]
        sig = g[g.significant]
        out = SWCC / ds / st; out.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
        # top: peak |r| vs time
        if ds == "ingv":
            for _, e in volc.iterrows():
                if pd.Timestamp(PERIODS[ds][0]) <= e.start <= pd.Timestamp(PERIODS[ds][1]):
                    ax[0].axvline(e.start, color="#16a34a", alpha=0.20, lw=1, zorder=1)
        ax[0].scatter(g.peak_time, g.abs_r, s=10, c="#9ca3af", alpha=0.5, label="all peaks", zorder=2)
        ax[0].scatter(sig.peak_time, sig.abs_r, s=45, c="#dc2626", zorder=3,
                      label=f"significant (n={len(sig)})")
        ax[0].axhline(THRESHOLD, ls="--", c="#6b7280", lw=1, label=f"old threshold {THRESHOLD}")
        ax[0].axhline(floor, ls="--", c="#16a34a", lw=2, label=f"null floor {floor:.2f}")
        ax[0].set_ylabel("peak |r|"); ax[0].legend(fontsize=8, ncol=2)
        ax[0].set_title(f"{ds} / {st} / {comp} — SWCC peaks"
                        + ("  (green = volcanic events)" if ds == "ingv" else ""))
        # bottom: cumulative significant count + cumulative |r|
        ax[1].plot(g.peak_time, np.arange(1, len(g)+1), color="#9333ea", lw=1.5,
                   label="cumulative all peaks")
        if len(sig):
            ax[1].plot(sig.peak_time, np.arange(1, len(sig)+1), color="#dc2626", lw=2.5,
                       label="cumulative significant")
        ax[1].set_ylabel("cumulative count"); ax[1].legend(fontsize=9)
        ax[1].set_xlabel("time"); ax[1].tick_params(axis="x", rotation=30)
        for a in ax: a.grid(alpha=0.3)
        fig.tight_layout(); fig.savefig(out / f"{st}_{comp}_swcc_overview.png", dpi=140)
        plt.close(fig)
    print(f"  swcc overviews → {SWCC}/<dataset>/<station>/")


# ── 7. dataset comparison ─────────────────────────────────────────────────────
def dataset_comparison(peaks):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for a, metric, lab in zip(ax, ["abs_r", "snr_db"], ["peak |r|", "SNR (dB)"]):
        data, labels = [], []
        for ds in ["ingv", "experiment"]:
            vals = peaks[peaks.dataset == ds][metric].dropna().to_numpy()
            if len(vals): data.append(vals); labels.append(f"{ds}\n(n={len(vals)})")
        if data:
            bp = a.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
            for patch, c in zip(bp["boxes"], ["#1f2937", "#dc2626"]):
                patch.set_facecolor(c); patch.set_alpha(0.4)
        a.set_ylabel(lab); a.grid(alpha=0.3)
        a.set_title(f"INGV vs experiment — {lab}")
    fig.tight_layout(); fig.savefig(SWCC / "dataset_comparison.png", dpi=140); plt.close(fig)
    print(f"  dataset comparison → {SWCC}/dataset_comparison.png")


# ── 7b. template / simulation comparison bar charts ──────────────────────────
def _topcorner_legend(ax, title, loc="upper right"):
    """Place legend in a top corner with opaque box + headroom so bars never clash."""
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.30)            # headroom above tallest bar
    ax.legend(loc=loc, title=title, framealpha=1.0, edgecolor="gray",
              fontsize=9, ncol=2)


def template_sim_bars(peaks):
    """Per-dataset grouped bar charts: peaks per station, grouped by template and by sim.
    (Modern equivalent of the old plot_station_template_grouped / _sim_grouped.)"""
    palette = plt.cm.tab10.colors
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for r, (dim, col, pre) in enumerate([("Template", "template", "T"),
                                         ("Simulation", "sim", "S")]):
        cats = sorted(peaks[col].dropna().unique())
        for c, ds in enumerate(["ingv", "experiment"]):
            ax = axes[r, c]
            sub = peaks[peaks.dataset == ds]
            stations = sorted(sub.station.unique())
            x = np.arange(len(stations)); n = len(cats); w = 0.8 / max(1, n)
            for j, cat in enumerate(cats):
                vals = [int(((sub.station == st) & (sub[col] == cat)).sum())
                        for st in stations]
                ax.bar(x + (j - (n-1)/2)*w, vals, w, zorder=2, color=palette[j],
                       label=f"{pre}{str(cat)[-1]}")
            ax.set_xticks(x); ax.set_xticklabels(stations, rotation=30, ha="right")
            ax.set_ylabel("number of peaks")
            ax.set_title(f"{ds}: peaks per station by {dim.lower()}")
            ax.grid(axis="y", alpha=0.3, zorder=1)
            _topcorner_legend(ax, dim)        # top-corner, never clashes
    fig.tight_layout()
    fig.savefig(SWCC / "template_sim_comparison.png", dpi=140); plt.close(fig)
    print(f"  template/sim bars → {SWCC}/template_sim_comparison.png")


# ── 8. candidate panels (waveform + spectrogram + Morlet wavelet) ─────────────
def candidate_panels(peaks):
    out = SWCC / "candidates"; out.mkdir(parents=True, exist_ok=True)
    sig = peaks[peaks.significant].copy()
    n = 0
    for _, p in sig.iterrows():
        d = load_clean(p.dataset, p.station, p.component)
        if d is None: continue
        seg = d[d.segment_id == p.segment_id]
        if len(seg) < 256: continue
        t = (seg.datetime - seg.datetime.iloc[0]).dt.total_seconds().to_numpy()
        x = seg.bandpassed.to_numpy()
        tpk = (p.peak_time - seg.datetime.iloc[0]).total_seconds()

        nrow = 3 if HAVE_PYWT else 2
        fig, ax = plt.subplots(nrow, 1, figsize=(12, 3*nrow))
        ax[0].plot(t/60, x, lw=0.7, color="#2563eb"); ax[0].axvline(tpk/60, color="#dc2626", ls="--")
        ax[0].set_ylabel("bandpassed"); ax[0].grid(alpha=0.3)
        ax[0].set_title(f"{p.dataset}/{p.station}/{p.component}  {p.sim}/{p.template}  "
                        f"|r|={p.abs_r:.2f}  {p.peak_time:%Y-%m-%d %H:%M}")
        f, tt, Sxx = spectrogram(x, fs=1.0, nperseg=min(256, len(x)//4))
        ax[1].pcolormesh(tt/60, f, 10*np.log10(Sxx+1e-20), shading="gouraud", cmap="viridis")
        ax[1].axhline(0.001, color="w", ls=":", lw=0.6); ax[1].axhline(0.01, color="w", ls=":", lw=0.6)
        ax[1].set_ylim(0, 0.02); ax[1].set_ylabel("freq (Hz)\nspectrogram")
        if HAVE_PYWT:
            freqs = np.logspace(np.log10(1e-4), np.log10(2e-2), 80)
            scales = pywt.central_frequency("morl") / freqs
            coef, _ = pywt.cwt(x, scales, "morl", sampling_period=1.0)
            ax[2].pcolormesh(t/60, freqs, np.abs(coef)**2, shading="gouraud", cmap="magma")
            ax[2].axvline(tpk/60, color="w", ls="--", lw=0.8)
            ax[2].set_ylim(1e-3, 1e-2); ax[2].set_ylabel("freq (Hz)\nMorlet CWT")
        ax[-1].set_xlabel("minutes from segment start")
        fig.tight_layout()
        fig.savefig(out / f"{p.dataset}_{p.station}_{p.component}_{p.sim}_{p.template}"
                    f"_{p.peak_time:%Y%m%d_%H%M}.png", dpi=120)
        plt.close(fig); n += 1
    print(f"  candidate panels: {n} → {out}")


# ── manifest ──────────────────────────────────────────────────────────────────
def write_manifest(peaks):
    sig = int(peaks.significant.sum())
    txt = f"""# Pipeline output manifest

Single entry point: `run_pipeline.py`. Re-run with `--plots-only` to rebuild figures
from existing data.

## Data products
- `clean_bandpassed/<dataset>/<station>_<dir|mag>_..._clean_bp.feather` (+ `.meta.json`)
- `SWCC_comprehensive/<dataset>/<station>_<comp>_peaks.csv`
- `SWCC_comprehensive/all_peaks.csv`, `all_peaks_flagged.csv` (significance + in_edge flags)
- `SWCC_comprehensive/summary_by_station_component.csv`, `significance_summary.txt`

## Figures
- `SWCC_comprehensive/<dataset>/<station>/<sim>/<station>_<comp>_<sim>_<template>_swcc.png`
      old-style per-template SWCC: blue |r| line, red shaded volcanic periods (INGV),
      0.2/0.5/0.7 + null-floor threshold lines, cyan-circle / red-star peak markers
- `SWCC_comprehensive/<dataset>/<station>/<station>_<comp>_swcc_overview.png`
      compact per-(station,component) overview: peaks vs time + null floor + cumulative
- `SWCC_comprehensive/dataset_comparison.png`        INGV vs experiment (|r|, SNR)
- `SWCC_comprehensive/component_comparison.png`      directional vs magnitude
- `SWCC_comprehensive/candidates/*.png`              waveform + spectrogram + Morlet CWT
      per significant detection (n={sig})
- `SWCC_comprehensive/characterization/`             timeline, synchronous events, station detections
- `SWCC_comprehensive/credibility/`                  filter response, before/after de-ringing, null test
- `pipeline_improvements/`                           before/after gallery of every processing change

## Old outputs intentionally dropped (obsolete under the new design)
- P-wave contamination-probability / impact plots, contaminated-vs-clean peak lists,
  station P-wave impact summaries. Reason: P-waves are now excised from the raw signal
  *before* analysis, so there is no contaminated/clean split to compare. The modern
  significance criterion is the phase-randomised **null floor** (`credibility/null_test.*`,
  `significance_summary.txt`), not the old r=0.2 threshold.

## Headline result
{sig} of {len(peaks)} peaks survive the null floor (~{100*sig/len(peaks):.1f}%); detections are
sparse and sit just above the floor (see PIPELINE_AUDIT.md §5 for the full, corrected discussion).
"""
    (BASE / "PIPELINE_MANIFEST.md").write_text(txt)
    print(f"  manifest → {BASE/'PIPELINE_MANIFEST.md'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plots-only", action="store_true",
                    help="skip recompute; rebuild figures from existing data")
    args = ap.parse_args()

    if not args.plots_only:
        run_compute()

    print("\n" + "="*70 + "\n[6-8] standard figures\n" + "="*70)
    peaks = pd.read_csv(SWCC / "all_peaks_flagged.csv", parse_dates=["peak_time"])
    volc = load_volc()
    print("  old-style per-template SWCC figures (blue line, red volcanic shading)…")
    swcc_oldstyle_plots.main()        # SWCC_comprehensive/<ds>/<st>/<sim>/..._swcc.png
    swcc_overview_plots(peaks, volc)  # compact per-(station,component) overview + cumulative
    dataset_comparison(peaks)
    template_sim_bars(peaks)
    candidate_panels(peaks)
    write_manifest(peaks)
    print("\n✅ pipeline complete.")


if __name__ == "__main__":
    main()
