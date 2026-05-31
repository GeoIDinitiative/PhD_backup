"""
swcc_analysis_plots.py  —  stage 7: clean-peak analysis plots (Design B)
=======================================================================
Consumes the per-detection list written by swcc_continuous.py
(SWCC_comprehensive/continuous/all_detections_continuous.csv) and regenerates the
analysis plots that the SWCC per-template figures don't cover:

  · detection_overview_<dataset>.png   score vs time per station (MAX & STACK),
        detection/significance floors, significant points highlighted, volcanic
        overlay for INGV
  · detections_by_station.png          bar chart: detections vs significant per station/method
  · candidates/<...>.png               waveform + spectrogram + Morlet CWT for each
        significant detection

Output: SWCC_comprehensive/analysis/
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import spectrogram

from swcc_oldstyle_plots import load_volcanic_events, plot_volcanic_events_on_swcc

try:
    import pywt
    HAVE_PYWT = True
except Exception:
    HAVE_PYWT = False

warnings.filterwarnings("ignore")
BASE = Path("/home/owen/etna_signals_phd")
CONT = BASE / "continuous_bandpassed"
SWCC = BASE / "SWCC_comprehensive"
OUT  = SWCC / "analysis"
DETS = SWCC / "continuous" / "all_detections_continuous.csv"
VOLC = load_volcanic_events(BASE / "etna_volcanic_events_cleaned.csv")
PERIODS = {"ingv": ("2022-11-14", "2023-03-01"), "experiment": ("2023-07-23", "2023-08-03")}
COL = {"max": "#2563eb", "stack": "#f59e0b"}


def load_cont(ds, st, comp):
    f = CONT / ds / f"{st}_{comp}_0p001-0p01Hz_cont_bp.feather"
    if not f.exists():
        return None
    d = pd.read_feather(f); d["datetime"] = pd.to_datetime(d["datetime"])
    return d


# ── 1. per-dataset detection overview ─────────────────────────────────────────
def overview(dets):
    for ds, g0 in dets.groupby("dataset"):
        stations = sorted(g0.station.unique())
        fig, axes = plt.subplots(len(stations), 1, figsize=(14, 2.4*len(stations)),
                                 squeeze=False, sharex=True)
        for ax, st in zip(axes[:, 0], stations):
            gst = g0[g0.station == st]
            if ds == "ingv":
                plot_volcanic_events_on_swcc(ax, VOLC, (pd.Timestamp(PERIODS[ds][0]), pd.Timestamp(PERIODS[ds][1])))
            for method in ["max", "stack"]:
                gm = gst[gst.method == method]
                if gm.empty:
                    continue
                ax.scatter(gm.peak_time, gm.score, s=14, c=COL[method], alpha=0.5, label=f"{method} detect")
                sig = gm[gm.significant]
                if len(sig):
                    ax.scatter(sig.peak_time, sig.score, s=70, marker="*", c=COL[method],
                               edgecolors="k", linewidths=0.6, label=f"{method} significant")
                ax.axhline(gm.floor_detect.iloc[0], ls="--", c=COL[method], alpha=0.5, lw=1)
                ax.axhline(gm.floor_signif.iloc[0], ls="-", c=COL[method], alpha=0.7, lw=1.2)
            ax.set_ylabel(f"{st}\nscore"); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2, loc="upper right")
            ax.set_xlim(pd.Timestamp(PERIODS[ds][0]), pd.Timestamp(PERIODS[ds][1]))
        axes[0, 0].set_title(f"{ds}: detections over time (dashed=detect floor, solid=significance floor"
                             + ("; green=volcanic events)" if ds == "ingv" else ")"))
        fig.tight_layout(); fig.savefig(OUT / f"detection_overview_{ds}.png", dpi=140); plt.close(fig)
    print(f"  overviews → {OUT}/detection_overview_*.png")


# ── 2. detections-by-station bar chart ────────────────────────────────────────
def by_station(dets):
    g = (dets.groupby(["dataset", "station", "method"])
         .agg(detect=("significant", "size"), signif=("significant", "sum")).reset_index())
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, ds in zip(axes, ["ingv", "experiment"]):
        sub = g[g.dataset == ds]
        stations = sorted(sub.station.unique())
        x = np.arange(len(stations)); w = 0.2
        for j, (method, hatch) in enumerate([("max", None), ("stack", None)]):
            sm = sub[sub.method == method].set_index("station").reindex(stations).fillna(0)
            ax.bar(x + (j-0.5)*w*2, sm.detect, w*2, color=COL[method], alpha=0.4, label=f"{method} detect")
            ax.bar(x + (j-0.5)*w*2, sm.signif, w*2, color=COL[method], label=f"{method} significant")
        ax.set_xticks(x); ax.set_xticklabels(stations, rotation=30, ha="right")
        ax.set_ylabel("count"); ax.set_title(f"{ds}: detections vs significant by station")
        ax.legend(fontsize=8, loc="upper right"); ax.grid(axis="y", alpha=0.3)
        ax.set_ylim(top=max(1, sub.detect.max())*1.3)
    fig.tight_layout(); fig.savefig(OUT / "detections_by_station.png", dpi=140); plt.close(fig)
    print(f"  bar chart → {OUT}/detections_by_station.png")


# ── 3. candidate panels for significant detections ────────────────────────────
CAP_CANDIDATES = 20   # per dataset — the rest are chance-level (synchrony test is null)


def candidates(dets):
    cdir = OUT / "candidates"; cdir.mkdir(parents=True, exist_ok=True)
    for old in cdir.glob("*.png"):     # clear stale panels
        old.unlink()
    # cap to the top-scoring significant detections per dataset (avoid hundreds of chance panels)
    sig = (dets[dets.significant].sort_values("score", ascending=False)
           .groupby("dataset").head(CAP_CANDIDATES).copy())
    n = 0
    for _, p in sig.iterrows():
        d = load_cont(p.dataset, p.station, p.component)
        if d is None:
            continue
        t0 = pd.Timestamp(p.peak_time)
        seg = d[(d.datetime >= t0 - pd.Timedelta(minutes=60)) & (d.datetime <= t0 + pd.Timedelta(minutes=90))]
        if len(seg) < 256:
            continue
        x = seg.bandpassed.to_numpy()
        tt = (seg.datetime - seg.datetime.iloc[0]).dt.total_seconds().to_numpy() / 60
        tpk = (t0 - seg.datetime.iloc[0]).total_seconds() / 60
        nrow = 3 if HAVE_PYWT else 2
        fig, ax = plt.subplots(nrow, 1, figsize=(11, 3*nrow))
        ax[0].plot(tt, x, lw=0.7, color=COL[p.method]); ax[0].axvline(tpk, ls="--", c="k")
        ax[0].set_ylabel("bandpassed"); ax[0].grid(alpha=0.3)
        ax[0].set_title(f"{p.dataset}/{p.station}/{p.component} · {p.method} · score={p.score:.2f} · {t0:%Y-%m-%d %H:%M}")
        f, ts, Sxx = spectrogram(x, fs=1.0, nperseg=min(256, len(x)//4))
        ax[1].pcolormesh(ts/60, f, 10*np.log10(Sxx+1e-20), shading="gouraud", cmap="viridis")
        ax[1].axhline(0.001, color="w", ls=":", lw=0.6); ax[1].axhline(0.01, color="w", ls=":", lw=0.6)
        ax[1].set_ylim(0, 0.02); ax[1].set_ylabel("freq (Hz)\nspectrogram")
        if HAVE_PYWT:
            freqs = np.logspace(np.log10(1e-4), np.log10(2e-2), 80)
            sc = pywt.central_frequency("morl") / freqs
            coef, _ = pywt.cwt(x, sc, "morl", sampling_period=1.0)
            ax[2].pcolormesh(tt, freqs, np.abs(coef)**2, shading="gouraud", cmap="magma")
            ax[2].axvline(tpk, color="w", ls="--", lw=0.8); ax[2].set_ylim(1e-3, 1e-2); ax[2].set_ylabel("freq (Hz)\nMorlet CWT")
        ax[-1].set_xlabel("minutes from window start")
        fig.tight_layout()
        fig.savefig(cdir / f"{p.dataset}_{p.station}_{p.component}_{p.method}_{t0:%Y%m%d_%H%M}.png", dpi=120)
        plt.close(fig); n += 1
    print(f"  candidate panels: {n} → {cdir}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not DETS.exists():
        print(f"  {DETS} not found — run swcc_continuous.py (stage 3) first."); return
    dets = pd.read_csv(DETS, parse_dates=["peak_time"])
    if dets.empty:
        print("  no detections to plot."); return
    overview(dets)
    by_station(dets)
    candidates(dets)
    n_sig = int(dets.significant.sum())
    print(f"\nclean-peak analysis: {len(dets)} detections, {n_sig} significant → {OUT}")


if __name__ == "__main__":
    main()
