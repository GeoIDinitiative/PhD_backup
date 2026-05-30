"""
UTC Alignment Diagnostic — RAW tilt signal
===========================================
Same layout as check_utc_alignment.py but using the unfiltered feather/CSV
tilt data instead of the bandpassed denoised signal.

Signal used:
  INGV      — 'east' component (µrad), from INGV feather files
  Experiment — 'mag' column for feather stations, 'x' column for EC1 CSV

Datetimes come directly from the 'datetime' column in each file (already
consistent with the corrected UTC epoch established in check_utc_alignment.py).

Outputs go to:  utc_check/outputs_raw/<station>/
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path

# pwave_buffer.py lives one directory up (etna_signals_phd/)
sys.path.insert(0, str(Path(__file__).parent.parent))
from pwave_buffer import calculate_pwave_buffers

# ── paths ──────────────────────────────────────────────────────────────────────
EQ_CSV   = Path("/home/owen/etna_signals_phd/earthquakes_merged_utc.csv")
OUT_ROOT = Path(__file__).parent / "outputs_raw"

RAW_SOURCES = {
    "ingv": {
        "ECPN": (Path("/home/owen/Signals/experiment/INGV/ECPN.feather"), "feather", "east"),
        "EEC1": (Path("/home/owen/Signals/experiment/INGV/EEC1.feather"),  "feather", "east"),
    },
    "experiment": {
        "EC1":  (Path("/home/owen/Signals/experiment/EC1.csv"),                                    "csv",     " x"),
        "EC10": (Path("/home/owen/Signals/experiment/school-data/INGV_feather/EC10.feather"),      "feather", "mag"),
        "ECIT": (Path("/home/owen/Signals/experiment/school-data/INGV_feather/ECIT.feather"),      "feather", "mag"),
        "ECOR": (Path("/home/owen/Signals/experiment/school-data/INGV_feather/ECOR.feather"),      "feather", "mag"),
        "EMAS": (Path("/home/owen/Signals/experiment/school-data/INGV_feather/EMAS.feather"),      "feather", "mag"),
    },
}

# ── plot parameters ────────────────────────────────────────────────────────────
MAG_THRESHOLD_OVERVIEW = 5.5
N_EVENT_PLOTS          = 60
CONTEXT_HOURS          = 3.0
DPI                    = 120


# buffer logic delegated to pwave_buffer.py (asymmetric pre/post)


# ── data loading ───────────────────────────────────────────────────────────────

def load_raw_signal(dataset, station):
    path, fmt, sig_col = RAW_SOURCES[dataset][station]
    if fmt == "feather":
        df = pd.read_feather(path)
    else:
        df = pd.read_csv(path)
    df["datetime"] = pd.to_datetime(df["datetime"])
    df = df.sort_values("datetime").reset_index(drop=True)
    return df["datetime"].values, df[sig_col].values


# ── helpers ────────────────────────────────────────────────────────────────────

def filter_eq_to_window(eq_df, t_start, t_end):
    eta = pd.to_datetime(eq_df["p_wave_eta"])
    mask = (eta >= t_start) & (eta <= t_end)
    out = eq_df[mask].copy()
    out["p_wave_eta_dt"] = eta[mask].values
    return out.sort_values("magnitude", ascending=False).reset_index(drop=True)


def signal_window(time_dt, signal, centre, half_hours):
    delta = pd.Timedelta(hours=half_hours)
    mask = (pd.DatetimeIndex(time_dt) >= centre - delta) & \
           (pd.DatetimeIndex(time_dt) <= centre + delta)
    return time_dt[mask], signal[mask]


def apply_full_datetime_ticks(ax, span_hours):
    if span_hours <= 2:
        locator = mdates.MinuteLocator(byminute=[0, 15, 30, 45])
        fmt     = "%Y-%m-%d\n%H:%M"
    elif span_hours <= 8:
        locator = mdates.HourLocator(interval=1)
        fmt     = "%Y-%m-%d\n%H:%M"
    elif span_hours <= 24:
        locator = mdates.HourLocator(interval=2)
        fmt     = "%Y-%m-%d\n%H:%M"
    elif span_hours <= 72:
        locator = mdates.HourLocator(interval=6)
        fmt     = "%Y-%m-%d\n%H:%M"
    else:
        locator = mdates.DayLocator(interval=7)
        fmt     = "%Y-%m-%d"
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    ax.tick_params(axis="x", which="major", labelsize=7.5)
    plt.setp(ax.xaxis.get_majorticklabels(), ha="center")
    ax.figure.autofmt_xdate(rotation=0, ha="center")


# ── plot 1: overview ───────────────────────────────────────────────────────────

def plot_overview(time_dt, signal, eq_in_window, station, dataset, sig_col, out_dir):
    big_eq = eq_in_window[eq_in_window["magnitude"] >= MAG_THRESHOLD_OVERVIEW]

    fig, ax = plt.subplots(figsize=(20, 4))
    ax.plot(pd.DatetimeIndex(time_dt), signal, lw=0.3, color="steelblue", rasterized=True)

    for _, row in big_eq.iterrows():
        eta  = row["p_wave_eta_dt"]
        mag  = row["magnitude"]
        dist = row.get("distance_km", None)
        pre_min, post_min = calculate_pwave_buffers(dist, mag)
        ax.axvline(eta, color="red", lw=0.6, alpha=0.7, zorder=3)
        ax.axvspan(eta - pd.Timedelta(minutes=pre_min),
                   eta + pd.Timedelta(minutes=post_min),
                   color="red", alpha=0.07, zorder=2)

    ax.set_title(
        f"{dataset.upper()} – {station}  |  raw tilt ({sig_col})  |  full time series\n"
        f"Red lines = p_wave_eta, shading = exclusion zone  "
        f"(showing M ≥ {MAG_THRESHOLD_OVERVIEW}, n={len(big_eq)})",
        fontsize=9
    )
    ax.set_xlabel("UTC datetime")
    ax.set_ylabel(f"Tilt raw ({sig_col})  [µrad]")
    ax.margins(x=0.002)
    ax.plot([], [], color="red", lw=1.2, label=f"p_wave_eta  (M≥{MAG_THRESHOLD_OVERVIEW})")
    ax.fill_between([], [], [], color="red", alpha=0.2, label="exclusion zone")
    ax.legend(fontsize=7, loc="upper right")

    span_h = (pd.Timestamp(time_dt[-1]) - pd.Timestamp(time_dt[0])).total_seconds() / 3600
    apply_full_datetime_ticks(ax, span_h)

    fig.tight_layout()
    fpath = out_dir / f"{station}_overview_raw.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)
    print(f"  saved: {fpath.name}")


# ── plot 2: event zoom ─────────────────────────────────────────────────────────

def plot_event_zoom(time_dt, signal, row, rank, station, dataset, sig_col, out_dir, eq_df_all):
    eta  = row["p_wave_eta_dt"]
    mag  = row["magnitude"]
    dist = row.get("distance_km", None)
    pre_min, post_min = calculate_pwave_buffers(dist, mag)

    eta_plus1h = eta + pd.Timedelta(hours=1)
    eta_plus2h = eta + pd.Timedelta(hours=2)

    t_win, s_win = signal_window(time_dt, signal, eta, CONTEXT_HOURS)
    if len(t_win) == 0:
        return

    # detrend within window (remove linear trend so seismic signal is visible)
    if len(s_win) > 2:
        coeffs = np.polyfit(np.arange(len(s_win)), s_win, 1)
        s_win  = s_win - np.polyval(coeffs, np.arange(len(s_win)))

    # other catalog arrivals M>=3 in window
    win_start = eta - pd.Timedelta(hours=CONTEXT_HOURS)
    win_end   = eta + pd.Timedelta(hours=CONTEXT_HOURS)
    others = eq_df_all[
        (eq_df_all["p_wave_eta"] >= win_start) &
        (eq_df_all["p_wave_eta"] <= win_end) &
        (eq_df_all["magnitude"] >= 3.0) &
        (abs(eq_df_all["p_wave_eta"] - eta) > pd.Timedelta(seconds=10))
    ]

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(pd.DatetimeIndex(t_win), s_win, "-", lw=0.7, color="steelblue", zorder=4)

    # other arrivals
    first_other = True
    for _, orow in others.iterrows():
        lbl = f"other p_wave_eta in window (n={len(others)}, M≥3)" if first_other else "_"
        ax.axvline(orow["p_wave_eta"], color="deeppink", lw=1.0, ls="--", alpha=0.7,
                   zorder=3, label=lbl)
        first_other = False

    # primary ETA + asymmetric exclusion zone
    ax.axvline(eta, color="red", lw=1.8, zorder=6, label=f"p_wave_eta  M{mag:.1f} (primary)")
    ax.axvspan(eta - pd.Timedelta(minutes=pre_min),
               eta + pd.Timedelta(minutes=post_min),
               color="red", alpha=0.10, zorder=2,
               label=f"exclusion  -{pre_min:.0f} / +{post_min:.0f} min")

    # hypothesis lines
    ax.axvline(eta_plus1h, color="darkorange", lw=1.8, ls="-", zorder=5,
               label="ETA + 1 h  (signal 1 h late)")
    ax.axvline(eta_plus2h, color="purple",     lw=1.8, ls="-", zorder=5,
               label="ETA + 2 h  (signal 2 h late)")

    ax.set_title(
        f"{dataset.upper()} {station}  |  raw tilt ({sig_col}, detrended in window)  |  "
        f"M{mag:.1f}  dist={dist:.0f} km  p_wave_eta = {eta.strftime('%Y-%m-%d %H:%M:%S')} UTC\n"
        f"Window = ETA ± {CONTEXT_HOURS:.0f} h    exclusion: -{pre_min:.0f} / +{post_min:.0f} min    "
        f"Pink dashed = other catalog arrivals (M≥3) in window",
        fontsize=8.5
    )
    ax.set_xlabel("UTC datetime")
    ax.set_ylabel(f"Tilt raw ({sig_col}, detrended)  [µrad]")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.margins(x=0.01)
    apply_full_datetime_ticks(ax, CONTEXT_HOURS * 2)

    fig.tight_layout()
    year_dir = out_dir / str(eta.year)
    year_dir.mkdir(exist_ok=True)
    fpath = year_dir / f"{station}_event_{rank:03d}_M{mag:.1f}_{eta.strftime('%Y%m%d_%H%M')}.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    eq_df = pd.read_csv(EQ_CSV)
    eq_df["p_wave_eta"] = pd.to_datetime(eq_df["p_wave_eta"])

    for dataset, stations in RAW_SOURCES.items():
        for station, (path, fmt, sig_col) in stations.items():
            print(f"\n{'='*60}")
            print(f"  {dataset.upper()}  {station}  ({sig_col})")

            out_dir = OUT_ROOT / station
            out_dir.mkdir(parents=True, exist_ok=True)

            try:
                time_dt, signal = load_raw_signal(dataset, station)
            except Exception as e:
                print(f"  SKIP – {e}")
                continue

            t_start = pd.Timestamp(time_dt[0])
            t_end   = pd.Timestamp(time_dt[-1])
            print(f"  Signal window: {t_start}  →  {t_end}  ({len(time_dt):,} samples)")

            eq_win = filter_eq_to_window(eq_df, t_start, t_end)
            print(f"  Earthquakes in window: {len(eq_win)}  "
                  f"(M≥{MAG_THRESHOLD_OVERVIEW}: {(eq_win['magnitude']>=MAG_THRESHOLD_OVERVIEW).sum()})")

            if len(eq_win) == 0:
                print("  No earthquakes in window – skipping")
                continue

            print("  Generating overview …")
            plot_overview(time_dt, signal, eq_win, station, dataset, sig_col, out_dir)

            top_n = eq_win.head(N_EVENT_PLOTS)
            print(f"  Generating {len(top_n)} event zoom plots …")
            for rank, (_, row) in enumerate(top_n.iterrows(), start=1):
                plot_event_zoom(time_dt, signal, row, rank, station, dataset,
                                sig_col, out_dir, eq_df)

            print(f"  Done – {station}")

    print(f"\n{'='*60}")
    print(f"All outputs in: {OUT_ROOT}")


if __name__ == "__main__":
    main()
