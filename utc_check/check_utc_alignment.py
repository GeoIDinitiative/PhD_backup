"""
UTC Alignment Diagnostic
========================
Visualises whether P-wave exclusion windows land on actual tilt disturbances.

Two outputs per station
  1. OVERVIEW  – full time series with a vertical line at every p_wave_eta
                 (magnitude >= MAG_THRESHOLD_OVERVIEW only, to avoid clutter)
  2. EVENTS    – individual zoom plots for the N largest earthquakes visible in
                 the signal window, each showing the tilt ± CONTEXT_HOURS around
                 p_wave_eta with the exclusion zone shaded

Outputs go to:  utc_check/outputs/<station>/
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
BASE          = Path("/home/owen/etna_signals_phd")
BP_ROOT       = BASE / "bandpassed_exports_001_01_csv_only_utc"
EQ_CSV        = BASE / "earthquakes_merged_utc.csv"
OUT_ROOT      = Path(__file__).parent / "outputs"

# ── signal epoch anchors (same as swcc_edit.py) ────────────────────────────────
START_DT = {
    "ingv":       pd.Timestamp("2022-11-14 22:00:00"),   # was 23:00; -1h: both datasets clocked on CET (UTC+1), not CEST
    "experiment": pd.Timestamp("2023-07-23 23:00:00"),   # was 2023-07-24 00:00; -1h: same reason
}

# ── stations to check ──────────────────────────────────────────────────────────
STATIONS = {
    "ingv":       ["ECPN", "EEC1"],       # EEC1 file, labelled EC1 in analysis
    "experiment": ["EC1", "EC10", "ECIT", "ECOR", "EMAS"],
}

# ── exclusion buffer logic (mirrors swcc_edit.py) ──────────────────────────────
def buffer_minutes(magnitude):
    if pd.isna(magnitude):
        return 10.0
    elif magnitude >= 5.0:
        return 15.0
    elif magnitude >= 4.0:
        return 10.0
    else:
        return 7.0

# ── plot parameters ────────────────────────────────────────────────────────────
MAG_THRESHOLD_OVERVIEW = 5.5   # only draw lines for M >= this on overview
N_EVENT_PLOTS          = 60    # number of individual event zoom plots per station
CONTEXT_HOURS          = 3.0   # ± hours shown around p_wave_eta in event plots
DPI                    = 120


# ── helpers ────────────────────────────────────────────────────────────────────

def load_signal(dataset, station):
    """Return (time_dt array, bandpassed array) using the same epoch as swcc_edit.py."""
    csv_path = BP_ROOT / dataset / f"{station}_0p001-0p01Hz_bp.csv"
    df = pd.read_csv(csv_path)
    time_dt = START_DT[dataset] + pd.to_timedelta(df["time_seconds"].values, unit="s")
    return time_dt, df["bandpassed"].values


def filter_eq_to_window(eq_df, t_start, t_end):
    """Keep only earthquakes whose p_wave_eta falls inside [t_start, t_end]."""
    eta = pd.to_datetime(eq_df["p_wave_eta"])
    mask = (eta >= t_start) & (eta <= t_end)
    out = eq_df[mask].copy()
    out["p_wave_eta_dt"] = eta[mask].values
    return out.sort_values("magnitude", ascending=False).reset_index(drop=True)


def signal_window(time_dt, signal, centre, half_hours):
    """Slice signal to centre ± half_hours."""
    delta = pd.Timedelta(hours=half_hours)
    mask = (time_dt >= centre - delta) & (time_dt <= centre + delta)
    return time_dt[mask], signal[mask]


def apply_full_datetime_ticks(ax, span_hours):
    """Format x-axis with full YYYY-MM-DD HH:MM labels, density scaled to window width."""
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
    ax.xaxis.set_minor_locator(mdates.AutoDateLocator())
    ax.tick_params(axis="x", which="major", labelsize=7.5)
    plt.setp(ax.xaxis.get_majorticklabels(), ha="center")
    # Force matplotlib to treat the axis as dates
    ax.figure.autofmt_xdate(rotation=0, ha="center")


# ── plot 1: overview ───────────────────────────────────────────────────────────

def plot_overview(time_dt, signal, eq_in_window, station, dataset, out_dir):
    big_eq = eq_in_window[eq_in_window["magnitude"] >= MAG_THRESHOLD_OVERVIEW]

    fig, ax = plt.subplots(figsize=(20, 4))
    ax.plot(pd.DatetimeIndex(time_dt), signal, "-", lw=0.4, color="steelblue", rasterized=True)

    for _, row in big_eq.iterrows():
        eta = row["p_wave_eta_dt"]
        mag = row["magnitude"]
        buf = pd.Timedelta(minutes=buffer_minutes(mag))
        ax.axvline(eta, color="red",    lw=0.6, alpha=0.7, zorder=3)
        ax.axvspan(eta - buf, eta + buf, color="red", alpha=0.07, zorder=2)

    ax.set_title(
        f"{dataset.upper()} – {station}  |  full time series\n"
        f"Red lines = p_wave_eta, shading = exclusion zone  "
        f"(showing M ≥ {MAG_THRESHOLD_OVERVIEW}, n={len(big_eq)})",
        fontsize=9
    )
    ax.set_xlabel("UTC datetime")
    ax.set_ylabel("Tilt (bandpassed)")
    ax.margins(x=0.002)

    # legend
    ax.plot([], [], color="red", lw=1.2, label=f"p_wave_eta  (M≥{MAG_THRESHOLD_OVERVIEW})")
    ax.fill_between([], [], [], color="red", alpha=0.2, label="exclusion zone")
    ax.legend(fontsize=7, loc="upper right")

    span_h = (pd.Timestamp(time_dt[-1]) - pd.Timestamp(time_dt[0])).total_seconds() / 3600
    apply_full_datetime_ticks(ax, span_h)

    fig.tight_layout()
    fpath = out_dir / f"{station}_overview.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)
    print(f"  saved: {fpath.name}")


# ── plot 2: individual event zoom ──────────────────────────────────────────────

def plot_event_zoom(time_dt, signal, row, rank, station, dataset, out_dir, eq_df_all):
    eta   = row["p_wave_eta_dt"]
    mag   = row["magnitude"]
    dist  = row.get("distance_km", float("nan"))
    buf   = pd.Timedelta(minutes=buffer_minutes(mag))

    eta_plus1h = eta + pd.Timedelta(hours=1)
    eta_plus2h = eta + pd.Timedelta(hours=2)

    t_win, s_win = signal_window(time_dt, signal, eta, CONTEXT_HOURS)
    if len(t_win) == 0:
        return

    # other catalog arrivals inside the plot window (excluding the primary event)
    win_start = eta - pd.Timedelta(hours=CONTEXT_HOURS)
    win_end   = eta + pd.Timedelta(hours=CONTEXT_HOURS)
    others = eq_df_all[
        (eq_df_all["p_wave_eta"] >= win_start) &
        (eq_df_all["p_wave_eta"] <= win_end) &
        (eq_df_all["magnitude"] >= 3.0) &
        (abs(eq_df_all["p_wave_eta"] - eta) > pd.Timedelta(seconds=10))
    ]

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(pd.DatetimeIndex(t_win), s_win, "-", lw=0.8, color="steelblue", zorder=4)

    # other P-wave arrivals — dashed pink
    first_other = True
    for _, orow in others.iterrows():
        o_eta = orow["p_wave_eta"]
        o_mag = orow["magnitude"]
        lbl = f"other p_wave_eta in window (n={len(others)})" if first_other else "_"
        ax.axvline(o_eta, color="deeppink", lw=1.0, ls="--", alpha=0.7, zorder=3, label=lbl)
        first_other = False

    # primary ETA + exclusion zone
    ax.axvline(eta, color="red", lw=1.8, zorder=6, label=f"p_wave_eta  M{mag:.1f} (primary)")
    ax.axvspan(eta - buf, eta + buf,
               color="red", alpha=0.10, zorder=2,
               label=f"exclusion ±{buffer_minutes(mag):.0f} min")

    # hypothesis lines — solid, thicker
    ax.axvline(eta_plus1h, color="darkorange", lw=1.8, ls="-", zorder=5,
               label="ETA + 1 h  (signal 1 h late)")
    ax.axvline(eta_plus2h, color="purple",     lw=1.8, ls="-", zorder=5,
               label="ETA + 2 h  (signal 2 h late)")

    ax.set_title(
        f"{dataset.upper()} {station}  |  "
        f"M{mag:.1f}  dist={dist:.0f} km  "
        f"p_wave_eta = {eta.strftime('%Y-%m-%d %H:%M:%S')} UTC\n"
        f"Window = ETA ± {CONTEXT_HOURS:.0f} h    buffer = ±{buffer_minutes(mag):.0f} min    "
        f"Pink dashed = other catalog arrivals in window",
        fontsize=8.5
    )
    ax.set_xlabel("UTC datetime")
    ax.set_ylabel("Tilt (bandpassed)")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.margins(x=0.01)

    apply_full_datetime_ticks(ax, CONTEXT_HOURS * 2)

    fig.tight_layout()
    fpath = out_dir / f"{station}_event_{rank:03d}_M{mag:.1f}_{eta.strftime('%Y%m%d_%H%M')}.png"
    fig.savefig(fpath, dpi=DPI)
    plt.close(fig)


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    eq_df = pd.read_csv(EQ_CSV)
    eq_df["p_wave_eta"] = pd.to_datetime(eq_df["p_wave_eta"])

    for dataset, stations in STATIONS.items():
        for station in stations:
            print(f"\n{'='*60}")
            print(f"  {dataset.upper()}  {station}")

            out_dir = OUT_ROOT / station
            out_dir.mkdir(parents=True, exist_ok=True)

            try:
                time_dt, signal = load_signal(dataset, station)
            except FileNotFoundError as e:
                print(f"  SKIP – {e}")
                continue

            t_start = pd.Timestamp(time_dt[0])
            t_end   = pd.Timestamp(time_dt[-1])
            print(f"  Signal window: {t_start}  →  {t_end}")

            eq_win = filter_eq_to_window(eq_df, t_start, t_end)
            print(f"  Earthquakes in window: {len(eq_win)}  "
                  f"(M≥{MAG_THRESHOLD_OVERVIEW}: {(eq_win['magnitude']>=MAG_THRESHOLD_OVERVIEW).sum()})")

            if len(eq_win) == 0:
                print("  No earthquakes in signal window – skipping")
                continue

            # ── overview ──
            print("  Generating overview …")
            plot_overview(time_dt, signal, eq_win, station, dataset, out_dir)

            # ── individual event zooms ──
            top_n = eq_win.head(N_EVENT_PLOTS)
            print(f"  Generating {len(top_n)} event zoom plots …")
            for rank, (_, row) in enumerate(top_n.iterrows(), start=1):
                plot_event_zoom(time_dt, signal, row, rank, station, dataset, out_dir, eq_df)

            print(f"  Done – {station}")

    print(f"\n{'='*60}")
    print(f"All outputs in: {OUT_ROOT}")


if __name__ == "__main__":
    main()
