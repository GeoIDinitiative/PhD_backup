# Tilt Denoising / Bandpass / SWCC — Pipeline Audit & Rebuild Record

_Audit date: 2026-05-29. Scope: observed-signal denoising, 0.001–0.01 Hz bandpass, and sliding-window
cross-correlation (SWCC) against simulation tilt templates, for the INGV (winter 2022–23) and
experiment (summer 2023) periods._

---

## 1. Data lineage

| Stage | Artifact | Notes |
|-------|----------|-------|
| Raw | `/home/owen/Signals/...` feathers + `EC1.csv` | 1 Hz nominal; INGV has east/north/mag, experiment has x/y(/mag) |
| P-wave excision | `tilt_raw_pwave_removed/<ST>_*.feather` | UTC datetimes; earthquake windows removed via `pwave_buffer.calculate_pwave_buffers` |
| **(rebuild input)** | ↑ these feathers | gaps = recording breaks **and** excised quake windows |
| Old bandpass | `bandpassed_exports_001_01_csv_only_utc/` | produced by `plot_fix/signals.py`; **kept for regression, not reused** |
| Templates | `tilt_templates_csv/<dataset>/` | 96 per dataset = 6 stations × 4 sims × 4 templates, single `x` component, ~3334 s long |

## 2. Findings

**A1 — Bandpass design (KEEP).** Butterworth order 4, `sos`, zero-phase `sosfiltfilt`, band
0.001–0.01 Hz, fs 1 Hz. Zero-phase is correct for preserving cross-correlation timing. Band verified
appropriate: periods 100–1000 s, above tidal/thermal noise, below the teleseismic body-wave band.

**A2 — Atmospheric/environmental cleaning ABSENT.** No barometric/thermal/tidal correction; the
declared `fc_highpass = 1e-5` is never applied. No barometer data exists on disk, so a pressure-
admittance correction is impossible. Physics mitigates this — environmental tilt noise is largely
at periods ≫1000 s, below the 0.001 Hz corner, so the bandpass removes most of it. Remaining
defects: **no detrend before filtering** (DC/drift → edge transients) and **no despiking**.
→ Rebuild adds per-segment detrend + Hampel/MAD despike. (Decision: spectral+detrend only.)

**A3 — Non-uniform time base filtered as uniform (DEFECT).** `seconds` is 1 Hz for >99.99 % of
samples but contains multi-million-second jumps at the recording breaks (~22 Nov, ~12 Dec).
`sosfiltfilt` across a gap rings for thousands of seconds. → Rebuild filters **per continuous
segment**.

**A4 — SWCC performance (BLOCKER).** Correct Pearson-r math, but a pure-Python loop:
≈2.8×10¹⁰ ops/template, run twice × 96 templates. → Rebuild vectorizes via cumsum running-moments +
FFT convolution (validated against the loop).

**A5 — Splice / through-quake artifacts (FLAW).** Old "cleaned" path concatenates non-contiguous
samples and correlates across splices (spurious peaks, mislabeled times); old "full" path filters
through earthquakes (ringing smears into clean windows). → Rebuild correlates **within segments
only**, on already-excised raw.

**A6 — Component mismatch.** Observed bandpass used `mag` (rectified ≥0); templates are bipolar
`x`. → Rebuild runs **both** a directional component (INGV `east`, experiment `x`) and magnitude,
and compares.

**A7 — UTC.** Old exports embed original time base + −1 h read-time anchor. Rebuild carries true UTC
datetimes end-to-end (input feathers already UTC); anchor bookkeeping retired.

**A8 — SNR & peaks (KEEP, document).** SNR guard 3000 s / half-window 10000 s; `find_peaks(|r|,
height=0.2, distance=1000)`. `|r|` counts anti-correlations as detections (polarity ambiguity) —
intended, now stated.

## 3. Rebuilt pipeline (this folder only; originals untouched)

1. `build_clean_bandpassed.py` — segment at gaps → despike (MAD) → linear detrend → per-segment
   zero-phase Butterworth bandpass 0.001–0.01 Hz → reassemble on UTC axis (gaps preserved).
   Output `clean_bandpassed/<dataset>/<ST>_<component>_0p001-0p01Hz_clean_bp.feather` + meta.
2. `swcc_comprehensive.py` — vectorized, segment-aware normalized SWCC for both components ×
   sims × templates → `SWCC_comprehensive/` peak CSVs, plots, mag-vs-directional comparison.
3. Credibility gates — numerical-equivalence vs old loop; filter magnitude/step response; old-vs-new
   before/after; phase-randomized null test for the r=0.2 false-peak floor.

## 4. Component map

| Station | Dataset | Directional | Magnitude |
|---------|---------|-------------|-----------|
| ECPN, EEC1 | ingv | `east` | `mag` |
| EC10, ECIT, ECOR, EMAS | experiment | `x` | `mag` |
| EC1 | experiment | `x` | √(x²+y²) (computed) |
