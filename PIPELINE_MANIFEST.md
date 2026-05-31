# Pipeline output manifest

Single entry point: `run_pipeline.py`. Re-run with `--plots-only` to rebuild figures
from existing data.

## Data products
- `clean_bandpassed/<dataset>/<station>_<dir|mag>_..._clean_bp.feather` (+ `.meta.json`)
- `SWCC_comprehensive/<dataset>/<station>_<comp>_peaks.csv`
- `SWCC_comprehensive/all_peaks.csv`, `all_peaks_flagged.csv` (significance + in_edge flags)
- `SWCC_comprehensive/summary_by_station_component.csv`, `significance_summary.txt`

## Figures
- `SWCC_comprehensive/<dataset>/<station>/<station>_<comp>_swcc_overview.png`
      peak |r| vs time + threshold + null floor + cumulative (volcanic overlay for INGV)
- `SWCC_comprehensive/dataset_comparison.png`        INGV vs experiment (|r|, SNR)
- `SWCC_comprehensive/component_comparison.png`      directional vs magnitude
- `SWCC_comprehensive/candidates/*.png`              waveform + spectrogram + Morlet CWT
      per significant detection (n=37)
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
37 of 14976 peaks survive the null floor (~0.2%); detections are
sparse and sit just above the floor (see PIPELINE_AUDIT.md §5 for the full, corrected discussion).
