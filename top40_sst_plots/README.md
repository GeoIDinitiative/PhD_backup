# Top 40 Station-Sim-Template Combination Plots

## 📊 Overview

These bar plots show the top 40 station-simulation-template combinations based on the number of matched peaks across different correlation thresholds, replicating the analysis from page 145 of the PhD thesis.

## 📄 Generated Plots

### 1. **01_top40_all_thresholds.png**
- **Correlation range:** r ≥ 0.2 (All thresholds)
- **Total unique SST combinations:** 112
- **Total peaks (top 40):** 49,839
- **Top performer:** ECPN-sim1-template3 (2,361 peaks)

### 2. **02_top40_moderate_02_05.png**
- **Correlation range:** 0.2 ≤ r < 0.5 (Moderate correlation)
- **Total unique SST combinations:** 112
- **Total peaks (top 40):** 48,711
- **Top performer:** ECPN-sim1-template3 (2,304 peaks)

### 3. **03_top40_high_quality_ge_05.png**
- **Correlation range:** r ≥ 0.5 (High quality)
- **Total unique SST combinations:** 75
- **Total peaks (top 40):** 1,210
- **Top performer:** ECPN-sim4-template1 (130 peaks)

## 📈 Key Statistics

### Overall Dataset (r ≥ 0.2):
- **Total peaks:** 61,554
- **INGV peaks:** 47,599 (77.3%)
- **IMPROVE peaks:** 13,955 (22.7%)
- **Unique SST combinations:** 112

### Moderate Correlation (0.2 ≤ r < 0.5):
- **Total peaks:** 60,265 (97.9% of all peaks)
- **Unique SST combinations:** 112

### High Correlation (r ≥ 0.5):
- **Total peaks:** 1,289 (2.1% of all peaks)
- **Unique SST combinations:** 75

## 🎨 Color Coding

- **Blue bars:** INGV dataset (eruptive period)
- **Red bars:** IMPROVE dataset (quiescence period)

## 🔍 Key Insights

### Top Performers (All Thresholds):
1. **ECPN-sim1-template3:** 2,361 peaks (INGV)
2. **ECPN-sim1-template1:** ~2,000+ peaks (INGV)
3. INGV station ECPN dominates the top positions

### Top Performers (High Quality r ≥ 0.5):
1. **ECPN-sim4-template1:** 130 peaks (INGV)
2. High-quality detections heavily favor INGV eruptive period
3. Template 1 and Template 3 perform best for high correlations

### Dataset Distribution:
- **INGV dominance:** 77.3% of all peaks (eruptive period)
- **IMPROVE contribution:** 22.7% of peaks (quiescence period)
- INGV's eruptive period shows much higher signal activity

### Template Performance:
- **Template 3:** Leads overall peak counts
- **Template 1:** Dominates high-quality (≥0.5) detections
- **Template 4:** Consistently weakest across all combinations

### Station Performance:
- **ECPN:** Primary detection station (INGV)
- **EC1:** Secondary strong performer (both datasets)
- IMPROVE stations show more distributed performance

## 📁 Data Files

Accompanying CSV files with detailed statistics:
- `top40_all_thresholds_data.csv`
- `top40_moderate_02_05_data.csv`
- `top40_high_quality_ge_05_data.csv`

Each file contains:
- Station-sim-template labels
- Peak counts
- Average correlation
- Average SNR (linear)

## 🔬 Methodology

### Data Source:
- Clean peaks from `/home/owen/Etna_signals/screening/SWCC_utc_fixed/`
- Both INGV and IMPROVE (experiment) datasets
- All stations, simulations (1-4), and templates (1-4)

### Processing:
1. Load all station-sim-template peak files
2. Filter by correlation threshold
3. Count peaks per SST combination
4. Rank by peak count
5. Display top 40 combinations

### Plot Features:
- Bars colored by dataset (INGV=blue, IMPROVE=red)
- X-axis labels show SST combinations
- Y-axis shows peak counts
- Statistics box shows totals and percentages
- High-resolution output (600 DPI)

## 🎯 Usage

These plots are designed for:
- PhD thesis integration (page 145 reference)
- Identifying best-performing station-sim-template combinations
- Comparing eruptive vs quiescence detection performance
- Understanding correlation distribution across combinations
- Selecting optimal templates for specific stations/periods

## 🛠️ Generation

Created using: `create_top40_sst_plots.py`

To regenerate:
```bash
python3 create_top40_sst_plots.py
```

## 📌 Notes

- Plots show clean peaks only (P-wave contaminated peaks excluded)
- Correlation thresholds based on PhD thesis methodology
- Top 40 limit matches original thesis analysis
- Some combinations may have < 40 entries if fewer combinations exist in that range
