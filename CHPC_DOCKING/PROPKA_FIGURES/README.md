# BRD4 protonation audit figures

PROPKA 3.5.1 predictions for the ligand-free 2YEL chain A (BD1) and 2YEM chain B (BD2) structures. Figure 1 shows single-site Henderson–Hasselbalch protonated fractions estimated from the predicted histidine pKa values at pH 7.0, 7.4 and 7.5. Figure 2 reports relative folding free-energy estimates from the PROPKA summaries over pH 0–14. These are model predictions, not experimental folding or binding measurements.

Reproduce both figures (PNG at 600 DPI, PDF, SVG) and both CSV files:

```bash
python -m pip install numpy matplotlib==3.10.8
python plot_propka.py
```

Figure size: 130 × 100 mm; white background; bold DejaVu Sans labels and legends; no grid or in-plot titles.
