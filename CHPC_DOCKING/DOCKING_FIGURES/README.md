# BRD4 docking QC and grid sensitivity figures

The executed notebook contains 11 figures generated in a Python kernel from the CHPC terminal summaries supplied on 26 November 2026. All 12 code cells executed without errors. Source values are embedded in the notebook. Run all cells to regenerate SVG, 600-dpi PNG, and numerical CSV exports in `figures/`.

## Provenance
- Primary docking: 15 ligands × two domains; expanded grid (+4 Å per dimension): another 30 calculations.
- Expanded BD2 scores reconstructed to three decimal places from the reported BD1 expanded scores and BD1−BD2 differences; Spearman correlations and MAEs independently reproduced from these values.
- Source is terminal summaries, not the original coordinate files or raw logs; coordinate-based analysis has not been rerun here. Reported structural QC does not constitute experimental validation. Original grids remain the primary protocol.
- Contact distances are geometric, not hydrogen-bond assignments. Single-seed sensitivity analysis does not establish multi-seed robustness.

## Run
Install `numpy pandas matplotlib scipy nbformat nbclient ipykernel` and execute the notebook with a Python 3 kernel. Existing outputs are retained in the committed notebook.
