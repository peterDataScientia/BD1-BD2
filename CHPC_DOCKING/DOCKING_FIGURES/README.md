# BRD4 docking QC and grid sensitivity figures

Executed notebook: 12 code cells, 11 figures. Revised to match the approved Figentra v0.1.2 visual settings (https://github.com/peterDataScientia/figentra/blob/main/src/figentra/style.py): DejaVu Sans, 15 pt bold labels, 12 pt bold ticks, 11 pt bold frameless legends, 2 pt axes and ticks, 6 pt outward ticks, no minor ticks or grid, white backgrounds and no axes titles. The default canvas is 130 × 100 mm; dense and multi-panel plots are enlarged to preserve legibility. Run all cells to generate 600 dpi PNG, PDF and SVG exports; the 11 committed SVGs are the executed notebook outputs. Existing numerical CSV tables remain unchanged.

## Provenance
- CHPC terminal summaries supplied on 26 November 2026; 15 ligands, two domains, two grid protocols, 60 calculations.
- Expanded BD2 scores reconstructed to three decimal places from reported BD1 expanded scores and BD1−BD2 score differences. Reported summary statistics independently reproduced.
- Coordinates and raw docking logs were not rerun. Original grids remain primary. Geometric contacts do not establish hydrogen bonds. Single-seed comparison does not demonstrate multi-seed robustness.

## Run
Install `numpy pandas matplotlib scipy nbformat nbclient ipykernel` and execute all cells using Python 3. Style settings are embedded to make the notebook self-contained.
