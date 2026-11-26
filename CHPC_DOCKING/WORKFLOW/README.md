# BRD4 docking workflow archive

This directory consolidates the commands actually shared during the CHPC receptor/protonation/QC stage into reusable scripts. Some wrappers and validation assertions are newly consolidated, not byte-identical copies of on-cluster scripts. The baseline PBS is reconstructed from workflow notes, not verified against the original. No new CHPC execution has been performed by this archive operation.

## Run order
From a CHPC login shell, with the pre-existing project directory ~/BRD4_NP_DOCKING:
1. On login1, source 01_environment.sh if creating the environment afresh (temporary compatibility pin for local preparation only).
2. bash 02_run_propka.sh
3. bash 03_initial_receptor_audit.sh
4. bash 04_prepare_ph75.sh
5. bash 05_pdb_qc.sh
6. ~/brd4_protonation_env/bin/python 06_heavy_atom_audit.py --stage pdb
7. ~/brd4_protonation_env/bin/python 07_sidechain_proximity.py
8. bash 08_convert_pdbqt.sh
9. ~/brd4_protonation_env/bin/python 06_heavy_atom_audit.py --stage pdbqt
10. ~/brd4_protonation_env/bin/python 09_restore_waters.py
11. ~/brd4_protonation_env/bin/python 10_water_qc.py
12. bash 11_prepare_redocking.sh, then manually submit REDOCK_PH75.pbs from the PH75/REDOCKING directory with DOMAIN=BD1 and DOMAIN=BD2.

Use the explicit Python 3.9 virtual-environment interpreter, especially after loading MGLTools. Inspect outputs before each next stage. Scripts may replace generated files in PH75 if rerun; preserve prior results separately. No baseline receptor files are overwritten.

## Existing archive
- [Ligand preparation script](../prepare_ligands.py) and [candidate manifest](../BRD4_NP15_ligand_manifest.csv): candidate structures must be checked against original screening rows and exact source identifiers before final docking; internal SMILES/InChIKey consistency alone is insufficient.
- [PROPKA figures, tables, and notebook](../PROPKA_FIGURES/).
- [QC observations](QC_observations.md).
- Other project scripts are already in CHPC_LOO, CHPC_NP_SCREEN, CHPC_STEREO, and CHPC_Y_RANDOM at repository root.

## Limits and original scripts
This archive requires the existing raw/clean receptors, baseline matched-wet PDBQT files, WSH reference/prepared ligand, and original grid configs. Earlier receptor extraction, water matching, crystal-ligand preparation, original grid generation and symmetry-aware RMSD validation scripts were not available in full in this chat. The actual CHPC originals are needed for an exact archive; use 13_archive_original_scripts.sh to collect those scripts/configs for later upload. Do not substitute reconstructed scripts for historical originals without checking.

pH 7.5 is a preparation setting, not experimentally proven optimal binding pH. PROPKA was run on apo structures; histidine tautomer assignments are model-based. The BD2 OXT is reconstructed on the last resolved residue, not the full crystallographic construct. Repeat redocking must be validated before final screening. Raw PDB/PQR/PDBQT and log files remain on CHPC; this is a script archive, not a full data backup.
