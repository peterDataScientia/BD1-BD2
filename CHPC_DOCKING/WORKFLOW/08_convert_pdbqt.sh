#!/bin/bash
set -euo pipefail
cd ~/BRD4_NP_DOCKING/08_PROTONATION
module load chpc/mgltools/1.5.6/gcc-7.3.0
MGL=/apps/chpc/chem/mgltools/1.5.6
ADT=$MGL/MGLToolsPckgs/AutoDockTools/Utilities24
mkdir -p PH75/PDBQT
for d in BD1 BD2; do
    "$MGL/bin/pythonsh" "$ADT/prepare_receptor4.py" \
        -r "PH75/${d}_pH75.pdb" -o "PH75/PDBQT/${d}_pH75_dry.pdbqt" \
        -A None -U nphs_lps_waters > "PH75/PDBQT/${d}_conversion.log" 2>&1
    echo "===== $d HISTIDINE QC ====="
    awk '/^ATOM/ && $4 ~ /^(HIS|HID|HIE|HIP)$/ && $3 ~ /^(ND1|NE2|HD1|HE2)$/' "PH75/PDBQT/${d}_pH75_dry.pdbqt"
    cat "PH75/PDBQT/${d}_conversion.log"
done
# MGLTools changes PATH to Python 2.5.6 on Lengau. Use the venv executable for QC.
