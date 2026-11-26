#!/bin/bash
set -euo pipefail
cd ~/BRD4_NP_DOCKING/02_RECEPTORS
for d in BD1 BD2; do
    echo "===== $d HISTIDINES ====="
    awk '/^ATOM/ && $4 ~ /^(HIS|HID|HIE|HIP)$/ {print}' "${d}_wet_matched_receptor.pdbqt"
    echo "===== $d C-TERMINUS ====="
    grep "^ATOM" "${d}_clean.pdb" | tail -n 15
done
grep "^REMARK 465" ../01_RAW/2YEM.pdb
