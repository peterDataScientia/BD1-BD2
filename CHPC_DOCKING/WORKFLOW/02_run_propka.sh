#!/bin/bash
set -euo pipefail
cd ~/BRD4_NP_DOCKING
mkdir -p 08_PROTONATION
cd 08_PROTONATION
for d in BD1 BD2; do
    cp "../02_RECEPTORS/${d}_clean.pdb" "${d}_apo.pdb"
    "$HOME/brd4_protonation_env/bin/propka3" "${d}_apo.pdb" > "${d}_propka.log" 2>&1
    echo "Completed: $d"
    sed -n '/SUMMARY OF THIS PREDICTION/,$p' "${d}_apo.pka" | head -n 80 || true
done
ls -lh *.pka *.log
