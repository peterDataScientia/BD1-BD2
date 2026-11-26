#!/bin/bash
set -euo pipefail
cd ~/BRD4_NP_DOCKING/08_PROTONATION
mkdir -p PH75
PQR="$HOME/brd4_protonation_env/bin/pdb2pqr"
"$PQR" --version
for d in BD1 BD2; do
    echo "===== Preparing $d at pH 7.5 ====="
    if "$PQR" --ff=AMBER --ffout=AMBER --keep-chain \
        --titration-state-method=propka --with-ph=7.5 \
        --pdb-output="PH75/${d}_pH75.pdb" "${d}_apo.pdb" \
        "PH75/${d}_pH75.pqr" > "PH75/${d}_preparation.log" 2>&1; then
        echo "PASS: $d"
    else
        tail -n 25 "PH75/${d}_preparation.log"
        exit 1
    fi
done
ls -lh PH75/
