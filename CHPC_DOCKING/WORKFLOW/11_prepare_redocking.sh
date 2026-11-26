#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd ~/BRD4_NP_DOCKING/08_PROTONATION
mkdir -p PH75/REDOCKING
cd PH75/REDOCKING
for d in BD1 BD2; do
    cp "../PDBQT/${d}_pH75_wet.pdbqt" .
    cp "../../../04_REDOCKING/${d}_WSH.pdbqt" .
    cp "../../../04_REDOCKING/${d}_grid.conf" .
done
cp "$HERE/REDOCK_PH75.pbs" .
ls -lh
# Explicit submission only, after checking inputs:
# qsub -v DOMAIN=BD1 REDOCK_PH75.pbs
# qsub -v DOMAIN=BD2 REDOCK_PH75.pbs
# qstat -u "$USER"
