#!/bin/bash
set -euo pipefail
cd ~/BRD4_NP_DOCKING/08_PROTONATION
for d in BD1 BD2; do
    echo "===== $d HISTIDINE STATES ====="
    awk '/^ATOM/ && $4 ~ /^(HIS|HID|HIE|HIP|HSD|HSE|HSP)$/ && $3 ~ /^(ND1|NE2|HD1|HE2)$/' "PH75/${d}_pH75.pdb"
    echo "===== $d TERMINAL OXYGEN ====="
    grep "^ATOM" "PH75/${d}_pH75.pdb" | grep " OXT "
    echo "===== $d WARNINGS ====="
    grep -Ei 'warning|error|missing|unable|failed|repair|rebuild' "PH75/${d}_preparation.log" | tail -n 40 || true
done
