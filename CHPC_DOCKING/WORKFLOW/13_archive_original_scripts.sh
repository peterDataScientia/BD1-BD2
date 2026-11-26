#!/bin/bash
# Collect actual on-cluster script files, not credentials or the entire home directory.
# Run locally on CHPC. This does NOT push to GitHub.
set -euo pipefail
cd ~/BRD4_NP_DOCKING
out="$HOME/BRD4_NP_DOCKING_original_scripts.tar.gz"
find . -type f \( -name '*.py' -o -name '*.sh' -o -name '*.pbs' -o -name '*.conf' -o -name '*.ipynb' \) -print0 | tar --null -T - -czf "$out"
tar -tzf "$out"
sha256sum "$out"
echo "Archive: $out"
