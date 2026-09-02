#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Guided restoration of real input data into data/.
#
# Real microscopy TIFFs (and any Tutorials 2-5 Excel workbooks) are intentionally
# NOT tracked by git (see .gitignore: data/ is ignored). This script:
#   1. Creates the expected data/ directory structure.
#   2. Tries to recover the real z-stacks from git history where they lived as
#      images/<file> before they were moved under data/.
#   3. Prints guidance for placing data that cannot be recovered automatically.
#
# Usage:  pixi run data  (or: bash scripts/download_data.sh)
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$REPO_ROOT/data/images"
mkdir -p "$DEST"

# Files historically tracked under images/ before the data/ reorganization.
REAL_STACKS=(
  HCT116-Cells-Monolayer-C1.tif
  Human-Colon-Organoids-C1.tif
  PDAC-C1.tif
  PDAC-C2.tif
  ZeroG-Breast-Cancer-Spheroid-C1.tif
)

restored=0
for f in "${REAL_STACKS[@]}"; do
  if [[ -f "$DEST/$f" ]]; then
    echo "already present: $f"
    continue
  fi
  if git -C "$REPO_ROOT" cat-file -e "HEAD:images/$f" 2>/dev/null; then
    git -C "$REPO_ROOT" show "HEAD:images/$f" > "$DEST/$f"
    echo "restored from git (HEAD): $f"
    restored=$((restored + 1))
  else
    echo "not found in git history: $f"
  fi
done

echo ""
if (( restored > 0 )); then
  echo "Restored $restored file(s) into $DEST."
fi
echo ""
echo "---"
echo "Real microscopy data is deliberately excluded from version control."
echo "Place your actual z-stack TIFFs under:  $DEST"
echo "Place any Tutorials 2-5 Excel datasets under:  $REPO_ROOT/data/tutorial_excel/"
echo "Run after restoring:  pixi run smoke   (smoke segmentation on PDAC-C1/C2)"