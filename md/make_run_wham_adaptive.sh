#!/bin/bash
# Generate a run_wham.sh for a given system using the per-window adaptive
# equilibration trim (Codes/md/adaptive_equil_trim.py) instead of a flat
# global percentage, plus the same wham-new histogram/binning as before
# (0-105 A, 420 bins).
# Usage: make_run_wham_adaptive.sh DIST_DIR OUT_DIR
set -e
DIST_DIR="$1"
OUT_DIR="$2"
mkdir -p "$OUT_DIR"
cat > "$OUT_DIR/run_wham.sh" <<EOF
#!/bin/bash
# WHAM analysis, PER-WINDOW adaptive equilibration trim (see
# Codes/md/adaptive_equil_trim.py for the algorithm/rationale), histogram
# 0-105 A / 420 bins.
# dist files: $DIST_DIR
set -e
python3 /gpfs/data/stippner/Codes/md/adaptive_equil_trim.py "$DIST_DIR" .

wham-new P 0 105 420 1e-6 300 0 metadata.dat pmf.dat

grep -v '^#' pmf.dat | awk '\$2==\$2+0 && \$2!="inf"{print \$1,\$2,\$4}' > _t
m=\$(awk 'NR==1||\$2<m{m=\$2}END{print m}' _t)
awk -v m=\$m '{printf "%.4f  %.4f  %.4e\n",\$1,\$2-m,\$3}' _t > pmf_shifted.dat && rm _t
EOF
chmod +x "$OUT_DIR/run_wham.sh"
echo "wrote $OUT_DIR/run_wham.sh"
