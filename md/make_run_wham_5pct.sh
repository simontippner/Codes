#!/bin/bash
# Generate a run_wham.sh for a given system that:
#  - discards the first 5% of every window as equilibration (EQUIL_FRAC=0.05,
#    per explicit user request -- narrower than the 10% trim used earlier this
#    session for MoS_in_brush RETRY_{full,half}_water)
#  - calls wham-new with histogram range 0-105 A, 420 bins (per explicit user
#    request, replacing whatever params/range the system's existing
#    wham_analysis/ used before)
# Usage: make_run_wham_5pct.sh DIST_DIR OUT_DIR
#   DIST_DIR : directory containing dist_XXX.X.dat (the ../dist_*.dat source)
#   OUT_DIR  : new wham_analysis_5pct directory to create and write into
set -e
DIST_DIR="$1"
OUT_DIR="$2"
mkdir -p "$OUT_DIR"
cat > "$OUT_DIR/run_wham.sh" <<EOF
#!/bin/bash
# WHAM analysis, 5% equilibration trim, histogram 0-105 A / 420 bins (user-specified).
# dist files: $DIST_DIR
set -e
mkdir -p ts
: > metadata.dat
EQUIL_FRAC=0.05
for f in $DIST_DIR/dist_*.dat; do
    b=\$(basename "\$f" .dat); c=\${b#dist_}
    center=\$(awk -v c="\$c" 'BEGIN{printf "%.4f", c+0}')
    n=\$(wc -l < "\$f")
    skip=\$(awk -v n="\$n" -v ef="\$EQUIL_FRAC" 'BEGIN{printf "%d", n*ef}')
    awk -v s="\$skip" 'NR>s{print \$1"\t"\$7}' "\$f" > "ts/win_\${c}.dat"
    printf "ts/win_%s.dat\t%s\t10.0\n" "\$c" "\$center" >> metadata.dat
done

wham-new P 0 105 420 1e-6 300 0 metadata.dat pmf.dat

grep -v '^#' pmf.dat | awk '\$2==\$2+0 && \$2!="inf"{print \$1,\$2,\$4}' > _t
m=\$(awk 'NR==1||\$2<m{m=\$2}END{print m}' _t)
awk -v m=\$m '{printf "%.4f  %.4f  %.4e\n",\$1,\$2-m,\$3}' _t > pmf_shifted.dat && rm _t
EOF
chmod +x "$OUT_DIR/run_wham.sh"
echo "wrote $OUT_DIR/run_wham.sh"
