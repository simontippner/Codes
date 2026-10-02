#!/bin/bash
# Exact commands used to generate the example images in this folder.
# Run from inside plotting/ (one level up from here):
#   cd plotting && bash examples/sanzo_wada/commands.sh

set -e
OUT=examples/sanzo_wada

python3 sanzo_wada_palette.py sequential \
    --color "#1c4286" --steps 8 \
    --out "$OUT/sequential"

python3 sanzo_wada_palette.py diverging \
    --colors "#cc1236,#00978d" --steps 9 \
    --out "$OUT/diverging"

python3 sanzo_wada_palette.py qualitative \
    --n 5 \
    --out "$OUT/qualitative"

python3 sanzo_wada_palette.py demo \
    --from-combination 121 \
    --out "$OUT/demo"
