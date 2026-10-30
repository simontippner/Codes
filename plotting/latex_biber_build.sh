#!/bin/bash
# Build 0main.pdf: lualatex -> biber -> lualatex x2. Usage: ./build.sh
cd "$(dirname "$0")"
lualatex -interaction=nonstopmode 0main.tex > /dev/null
biber 0main > /dev/null
lualatex -interaction=nonstopmode 0main.tex > /dev/null
lualatex -interaction=nonstopmode 0main.tex > /dev/null
echo "errors:   $(grep -c '^!' 0main.log)"
echo "undefined: $(grep -cE 'undefined|Citation .* undefined' 0main.log)"
