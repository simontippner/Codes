#!/usr/bin/env python3
"""Convert a VMD TachyonInternal .tga render to .png.

Usage: python3 tga2png.py IN.tga OUT.png
"""
import sys
from PIL import Image

Image.open(sys.argv[1]).save(sys.argv[2])
print(f"wrote {sys.argv[2]}")
