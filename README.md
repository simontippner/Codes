# Codes

Personal code repository for computational chemistry / data analysis scripts
(MD analysis, QM job setup, plotting, etc.).

## qm/ — ORCA job setup and TD-DFT UV-vis post-processing

Scripts built around an ORCA 6 workflow: ground-state geometry optimization
(B3LYP-D3(BJ)/def2-SVP, CPCM solvent, RIJCOSX) followed by a single-point
TD-DFT calculation (def2-TZVP) on the optimized geometry to get vertical
excitation energies and oscillator strengths, from which a UV-vis spectrum
is simulated.

- **`orca_submit.slurm`** — generic SLURM submission script for ORCA 6 jobs
  on the cluster. Usage: `sbatch orca_submit.slurm <input>.inp`. Loads the
  `orca6` module (adding `/usr/license/modulefiles` to `MODULEPATH` first,
  since it isn't there by default), stages the job (input + all `.xyz`
  files in the submit directory) on node-local scratch
  (`/public/<host>/scratch/tmp/...`), runs ORCA there, and copies every
  output file back to the submit directory on success. On failure it
  appends an entry to an `I_CRASHED` file (timestamp, host, scratch path)
  instead of silently losing the job.

All the Python scripts below follow the same convention: `python3
script.py FILENAME [options]`, where `FILENAME` is the ORCA log/output file
(or, for `make_combined_xyz.py`, an xyz file) to read — no editing the
script or hardcoding paths required. Every one also works on a log file
from a job that's still running, since they only read the shared-filesystem
`.log` (never the job's node-local scratch directory).

- **`plot_uvvis.py FILENAME [--out <prefix>] [--fwhm 0.4] [--xmin 200]
  [--xmax 700] [--exp-ref nm nm ...]`** — parses the `ABSORPTION SPECTRUM
  VIA TRANSITION ELECTRIC DIPOLE MOMENTS` table out of an ORCA TD-DFT
  output/log file (state, energy, wavelength, oscillator strength for each
  root) and converts the discrete stick transitions into a continuous,
  Gaussian-broadened UV-vis spectrum using the standard
  `eps(E) = 1.3062974e8 * sum_i [f_i/FWHM] * exp(-4 ln2 ((E-E_i)/FWHM)^2)`
  convention (same one used by Multiwfn/GaussView-style UV-vis tools).
  Writes the raw transitions and the broadened spectrum to CSV, plus a
  PNG/PDF plot (broadened curve + stick spectrum overlay). `--out` defaults
  to the input filename's basename; pass `--exp-ref 286 450` to mark
  experimental reference wavelengths if you have them (none marked by
  default).
  Example: `python3 plot_uvvis.py rubpy3_uvvis.log`

- **`plot_opt_energy.py FILENAME [--out <prefix>]`** — plots the SCF energy
  at every cycle of an ORCA geometry optimization (parsed from each `TOTAL
  SCF ENERGY` block in the log), relative to the first cycle, in kcal/mol
  vs. cycle number. Useful for spotting real conformational jumps vs. just
  slow convergence — e.g. a sudden multi-kcal/mol step partway through
  usually means a genuine structural rearrangement, not numerical noise.
  Writes a CSV (cycle, energy in Eh/eV, relative kcal/mol) plus a PNG/PDF.
  Example: `python3 plot_opt_energy.py rubpy_mos_opt.log`

- **`make_trj_xyz.py FILENAME [--out trajectory.xyz]`** — extracts every
  `CARTESIAN COORDINATES (ANGSTROEM)` block from an ORCA geometry-
  optimization log file and writes them out as a standard multi-frame xyz
  trajectory (one frame per optimization cycle), readable directly in VMD
  or any other xyz-trajectory viewer.
  Example: `python3 make_trj_xyz.py rubpy_mos_opt.log`

- **`make_combined_xyz.py FRAGMENT1.xyz FRAGMENT2.xyz [--gap 6.0] [--out
  combined.xyz]`** — combines two independently-optimized molecular
  fragments into a single multi-fragment `.xyz` file, placing the second
  fragment along +x at a safe, non-clashing center-to-center distance (sum
  of each fragment's own bounding-sphere radius, plus `--gap`). Used to
  build a starting geometry for a combined/supramolecular system (e.g. a
  photosensitizer + a catalyst) out of two separately-optimized pieces,
  ready for a joint geometry optimization.
  Example: `python3 make_combined_xyz.py rubpy3_opt.xyz mo3s13_opt.xyz --out rubpy_mos.xyz`

## plotting/ — figure-making helpers

- **`sanzo_wada_palette.py`** — browse and export color palettes from Sanzo
  Wada's 1933 *A Dictionary of Colour Combinations* for use in publication
  figures. Uses the open-sourced digitization of the book (159 named colors,
  348 curated 2-4 color combinations; MIT-licensed dataset from
  [mattdesl/dictionary-of-colour-combinations](https://github.com/mattdesl/dictionary-of-colour-combinations),
  cached locally as `sanzo_wada_colors.json` so the script works offline
  after the first run). Subcommands:
  - `list [--n-colors 2|3|4] [--search <name>]` — list/filter the 348 combinations
  - `show <id>` — print one combination's hex/RGB/name in the terminal
  - `export <id> --out <prefix> [--cvd-preview]` — render a clean swatch
    figure (PNG+PDF) plus ready-to-paste code: a Python hex list, a LaTeX
    `xcolor` `\definecolor` block, and a raw JSON dump. `--cvd-preview` adds
    protanopia/deuteranopia/tritanopia simulated rows underneath the
    original swatch as a quick (approximate) accessibility sanity check —
    these are historical aesthetic combinations, not validated
    colorblind-safe data-viz palettes, so treat this as a first-pass check,
    not a substitute for a proper categorical-palette validator.
  - `random [--n-colors N] [--out <prefix>]` — pick (and optionally export)
    a random combination, e.g. for picking an accent pairing quickly.
  - `gallery [--n-colors N] [--search <name>] [--limit 30] [--out <prefix>]`
    — render many combinations at once as a single browsable grid figure
    (one row per combination, labeled by id on the left), for visually
    scanning a filtered set instead of reading hex codes off `list`.

  Three commands follow the standard data-viz palette taxonomy (the same
  "job the color does" split used by the `dataviz` skill: sequential for
  magnitude, diverging for polarity, qualitative/categorical for identity):

  - `sequential (--color '#hex' | --from-combination <id>) [--steps 8]
    [--out <prefix>]` — one hue, light→dark: a magnitude colormap (e.g. for
    a heatmap), built by interpolating a single base hue between a near-
    white tint and a dark shade in CIELAB space (perceptually smoother than
    interpolating raw RGB).
  - `diverging (--colors '#hex1,#hex2' | --from-combination <id>)
    [--midpoint '#hex'] [--steps 9] [--out <prefix>]` — two hues + a
    neutral gray midpoint: a polarity colormap (e.g. positive/negative, or
    above/below a reference value). Never puts a hue at the midpoint.
  - `qualitative --n N [--out <prefix>] [--cvd-preview]` — pick N
    maximally-separated colors out of the full 159-color set via greedy
    farthest-point sampling in CIELAB space, for a categorical/identity
    palette rather than an aesthetically-coordinated combination.

  Both `sequential`/`diverging` render a continuous 256-step colorbar plus
  the requested discrete steps, and write a ready-to-use
  `LinearSegmentedColormap` to the Python snippet.

  - `cvd-test (--colors '#hex,...' | --from-combination <id>) [--out
    <prefix>]` — test **any** palette's real distinguishability under color
    vision deficiency: renders the protanopia/deuteranopia/tritanopia
    preview rows and prints the minimum pairwise CIE76 Delta-E, both for
    normal vision and under each simulated CVD type. This is a genuine
    numeric diagnostic, not just a picture to eyeball — it will explicitly
    warn when two colors are likely to collide (e.g. a red/green pair
    collapsing under deuteranopia), including for `qualitative` picks,
    which are CIELAB-separated under *normal* vision but aren't
    automatically guaranteed to stay separated under every CVD type.
  - `demo (--from-combination <id> | --colors '#hex,...') [--out <prefix>]`
    — render example line, scatter, and heatmap plots side by side, all
    using the same chosen palette (categorical colors for the lines/dots,
    the palette Lab-interpolated into a continuous colormap for the
    heatmap) — a quick way to see how a palette actually reads across
    different plot types before committing to it in a real figure.

  ### Example output

  All commands below are run from inside `plotting/`; the exact commands
  that produced every image here are also saved as a runnable script at
  [`examples/sanzo_wada/commands.sh`](plotting/examples/sanzo_wada/commands.sh)
  (`cd plotting && bash examples/sanzo_wada/commands.sh` reproduces all four).

  ```bash
  python3 sanzo_wada_palette.py sequential --color "#1c4286" --steps 8 \
      --out examples/sanzo_wada/sequential
  ```
  ![sequential example](plotting/examples/sanzo_wada/sequential.png)

  ```bash
  python3 sanzo_wada_palette.py diverging --colors "#cc1236,#00978d" --steps 9 \
      --out examples/sanzo_wada/diverging
  ```
  ![diverging example](plotting/examples/sanzo_wada/diverging.png)

  ```bash
  python3 sanzo_wada_palette.py qualitative --n 5 \
      --out examples/sanzo_wada/qualitative
  ```
  (note the honest deuteranopia/tritanopia warnings this particular pick
  gets — greedy CIELAB separation under normal vision doesn't guarantee CVD
  safety; see `cvd-test` above)

  ![qualitative example](plotting/examples/sanzo_wada/qualitative.png)

  ```bash
  python3 sanzo_wada_palette.py demo --from-combination 121 \
      --out examples/sanzo_wada/demo
  ```
  ![demo example](plotting/examples/sanzo_wada/demo.png)

  Each example folder also contains the matching `.pdf`, `_python.txt`
  (ready-to-paste hex list / colormap), `_latex.tex` (`xcolor`
  `\definecolor` block), and `.json` that the command itself writes.
