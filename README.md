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

- **`plot_uvvis.py`** — parses the `ABSORPTION SPECTRUM VIA TRANSITION
  ELECTRIC DIPOLE MOMENTS` table out of an ORCA TD-DFT output/log file
  (state, energy, wavelength, oscillator strength for each root) and
  converts the discrete stick transitions into a continuous, Gaussian-
  broadened UV-vis spectrum using the standard
  `eps(E) = 1.3062974e8 * sum_i [f_i/FWHM] * exp(-4 ln2 ((E-E_i)/FWHM)^2)`
  convention (same one used by Multiwfn/GaussView-style UV-vis tools).
  Writes the raw transitions and the broadened spectrum to CSV, plus a
  PNG/PDF plot (broadened curve + stick spectrum overlay, with optional
  experimental reference wavelengths marked).
  Usage: `python3 plot_uvvis.py --log <job>.log --out <prefix> [--fwhm 0.4]
  [--xmin 200] [--xmax 700] [--exp-ref 286 450]`.

- **`make_combined_xyz.py`** — combines two independently-optimized
  molecular fragments (each read from its own ORCA-optimized `.xyz`) into
  a single multi-fragment `.xyz` file, placing the second fragment along
  +x at a safe, non-clashing center-to-center distance (sum of each
  fragment's own bounding-sphere radius, plus a configurable gap). Used to
  build a starting geometry for a combined/supramolecular system (e.g. a
  photosensitizer + a catalyst) out of two separately-optimized pieces,
  ready for a joint geometry optimization. Paths to the two input
  fragments, the gap distance, and the output path are set as constants at
  the top of the script.

- **`make_trj_xyz.py`** — extracts every `CARTESIAN COORDINATES (ANGSTROEM)`
  block from an ORCA geometry-optimization log file and writes them out as a
  standard multi-frame xyz trajectory (one frame per optimization cycle,
  readable directly in VMD or any other xyz-trajectory viewer). Works even
  while the optimization job is still running, since it only needs the
  live-written `.log` file on the shared filesystem — not the job's
  node-local scratch directory, which usually isn't accessible until the
  job finishes. Usage: `python3 make_trj_xyz.py <job>.log <output>.xyz`.

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
  - `gradient (--colors '#hex,#hex,...' | --from-combination <id>)
    [--steps 8] [--out <prefix>]` — build a smooth, perceptually-uniform
    gradient/colormap by interpolating through 2+ anchor colors in CIELAB
    space (much smoother than interpolating raw RGB). Renders a continuous
    256-step colorbar plus the requested number of discrete steps, and
    writes a ready-to-use `LinearSegmentedColormap` to the Python snippet.
  - `distinct --n N [--out <prefix>] [--cvd-preview]` — pick N
    maximally-separated colors out of the full 159-color set via greedy
    farthest-point sampling in CIELAB space, for a categorical/identity
    palette rather than an aesthetically-coordinated combination. Prints
    the minimum pairwise CIE76 Delta-E among the picks, both for normal
    vision and simulated under each CVD type — this is a real diagnostic
    (it will honestly report a low Delta-E if two picks collide under e.g.
    deuteranopia), not just a picture to eyeball.
