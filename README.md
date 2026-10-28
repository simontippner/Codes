# Codes

Personal code repository for computational chemistry / data analysis scripts
(MD analysis, QM job setup, plotting, etc.).

## md/ — AMBER/cpptraj trajectory analysis

- **`pmf_from_pairdist.py`** — approximate CAT...PS binding free energy
  (potential of mean force) from an EXISTING unbiased MD trajectory's
  pairwise center-of-mass distance distribution -- no new simulation
  needed, just re-analysis of the `com_distances.dat` files this project's
  aggregation pipeline (`01_generate_cpptraj/gen_mix_cpptraj.py` /
  `gen_mix10_cpptraj.py`) already writes. For N independent CAT...PS pairs
  sampled every frame, the expected pair count in a shell under a uniform
  ("ideal gas") reference is `n_frames * N_pairs * shell_volume / V_box`;
  the ratio of observed to expected count is g(r) (a body-center-of-mass
  analogue of a site-site RDF), and PMF(r) = -RT ln[g(r)], shifted so the
  plateau at large r (dissociated reference state) sits at zero. The depth
  of the resulting well at the contact-distance minimum is the pairwise
  association free energy, estimated directly from how often the unbiased
  trajectory visits close vs. far separations -- legitimate when (as here)
  the trajectory already samples both bound and unbound configurations,
  but **not** a substitute for a proper umbrella-sampling/WHAM PMF along a
  defined reaction coordinate, nor a true single-pair/infinite-dilution
  binding free energy (this project's many-body 10:10 / 10:1 systems have
  competing pairs, not one isolated pair). Usage:
  ```
  python3 pmf_from_pairdist.py com_distances.dat BOX_VOLUME_A3 \
      --pair-regex 'd_(CAT\d*)_(RUB\d+)$' --r-max 25 --out-csv pmf.csv
  ```

- **`com_density_single_anchor.py`** — single-colormap solute density
  heatmap (the plain style, before `com_density_timeseries.py` grew
  species-coloring), but centered on **one specific body** (e.g. one chosen
  PS molecule) rather than the solute assembly's collective centroid.
  `com_density_timeseries.py` deliberately avoids pinning any one body at
  the origin (diagnosed as a reference-choice artifact when it happened
  implicitly via cpptraj's `autoimage anchor`); this script does that ON
  PURPOSE -- no centroid correction afterward, so the chosen anchor body
  sits fixed at the origin (marked with a star) in every panel by
  construction. That answers a genuinely different question ("what does
  the local environment around this one specific molecule look like") from
  "where does the assembly condense relative to its own center" -- pick
  the right script for the question being asked. Defaults to short (e.g.
  0.1 ns) windows for near-instantaneous snapshots rather than long (e.g.
  40 ns) pooled ones -- note short windows are necessarily sparse (few
  frames x few bodies), showing up as a speckled rather than smooth
  density; that's the expected tradeoff for temporal sharpness, not a bug.
  Usage:
  ```
  python3 com_density_single_anchor.py PRMTOP TRAJ \
      --body-mask ":1-10" --body-mask ":11-20" ... \
      --anchor-mask ":101" --n-panels 6 --window-ns 0.1 --out out.png
  ```

- **`com_density_timeseries.py`** — visualizes solute aggregation over an MD
  trajectory as a 2D spatial density heatmap of solute bodies' centers of
  mass, binned into time windows, one or more systems as separate rows (e.g.
  several compositions stacked on the same time-window and color axes). Each
  body is tagged with a species (`--body-mask SPECIES MASK`); one species is
  the **field** (`--field-species NAME`, a standard single-colormap density
  heatmap with a normal log-scale colorbar) and the other(s) are **marked**
  (`--mark-species NAME --mark-color HEX`, each body's actual per-frame
  position plotted directly as scatter points on top). An earlier version
  tried encoding both species as a single blended hue (density-weighted
  color mix) — technically correct but hard to read at a glance, since a
  flat region of one blended tone doesn't obviously parse as "a mix of two
  things" the way a colormap + distinct marker does. Marking is also the
  natural choice whenever one species has very few bodies (a single catalyst
  cluster against ten photosensitizer copies, say): a "density" of one body
  is a weak concept, whereas marking its actual position is direct and
  unambiguous.

  Extracts each body's per-frame COM via cpptraj (`vector ... center`), then
  **re-centers every body on the solute assembly's own collective centroid
  each frame** before histogramming — important: naively plotting raw (or
  single-body-anchored) COM positions conflates "where is the assembly
  relative to one arbitrarily pinned reference body" with genuine spread/
  condensation, since cpptraj's `autoimage anchor <mask> origin` (needed to
  resolve periodic-boundary consistency across bodies) pins that anchor body
  at the origin in every frame by construction — recentering on the group's
  own centroid afterward removes that artifact. Requires cpptraj (AmberTools)
  on PATH. Usage (PS density as the field, CAT positions marked — this
  project's convention):
  ```
  python3 com_density_timeseries.py \
      --system "0% water" PRMTOP1 TRAJ1 \
      --system "100% water" PRMTOP2 TRAJ2 \
      --field-species PS --mark-species CAT --mark-color "#386641" \
      --body-mask CAT ":1-10" --body-mask CAT ":11-20" ... \
      --body-mask PS ":101" --body-mask PS ":102" ... \
      --anchor-mask ":1-10" --windows 5 --dt-ps 10 --out out.png
  ```
  `--mark-species`/`--mark-color` are optional — omit them for a plain
  single-species density heatmap with no marker overlay. (Fixed bug: the
  marker legend was built unconditionally with `ncol=len(mark_species)`,
  which crashes `numpy.array_split` when that's 0 — i.e. exactly the
  omitted-marker case the option is meant to support. Guarded on
  `legend_handles` being non-empty.)

  The per-frame COM vectors and the cpptraj input that produced them used
  to live only in a `tempfile.TemporaryDirectory()` and vanish the moment
  the script exited, with no way to audit or re-plot a heatmap without
  re-running cpptraj from scratch. Now persisted for each system into
  `<that system's own MD dir>/analysis/density_vectors/<tag>.dat` and
  `.in` (`<tag>` = the basename of `--out`, so e.g. several different
  heatmap runs against the same system — `com_density_species_all7`,
  `solvent_density_all7`, `nh4_density_all7`, ... — don't overwrite each
  other, each landing right next to that system's other analysis output
  like `rdf_*.dat`/`com_distances.dat`). The saved `.in` file uses a bare
  relative `com.dat` output filename (cpptraj run with `cwd=workdir`), not
  the throwaway tempdir's absolute path, so it's directly re-runnable as-is
  — `cd` into the saved `density_vectors/` folder and `cpptraj -i <tag>.in`
  regenerates `com.dat` right there, no path editing needed.

- **`solvent_density_timeseries.py`** — the solvent-side companion: same
  time-windowed, centroid-recentered 2D density heatmap, but for solvent
  molecules (e.g. water/methanol O) instead of solute bodies, to see
  whether solvent is depleted from the region where solutes condense
  (classic hydrophobic-collapse signature) or stays uniform. Extracting a
  per-molecule COM the way the solute script does doesn't scale to
  thousands of solvent molecules, so instead one cpptraj pass computes the
  solute centroid (same `vector ... center` calls, used only for
  recentering, not plotted) **and** strips the trajectory down to just the
  solvent's representative atom, writing it out as a small NetCDF read
  directly with `scipy.io.netcdf_file` — far cheaper than thousands of
  individual `vector` calls. Defaults to a **linear** color scale, not log:
  solvent density only varies by a modest factor between bulk and a
  depleted pocket near a solute cluster (not orders of magnitude the way
  the sparse solute case does), so log washes the contrast out under the
  huge uniform bulk value. Hit (and fixed) a real cpptraj mask-parser bug
  along the way: `strip !((:WAT@O)|(:MOH@O1))` throws `Mask::ToRPN:
  unbalanced parentheses in expression` even though the parens are balanced
  by any normal count — wrapping `!(...)` around an inner expression that
  itself has parens around each term breaks the parser; `strip
  !(:WAT@O|:MOH@O1)` (no inner parens) works. Requires cpptraj and scipy on
  PATH. Usage:
  ```
  python3 solvent_density_timeseries.py \
      --system "0% water" PRMTOP1 TRAJ1 \
      --solute-body-mask ":1-10" --solute-body-mask ":11-20" ... \
      --solvent-mask ":WAT@O" --solvent-mask ":MOH@O1" \
      --anchor-mask ":1-10" --windows 5 --box-half 27 --out out.png
  ```
  Same vector-persistence behavior as `com_density_timeseries.py` above:
  the solute-body COM vectors (`com.dat`) and the cpptraj input that
  produced them are saved into each system's own `MD/analysis/
  density_vectors/<tag>.dat`/`.in`. The large stripped solvent trajectory
  itself (`solvent.nc`, one frame × every solvent/ion atom — can be
  hundreds of MB to GB-scale at high water content) is deliberately **not**
  kept; it's deleted right after use and is fully regenerable from the
  saved `.in` file run against the same production trajectory, so nothing
  is lost for reproducibility, just disk space.

- **`com_3d_snapshots.py`** — the 3D companion: matplotlib has no good
  volumetric/density rendering, so rather than a pooled 3D density this
  plots a handful of representative single frames (e.g. early/mid/late)
  side by side in 3D, colored by species — far more legible than a dense 3D
  point cloud would be. Same cpptraj COM extraction and centroid-
  recentering as `com_density_timeseries.py`. Usage:
  ```
  python3 com_3d_snapshots.py PRMTOP TRAJ \
      --body-mask CAT ":1-10" --body-mask CAT ":11-20" ... \
      --body-mask PS ":101" --body-mask PS ":102" ... \
      --species-color CAT "#386641" --species-color PS "#7798ab" \
      --anchor-mask ":1-10" --n-snapshots 4 --out out.png
  ```

- **`compare_pmf.py LABEL=FILE [LABEL=FILE ...] [--out prefix] [--no-shift]`**
  — overlays several WHAM PMF curves (column 1 = z, column 2 = PMF, e.g.
  Grossfield `pmf.dat` / `pmf_shifted.dat`) on one plot, without error bands,
  each shifted to its own minimum. Handy for comparing an old analysis
  with a re-run on new data.

  ```bash
  python3 md/compare_pmf.py "old=wham_analysis/pmf_shifted.dat" \
      "new=wham_analysis_2026-10-23/pmf_shifted.dat" --out pmf_old_vs_new
  ```

- **`combined_orientation_allsolvents.py`** — umbrella-window orientation of
  BODIPY (dipyrrin core and phenyl planes) and Mo3S13 (Mo3 plane) vs height in
  the brush, for 4 conditions (partially/fully protonated brush x 50:50
  H2O:MeOH / 100% H2O). Reads cpptraj `vector ... corrplane` outputs
  (`rc.dat`, `core.dat`, `phenyl.dat`, `norm.dat`) and writes P(theta|z) heat
  maps, mean-angle and order-parameter S(z) profiles, and a text summary of
  brush/solvent averages with window-level SE (200 frames per window).
  System paths are hardcoded at the top.

  ```bash
  python3 md/combined_orientation_allsolvents.py   # run in the output folder
  ```

- **`combined_orientation_v3_both.py`** — publication-style version of the
  above: P(theta|z) heat maps with structure icons, 4 rows (partially/fully
  protonated brush x 50:50 H2O:MeOH / 100% H2O) x 3 planes, one shared magma
  colour scale so mixture and water panels are directly comparable. Panels
  are lettered a-l (row-major); row conditions go in the caption.

  ```bash
  python3 md/combined_orientation_v3_both.py   # run in the output folder
  ```

- **`combined_v2_inbrush.py`** — PS (BODIPY) + CAT (Mo3S13) donor-H solvation
  figure (brush-solvent, solute-brush and solute-solvent RDFs plus a
  preferential-solvation map) using only the umbrella windows with the solute
  inside the brush (z < 58 A = frames 1-23200 of all.nc). Reads cpptraj
  `radial` outputs from `.../traj/analysis/inbrush/{solute_solvent_H,solute_brush}`.

  ```bash
  python3 md/combined_v2_inbrush.py   # run in the output folder
  ```

- **`combined_v2_water_half_inbrush.py`, `combined_v2_water_full_inbrush.py`** —
  the same brush-only RDF figure for the water-only runs (partially / fully
  protonated brush; no methanol, so no preferential-solvation map). Reads
  `.../traj/analysis/inbrush/{solvent_brush,solute_brush,solute_solvent}`.

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

- **`make_qmmm_frame.py PRMTOP TRAJ FRAME --qm-mask MASK --out-prefix
  PREFIX [--shell 18.0] [--anchor-mask] [--fixed-mask]`** — electrostatic-
  embedding QM/MM frame extraction from an existing AMBER MD trajectory:
  autoimages/recenters on the QM region, truncates the rest to a finite MM
  shell (default 18 Å — ORCA has no PBC, so the whole periodic box can't be
  embedded), then writes the QM-region atoms as a plain `.xyz` and every
  kept MM atom's real `system.prmtop` charge + coordinate as an ORCA
  point-charge `.pc` file. No link-atom support — only for a QM region with
  nothing covalently bonded across the boundary (a non-covalent donor...
  acceptor pair in explicit solvent, e.g.). Use a `^`-prefixed (molecule
  number) `--qm-mask` when the QM region is one or more whole AMBER
  "molecules" with multi-residue internal topology (e.g. a metal cluster
  split into several residues) — a `:`-residue mask would only grab one
  sub-residue. See `Luca/05_qmmm_ensemble/qmmm_ensemble_workflow.pdf` for
  the full design rationale.
  Example: `python3 make_qmmm_frame.py system.prmtop 04_prod.nc 4990
  --qm-mask "^1,2" --out-prefix frame499`

- **`setup_qmmm_ensemble.py --system LABEL PRMTOP TRAJ [--system ...]
  --qm-mask MASK [--qm-resmask MASK] [--shell 18.0] [--qm-charge 0]
  [--qm-mult 1] [--n-frames 2000] [--nprocs 16] [--out-dir QM_ensemble_SPE]
  [--workers 16] [--sander-only]`** — batch-builds a QM/MM snapshot
  ensemble from an existing AMBER MD trajectory, two stages per `--system`:
  Stage A is one single cpptraj pass — `trajin` every `--n-frames`
  equally-spaced frame, `autoimage anchor <CAT/MoS> origin`, then
  `trajout ... ncrestart multi` to dump each selected frame straight to
  its own full-system restart (cpptraj numbers these sequentially by
  selection order, matching the 1-based `FRAMES/<i>` layout used
  throughout) — persisted, not a throwaway tempdir. Stage B runs
  `make_qmmm_frame.py` on each small per-frame restart (not the original
  multi-GB trajectory) in parallel (`ProcessPoolExecutor`, `--workers`
  processes) to also write the shell-truncated QM `.xyz` + MM point-charge
  `.pc`, plus a ready-to-run plain-ORCA `.inp` (`%pointcharges`,
  B3LYP-D3(BJ)/def2-TZVP/RIJCOSX/def2-J — same level as the rest of
  `04_uvvis` but WITHOUT CPCM, since the point charges already supply the
  solvent electrostatics). Passing `--qm-resmask` (a standard AMBER
  residue mask, e.g. `":1-11"` — sander doesn't understand `--qm-mask`'s
  `^`-molecule-number syntax) additionally writes a shared `qmmm.in`
  (sander `&cntrl`/`&qmmm`/`&orc`, `ifqnt=1`, `qm_theory='EXTERN'`) +
  `orc_job.tpl` (the ORCA method block, via `use_template=1`) per system,
  symlinked together with `system.prmtop` into every `FRAMES/<i>/` — the
  full-box electrostatic-embedding alternative to the point-charge/.pc
  path, matching this group's established Redoxpotentials QM-ensemble
  project layout. `--sander-only` (re)writes just that qmmm.in/orc_job.tpl
  step against an already-extracted `FRAMES/` tree, skipping Stage A/B
  entirely. Every frame folder's provenance (which original trajectory
  frame `i` maps to, and the exact mask/charge/method used) is recorded in
  that system's `manifest.txt`. Also writes (but does NOT submit) a
  `submit_array.slurm` SLURM array for the plain-ORCA `.inp` path, reusing
  `orca_submit.slurm` unmodified per array task — this script only ever
  touches the filesystem, never the queue. Run the RMSD autocorrelation of
  the system first (see `Luca/05_qmmm_ensemble/decorrelation/`) to know
  how many of `--n-frames` are actually statistically independent — equal
  spacing alone does not guarantee that.
  Example: `python3 setup_qmmm_ensemble.py --system "0% water, buf=exp"
  system.prmtop 04_prod.nc --qm-mask "^1,2" --qm-resmask ":1-11" --n-frames 2000`

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
