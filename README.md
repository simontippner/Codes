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
