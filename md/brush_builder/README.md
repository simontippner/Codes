# brush_builder

All-atom models of **randomly grafted, atactic polyelectrolyte brushes** for Amber, with periodic boundary
conditions imposed while the chains are packed. One Python tool, one configuration file.

Method: S. Tippner, D. Hernández-Castillo, F. H. Schacher, L. González,
*All-atom molecular dynamics simulations of grafted poly(N,N-dimethylaminoethyl methacrylate) brushes*,
J. Phys. Chem. B **129**, 2105–2114 (2025), https://doi.org/10.1021/acs.jpcb.4c07928.
This is a rewrite of the scripts used there; the procedure is the same (see *Differences* below).

## Steps

| Step | Subcommand | What it does |
|---|---|---|
| 1 | `params` | head, middle and tail residue of every monomer type from a RESP-charged antechamber file (`prepgen`), GAFF2 parameters (`parmchk2`) |
| 2 | `chains` | `n_chains` random sequences: `round(alpha * n_units)` protonated units at random positions, random stereochemistry; built with `tleap`, aligned along z with the grafted unit at the bottom |
| 3 | `pack` | (a) smallest square box holding all chains upright at the PACKMOL `tolerance` (minimum distance between chains, which sets the grafting density), widened by `box_step`; (b) packing in a region of 3 × 3 cells: each chain is placed with its grafted atom on the plane `anchor_z` inside the central cell, against all chains placed before and their copies, then copied into the 8 neighbouring cells (copies kept if any atom lies within half a cell of the central cell). After `max_attempts` failures for one chain the packing restarts with chain 1; after `max_fail` failed packings the region grows by `region_step` |
| 4 | `solvate` | the whole region, copies included, is filled with water (`shell_xy` beyond the polymer laterally, `shell_z` above it; number of water molecules from the `density` and the solute mass), counterions go into the central cell; then the copies and the water outside the central cell are removed and `tleap` writes the periodic system (GAFF2, TIP3P) |
| 5 | `single` | optional: every chain solvated on its own with its counterions (reference systems) |

`all` runs steps 1 to 4. Because chains and water are packed around the copies, the system is consistent across
the box boundaries from the start: no overlaps and no gaps at the faces of the simulation box.

## Requirements

- Python ≥ 3.11 with numpy
- AmberTools (`prepgen`, `parmchk2`, `tleap`) on `PATH`, e.g. `source $AMBERHOME/amber.sh`
- PACKMOL 20.x (path in the configuration, `[tools] packmol`, or on `PATH`)

## Usage

```bash
cd examples/pdmaema
python3 ../../brush_builder.py config_quick.toml all     # 4 short chains, a few minutes
python3 ../../brush_builder.py config.toml all           # 16 x 30-mers as in the paper (packing can take hours)
python3 ../../brush_builder.py config.toml single        # optional single-chain systems
```

Steps can be run one at a time (`params`, `chains`, `pack`, `solvate`); each reads the output of the previous one
from the output folder.

## Output (`<output>/`)

| Path | Content |
|---|---|
| `params/` | `*.prepi`, `*.frcmod`, generated `mainchain.*` files |
| `chains/chain_<i>.pdb`, `chains/sequences.json` | aligned chains and their sequences |
| `pack/region.pdb`, `pack/central.pdb`, `pack/pack.json` | packed chains with copies, central chains, cell size and grafting density |
| `solvate/system.prmtop`, `system.inpcrd`, `system.pdb`, `system.json` | the periodic simulation system (box in `system.json`) |
| `single/chain_<i>/system.prmtop`, ... | single-chain systems |

## Example: PDMAEMA (`examples/pdmaema`)

- `monomers/dda.ac`, `ddb.ac` (unprotonated), `dpa.ac`, `dpb.ac` (protonated): the two stereoisomers of the DMAEMA
  unit with RESP charges from B3LYP-D3(BJ)/def2-SVP Merz–Kollman potentials (Gaussian inputs in `qm/`; charges
  fitted with `antechamber -fi gesp -c resp -at gaff2`)
- `config.toml`: settings of the paper (16 chains of 30 units, α = 0.5). The tolerance sets the grafting density:
  2, 5 and 10 Å gave 0.68, 0.39 and 0.21 chains per nm²
- `config_quick.toml`: 4 chains of 8 units at a tolerance of 5 Å, to test the whole pipeline

Other monomers only need their antechamber files and the atom names in `[residues]` (backbone head and tail atom,
end group kept by the head unit, H cap kept by the tail unit).

## Differences to the original scripts

- One tool and one configuration instead of separate scripts per step and per system; the 12 `mainchain` files are
  generated from the atom names in `[residues]`. The residues are identical to the published ones (checked: atoms,
  types and charges).
- Chains are aligned along z by rotating their longest axis onto z. The original fitted each chain into a narrow,
  tall PACKMOL box, widened in steps of 0.1 Å; since PACKMOL moves molecules as rigid bodies, the result is the same
  rotation, but it took many PACKMOL runs.
- The packing no longer adds a duplicate copy of each chain at zero shift (it was removed again before the
  simulation and did not change the systems).
- A copy is kept if any of its atoms lies within half a cell of the central cell (the original tested only the last
  atom of the chain), and the plane of the grafted atoms in step 3a spans the full box (the original used the box
  height as its y limit).
- Solvation is done directly with PACKMOL (water and counterions placed around the fixed polymer), as the
  `solvate.tcl` script of PACKMOL did in the original workflow, with the same first, mass-based water count.
- New: the water count is corrected until the density of the free water layer above the brush is within 3 % of the
  target (the polymer excludes more volume than its mass suggests; in the 4-chain example the first pass gave
  1.26 g/cm³, the second 1.00 g/cm³).
- New: no clashes across the periodic boundaries. The water of the neighbouring cells is not a periodic copy of the
  central cell, so after the cut water molecules near opposite faces can overlap; they are removed (closer than 2 Å
  to another molecule under periodic boundaries), the box gets a gap of 2 Å at the top, and counterions stay 1 Å
  from the cell faces. In the original workflow such contacts (including Cl⁻ on the polymer) were left to the first
  minimisation.

## Checks

`examples/pdmaema/config_quick.toml` (4 chains of 8 units, tolerance 5 Å, AmberTools 23, PACKMOL 20.4.0):
residues identical to the published ones; final box 28.0 × 28.0 × 58.0 Å with 4 chains, 16 Cl⁻ and 853 water
molecules, total charge 0, no heavy atoms of different molecules closer than 2 Å (periodic boundaries included),
water density above the brush 1.01 g/cm³.
