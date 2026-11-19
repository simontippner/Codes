# brush_builder

All-atom models of **randomly grafted, atactic polyelectrolyte brushes** for Amber, in water or solvent mixtures,
with periodic boundary conditions imposed while the chains and counterions are packed. One Python tool, one
configuration file.

Method: S. Tippner, D. Hernández-Castillo, F. H. Schacher, L. González,
*All-atom molecular dynamics simulations of grafted poly(N,N-dimethylaminoethyl methacrylate) brushes*,
J. Phys. Chem. B **129**, 2105–2114 (2025), https://doi.org/10.1021/acs.jpcb.4c07928.
This is a rewrite and extension of the scripts used there (see *Differences* below).

**Tutorial with images:** [docs/TUTORIAL.md](docs/TUTORIAL.md).

## What you set and what the tool works out

| You set | The tool works out |
|---|---|
| `grafting_density` (chains per nm²) and `chains_per_box` | box side L = √(chains / density) |
| `min_distance` between chains | lowers it (to `min_distance_floor`) if the density does not pack, otherwise stops with a message |
| solvent `layer` above the brush (Å, or `"auto"`) | box height = brush height + layer + 2 Å |
| solvent components with vol% or mole fractions | number of molecules of every component (from the density, see below) |
| counterion | one per charged unit, placed with periodic copies |

Without `grafting_density`, the tool searches the densest packing at `min_distance`, as in the paper.

## Steps

| Step | Subcommand | What it does |
|---|---|---|
| 1 | `params` | head, middle and tail residue of every monomer type from a RESP-charged antechamber file (`prepgen`), GAFF2 parameters (`parmchk2`) |
| 2 | `chains` | random sequences: `round(alpha * n_units)` protonated units at random positions, random stereochemistry; built with `tleap`, aligned along z with the grafted unit at the bottom |
| 3 | `pack` | packing in a region of 3 × 3 cells of side L: each chain is placed with its grafted atom on the plane `anchor_z` inside the central cell, against all chains placed before *and their copies*, then copied into the 8 neighbouring cells (copies kept within half a cell of the central cell). After `max_attempts` failures for one chain the packing restarts with chain 1; after `max_fail` failed packings `min_distance` is lowered by `min_distance_step` (fixed density) or the region enlarged by `region_step` (densest mode) |
| 4 | `solvate` | counterions are placed in the central cell against the chains and their copies and copied into the neighbouring cells like the chains; the central cell is then filled with the solvent around all of them, in two zones: the brush zone (solvent mass = target density × volume − polymer mass, as `solvate.tcl` of PACKMOL, started at `brush_fill`) and the layer above (corrected until its measured density is within 3 % of the target); solvent molecules clashing across the periodic boundaries are removed, the mixture ratio restored, and `tleap` writes the periodic system |
| 5 | `single` | optional: every chain in the same solvent on its own, with its counterions |

`all` runs steps 1 to 4.

## Solvents and mixtures

- **Water** uses the models of AmberTools: `water_model = "tip3p"`, `"spce"`, `"opc"` or `"spcfw"` (flexible SPC/Fw).
- **Other solvents** are given by you: a one-molecule structure (PDB with the residue and atom names of the
  parameters), the residue file (`.prepi` or `.lib`), the `.frcmod`, and the density of the pure liquid.
  Examples: acetonitrile (six-site model) and methanol (GAFF2) in `examples/pdmaema/solvents/`.
- **Composition** by volume (`mode = "volume"`, fractions in vol%, converted with nᵢ ∝ φᵢ ρᵢ / Mᵢ) or by mole
  fraction (`mode = "mole"`). The mixture density is that of ideal mixing of the pure volumes, or `density_mix` if
  you give the measured value; NPT equilibration corrects the small difference.
- The final numbers, mole fractions and vol% of every component are printed and written to `system.json`.

## Requirements

- Python ≥ 3.11 with numpy and scipy
- AmberTools (`prepgen`, `parmchk2`, `tleap`) on `PATH`, e.g. `source $AMBERHOME/amber.sh`
- PACKMOL 20.x (path in `[tools] packmol`, or on `PATH`)

## Usage

```bash
cd examples/pdmaema
python3 ../../brush_builder.py config_quick.toml all      # 4 chains, 0.5 nm^-2, water: a few minutes
python3 ../../brush_builder.py config_mixture.toml all    # same brush in water/acetonitrile 95/5 vol%
python3 ../../brush_builder.py config_methanol.toml all   # same brush in water/methanol 50/50 mol%
python3 ../../brush_builder.py config_paper.toml all      # 16 x 30-mers, densest packing as in the paper (hours)
```

Steps can be run one at a time (`params`, `chains`, `pack`, `solvate`, `single`); each reads the output of the
previous one from the output folder.

## Configuration

| Section | Key | Meaning |
|---|---|---|
| `[monomers]` | `<key> = { ac, protonated }` | monomer types; residues `H<key>`, `D<key>`, `T<key>` |
| `[residues]` | `head_name`, `tail_name`, `end_group_atoms`, `cap_atom`, `pre_head_type`, `post_tail_type` | prepgen definitions from atom names |
| `[chains]` | `n_units`, `alpha`, `seed` | chain length, fraction of protonated units, random seed |
| `[brush]` | `chains_per_box`, `grafting_density` | chains in the box; chains per nm² (omit for the densest mode) |
| | `min_distance`, `min_distance_floor`, `min_distance_step` | minimum distance between chains (Å) and how far it may be lowered |
| | `anchor_atom`, `anchor_z` | grafted atom (index in the head unit) and its plane |
| | `max_attempts`, `max_fail`, `box_step`, `region_step` | restart and enlargement settings |
| `[solvent]` | `water_model`, `mode`, `density_mix` | water model; vol% or mole fractions; measured mixture density (optional) |
| | `layer`, `brush_fill`, `margin` | solvent above the brush (Å or `"auto"`); start fill of the brush zone (0.9); solvent beyond the cell (0) |
| | `counterion`, `shell_single` | counterion name; solvent shell of single chains |
| `[[solvent.components]]` | `name`, `fraction`, (`density`, `structure`, `params`) | one table per component; water needs only name and fraction |

## Output (`<output>/`)

| Path | Content |
|---|---|
| `params/` | `*.prepi`, `*.frcmod`, generated `mainchain.*` |
| `chains/chain_<i>.pdb`, `chains/sequences.json` | aligned chains and their sequences |
| `pack/region.pdb`, `pack/central.pdb`, `pack/pack.json` | chains with copies, central chains, box side, grafting density, final minimum distance |
| `solvate/system.prmtop`, `system.inpcrd`, `system.pdb`, `system.json` | the periodic system; `system.json` holds the box, the restraint mask of the grafted units and the solvent composition |
| `single/chain_<i>/system.prmtop`, ... | single-chain systems |

## Example inputs (`examples/pdmaema`)

- `monomers/`: the two stereoisomers of the DMAEMA unit, unprotonated (`dda.ac`, `ddb.ac`) and protonated
  (`dpa.ac`, `dpb.ac`), with RESP charges from B3LYP-D3(BJ)/def2-SVP Merz–Kollman potentials (Gaussian inputs in `qm/`)
- `solvents/acetonitrile/`, `solvents/methanol/`: structure, residue file and parameters
- `config_quick.toml`, `config_mixture.toml`, `config_methanol.toml`, `config_paper.toml` (see *Usage*)

## Differences to the original scripts

- One tool and one configuration; the 12 `mainchain` files are generated from atom names. Residues are identical to
  the published ones (atoms, types, charges).
- New: grafting density and chains per box as input, box dimensions and solvent amounts derived; solvent mixtures by
  vol% or mole fraction; choice of water model.
- New: the counterions are packed with periodic copies, like the chains, so they also see their periodic neighbours.
  In the original only the chains had copies, and Cl⁻ near the box faces could overlap with polymer or water across
  the boundary (67 Cl⁻–polymer contacts below 2 Å in the published 0.68 nm⁻² start structure).
- The solvent is packed in the central cell only, around the chains, the ions and their copies (the original
  solvated the whole 3 × 3 region and cut it). Solvent molecules that clash across the boundaries are removed; the
  box has a gap of 2 Å at the top.
- Chains are aligned along z by rotating their longest axis onto z (PACKMOL moves molecules as rigid bodies, so the
  original box fitting gave the same rotation with many PACKMOL runs).
- No duplicate zero-shift copy of each chain; copies kept if any atom lies within half a cell of the central cell.

## Checks

With the quick examples (4 chains of 8 units at 0.5 nm⁻², AmberTools 23, PACKMOL 20.4.0): residues identical to
the published ones; box 28.28 × 28.28 × 58.1 Å; total charge 0; no heavy atoms of different molecules closer than
2 Å, periodic boundaries included; in the brush about 0.4 to 0.5 g/cm³ of solvent next to 0.5 g/cm³ of polymer,
about 0.95 g/cm³ in the water layer above; mixtures at 4.8 vol% acetonitrile (target 5) and 50.0 mol% methanol.
Requesting 2.5 nm⁻² for the same chains lowers the minimum distance to 1.5 Å and stops with a message; the densest
mode reproduces the result of version 1.
