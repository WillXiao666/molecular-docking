# Molecular Docking Experiment

AutoDock Vina workflows for a drug-design molecular docking experiment. The package runs CPU docking for a known binding pocket and evaluates docked poses against the experimental ligand coordinates with heavy-atom RMSD.

## Features

- Single-run Vina docking with receptor/ligand preparation
- Heavy-atom RMSD evaluation against `ligand.sdf`
- Exhaustiveness ablation with repeated random seeds
- Num modes ablation
- Windows-friendly scratch workflow for non-ASCII data paths

## Installation

Create or activate an environment with AutoDock Vina, Meeko, and RDKit:

```powershell
conda activate docking-vina
```

Install this package:

```powershell
pip install -e . --no-deps
```

After installation, five commands are available:

```text
dock-all
dock-run
dock-exhaustiveness-ablation
dock-seed-ablation
dock-num-modes-ablation
```

## Data Layout

Place input files in one data directory:

```text
data/
  protein_with_h.pdb   # preferred receptor; non-standard residues are removed by default
  protein.pdb          # fallback receptor
  ligand.sdf           # experimental ligand coordinates; explicit H atoms are added automatically
  box_config.txt       # Vina pocket box
```

Example `box_config.txt`:

```text
--center_x 13.3 --center_y 3.0 --center_z 0.2 --size_x 18.4 --size_y 16.4 --size_z 19.0
```

## Quick Start

Run from the data directory:

```powershell
dock-all
```

Or specify a data directory:

```powershell
dock-all --data-dir path\to\data
```

This runs all experiments in sequence:

```text
data/vina_runs/results/
  single_docking/
  exhaustiveness_ablation/
  seed_ablation/
  num_modes_ablation/
```

To run only the main docking experiment:

```powershell
dock-run
```

Recommended single-run setting:

```powershell
dock-run --exhaustiveness 16 --num-modes 10 --energy-range 4 --seed 12345 --cpu 20
```

Outputs are written to:

```text
data/vina_runs/single_docking/
```

## Experiments

### Exhaustiveness Ablation

```powershell
dock-exhaustiveness-ablation `
  --exhaustiveness-values 4,8,12,16,24,32,48,64 `
  --seed 12345 `
  --num-modes 10 `
  --energy-range 4 `
  --cpu 20
```

Outputs:

```text
data/vina_runs/exhaustiveness_ablation/
  summary_by_exhaustiveness.csv
```

### Seed Ablation

```powershell
dock-seed-ablation `
  --seeds 1,2,3,4,5,6,7,8,9,10 `
  --exhaustiveness 16 `
  --num-modes 10 `
  --energy-range 4 `
  --cpu 20
```

Outputs:

```text
data/vina_runs/seed_ablation/
  summary_by_seed.csv
```

### Num Modes Ablation

```powershell
dock-num-modes-ablation `
  --num-modes-values 10,20,30,40 `
  --exhaustiveness 16 `
  --energy-range 6 `
  --seed 12345 `
  --cpu 20
```

Outputs:

```text
data/vina_runs/num_modes_ablation/
  summary_by_num_modes.csv
```

## Metrics

Ablation summaries report three RMSD metrics:

```text
rmsd_of_best_affinity_A
best_rmsd_A
average_rmsd_of_top5_affinity_A
```

RMSD is computed on heavy atoms in the receptor coordinate frame without ligand superposition. Symmetry-equivalent atom mappings are considered.

Typical interpretation:

| RMSD | Interpretation |
|---:|---|
| <= 2 A | Successful docking |
| 2-3 A | Acceptable, inspect interactions |
| > 3 A | Poor pose recovery |

## Notes

- `protein_with_h.pdb` is used by default when available.
- Non-standard receptor residues such as glycans, salts, solvent, and co-crystallized ligands are removed by default before receptor preparation.
- Use `--keep-nonstandard-residues` to disable receptor cleanup.
- `ligand.sdf` is automatically converted to an explicit-H SDF before Meeko ligand preparation.
- Receptor preparation uses `default_altloc=A` unless changed by `--default-altloc`.
- Use `--no-allow-bad-res` to disable Meeko's permissive residue handling.
