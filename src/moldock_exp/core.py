from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class BoxConfig:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


@dataclass
class DockingParams:
    exhaustiveness: int
    num_modes: int
    energy_range: int
    seed: int
    cpu: int


@dataclass
class PoseResult:
    mode_rank: int
    affinity_kcal_mol: float
    heavy_atom_rmsd_A: float
    is_success_2A: int
    is_success_3A: int


STANDARD_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HID",
    "HIE",
    "HIP",
    "CYX",
    "ASH",
    "GLH",
}


def make_scratch(work_root: str | None, name: str) -> Path:
    scratch = Path(work_root).resolve() if work_root else Path(tempfile.gettempdir()) / name
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    return scratch


def env_scripts_dir() -> Path:
    exe = Path(sys.executable).resolve()
    return exe.parent / "Scripts" if os.name == "nt" else exe.parent


def find_tool(name: str, *, python_script: bool = False) -> Path:
    candidates: list[Path] = []
    found = shutil.which(name)
    if found:
        candidates.append(Path(found))
    scripts = env_scripts_dir()
    candidates.append(scripts / name)
    if os.name == "nt" and not python_script and not name.endswith(".exe"):
        candidates.append(scripts / f"{name}.exe")
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Cannot find {name}. Activate the conda environment containing Vina and Meeko.")


def run_command(command: list[str], cwd: Path) -> None:
    print("  " + " ".join(f'"{x}"' if " " in x else x for x in command))
    subprocess.run(command, cwd=str(cwd), check=True)


def parse_box_config(path: Path) -> BoxConfig:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, float] = {}
    for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
        match = re.search(rf"--{key}\s+(-?\d+(?:\.\d+)?)", text)
        if not match:
            raise ValueError(f"Missing --{key} in {path}")
        values[key] = float(match.group(1))
    return BoxConfig(**values)


def resolve_path(data_dir: Path, value: str | None, default_name: str) -> Path:
    path = Path(value) if value else data_dir / default_name
    return path if path.is_absolute() else data_dir / path


def resolve_protein(data_dir: Path, value: str | None) -> Path:
    if value:
        return resolve_path(data_dir, value, value)
    protein_with_h = data_dir / "protein_with_h.pdb"
    return protein_with_h if protein_with_h.exists() else data_dir / "protein.pdb"


def copy_inputs_to_scratch(
    data_dir: Path,
    scratch: Path,
    *,
    protein: str | None,
    ligand: str | None,
    box_config: str | None,
) -> tuple[Path, Path, Path]:
    protein_src = resolve_protein(data_dir, protein)
    ligand_src = resolve_path(data_dir, ligand, "ligand.sdf")
    box_src = resolve_path(data_dir, box_config, "box_config.txt")
    for path in (protein_src, ligand_src, box_src):
        if not path.exists():
            raise FileNotFoundError(path)
    protein_dst = scratch / protein_src.name
    ligand_dst = scratch / "ligand.sdf"
    box_dst = scratch / "box_config.txt"
    shutil.copy2(protein_src, protein_dst)
    shutil.copy2(ligand_src, ligand_dst)
    shutil.copy2(box_src, box_dst)
    return protein_dst, ligand_dst, box_dst


def clean_receptor_to_standard_residues(protein: Path, scratch: Path) -> Path:
    cleaned = scratch / "receptor_standard_residues.pdb"
    kept_atoms = 0
    skipped_atoms = 0
    with protein.open(encoding="utf-8", errors="ignore") as source, cleaned.open("w", encoding="ascii") as target:
        for line in source:
            record = line[:6].strip()
            if record in {"ATOM", "HETATM"}:
                residue_name = line[17:20].strip()
                if record == "ATOM" or residue_name in STANDARD_RESIDUES:
                    target.write(line)
                    kept_atoms += 1
                else:
                    skipped_atoms += 1
            elif record in {"TER", "END"}:
                target.write(line)
    print(f"\nReceptor cleanup: kept {kept_atoms} protein atoms, skipped {skipped_atoms} non-standard atoms")
    return cleaned


def prepare_receptor(
    protein: Path,
    scratch: Path,
    *,
    default_altloc: str,
    allow_bad_res: bool,
    clean_receptor: bool,
) -> Path:
    receptor_input = clean_receptor_to_standard_residues(protein, scratch) if clean_receptor else protein
    command = [
        sys.executable,
        str(find_tool("mk_prepare_receptor.py", python_script=True)),
        "--read_pdb",
        str(receptor_input),
        "-o",
        str(scratch / "receptor"),
        "-p",
        "--default_altloc",
        default_altloc,
    ]
    if allow_bad_res:
        command.append("--allow_bad_res")
    print("\nPreparing receptor")
    run_command(command, scratch)
    return scratch / "receptor.pdbqt"


def ensure_explicit_hydrogens(ligand_sdf: Path, scratch: Path) -> Path:
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
    molecules = [mol for mol in supplier if mol is not None]
    if not molecules:
        raise ValueError(f"RDKit could not read any molecules from {ligand_sdf}")

    explicit_h_sdf = scratch / "ligand_explicit_h.sdf"
    writer = Chem.SDWriter(str(explicit_h_sdf))
    for mol in molecules:
        writer.write(Chem.AddHs(mol, addCoords=True))
    writer.close()
    print(f"\nLigand cleanup: wrote explicit hydrogens for {len(molecules)} molecule(s)")
    return explicit_h_sdf


def prepare_ligand(ligand_sdf: Path, scratch: Path) -> Path:
    ligand_input = ensure_explicit_hydrogens(ligand_sdf, scratch)
    ligand_pdbqt = scratch / "ligand.pdbqt"
    command = [
        sys.executable,
        str(find_tool("mk_prepare_ligand.py", python_script=True)),
        "-i",
        str(ligand_input),
        "-o",
        str(ligand_pdbqt),
    ]
    print("\nPreparing ligand")
    run_command(command, scratch)
    return ligand_pdbqt


def run_vina(receptor: Path, ligand: Path, box: BoxConfig, out_dir: Path, params: DockingParams) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    docked_pdbqt = out_dir / "docked_poses.pdbqt"
    log_file = out_dir / "vina.log"
    command = [
        str(find_tool("vina")),
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(box.center_x),
        "--center_y",
        str(box.center_y),
        "--center_z",
        str(box.center_z),
        "--size_x",
        str(box.size_x),
        "--size_y",
        str(box.size_y),
        "--size_z",
        str(box.size_z),
        "--exhaustiveness",
        str(params.exhaustiveness),
        "--num_modes",
        str(params.num_modes),
        "--energy_range",
        str(params.energy_range),
        "--seed",
        str(params.seed),
        "--cpu",
        str(params.cpu),
        "--out",
        str(docked_pdbqt),
        "--log",
        str(log_file),
    ]
    print(f"\nRunning Vina: exhaustiveness={params.exhaustiveness}, seed={params.seed}, num_modes={params.num_modes}")
    run_command(command, out_dir)
    return docked_pdbqt, log_file


def export_poses(docked_pdbqt: Path, out_dir: Path) -> Path:
    docked_sdf = out_dir / "docked_poses.sdf"
    command = [
        sys.executable,
        str(find_tool("mk_export.py", python_script=True)),
        str(docked_pdbqt),
        "-s",
        str(docked_sdf),
    ]
    print("\nExporting docked poses")
    run_command(command, out_dir)
    return docked_sdf


def parse_affinities(pdbqt_path: Path) -> list[float]:
    pattern = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")
    values: list[float] = []
    with pdbqt_path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                values.append(float(match.group(1)))
    return values


def _remove_hs(mol):
    from rdkit import Chem

    return Chem.RemoveHs(mol, sanitize=True)


def direct_symmetric_heavy_rmsd(reference, pose) -> float:
    ref = _remove_hs(reference)
    dock = _remove_hs(pose)
    if ref.GetNumAtoms() != dock.GetNumAtoms():
        raise ValueError(f"Heavy atom count mismatch: reference={ref.GetNumAtoms()}, pose={dock.GetNumAtoms()}")
    ref_conf = ref.GetConformer()
    dock_conf = dock.GetConformer()
    matches = dock.GetSubstructMatches(ref, uniquify=False, maxMatches=10000)
    if not matches:
        raise ValueError("Could not match docked pose graph to reference ligand graph")

    best = math.inf
    for match in matches:
        total = 0.0
        for ref_idx, dock_idx in enumerate(match):
            rp = ref_conf.GetAtomPosition(ref_idx)
            dp = dock_conf.GetAtomPosition(dock_idx)
            total += (rp.x - dp.x) ** 2 + (rp.y - dp.y) ** 2 + (rp.z - dp.z) ** 2
        best = min(best, math.sqrt(total / ref.GetNumAtoms()))
    return best


def evaluate_poses(reference_sdf: Path, docked_sdf: Path, docked_pdbqt: Path) -> list[PoseResult]:
    from rdkit import Chem

    reference = Chem.SDMolSupplier(str(reference_sdf), removeHs=False)[0]
    if reference is None:
        raise ValueError(f"RDKit could not read reference ligand: {reference_sdf}")
    affinities = parse_affinities(docked_pdbqt)
    results: list[PoseResult] = []
    for idx, pose in enumerate(Chem.SDMolSupplier(str(docked_sdf), removeHs=False), start=1):
        if pose is None:
            continue
        rmsd = direct_symmetric_heavy_rmsd(reference, pose)
        affinity = affinities[idx - 1] if idx - 1 < len(affinities) else float("nan")
        results.append(PoseResult(idx, affinity, rmsd, int(rmsd <= 2.0), int(rmsd <= 3.0)))
    return results


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def pose_results_to_rows(results: list[PoseResult], extra: dict[str, object] | None = None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        row = asdict(result)
        if extra:
            row = {**extra, **row}
        rows.append(row)
    return rows


def summarize_pose_results(results: list[PoseResult], extra: dict[str, object] | None = None) -> dict[str, object]:
    best_affinity = min(results, key=lambda r: r.affinity_kcal_mol)
    best_rmsd = min(results, key=lambda r: r.heavy_atom_rmsd_A)
    row: dict[str, object] = {
        "poses_returned": len(results),
        "best_affinity_kcal_mol": best_affinity.affinity_kcal_mol,
        "best_affinity_mode": best_affinity.mode_rank,
        "rmsd_of_best_affinity_A": best_affinity.heavy_atom_rmsd_A,
        "best_rmsd_A": best_rmsd.heavy_atom_rmsd_A,
        "best_rmsd_mode": best_rmsd.mode_rank,
        "affinity_at_best_rmsd_kcal_mol": best_rmsd.affinity_kcal_mol,
        "success_2A_count": sum(r.is_success_2A for r in results),
        "success_3A_count": sum(r.is_success_3A for r in results),
    }
    return {**extra, **row} if extra else row


def summarize_ablation_results(results: list[PoseResult], extra: dict[str, object] | None = None) -> dict[str, object]:
    best_affinity = min(results, key=lambda r: r.affinity_kcal_mol)
    best_rmsd = min(results, key=lambda r: r.heavy_atom_rmsd_A)
    top5_by_affinity = sorted(results, key=lambda r: r.affinity_kcal_mol)[:5]
    average_top5_rmsd = sum(r.heavy_atom_rmsd_A for r in top5_by_affinity) / len(top5_by_affinity)
    row: dict[str, object] = {
        "rmsd_of_best_affinity_A": best_affinity.heavy_atom_rmsd_A,
        "best_rmsd_A": best_rmsd.heavy_atom_rmsd_A,
        "average_rmsd_of_top5_affinity_A": average_top5_rmsd,
    }
    return {**extra, **row} if extra else row


def copy_run_files(run_dir: Path, output_dir: Path, include_prepared: tuple[Path, Path] | None = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("docked_poses.pdbqt", "docked_poses.sdf", "vina.log"):
        source = run_dir / name
        if source.exists():
            shutil.copy2(source, output_dir / name)
    if include_prepared:
        receptor, ligand = include_prepared
        shutil.copy2(receptor, output_dir / "receptor.pdbqt")
        shutil.copy2(ligand, output_dir / "ligand.pdbqt")


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
