from __future__ import annotations

import argparse
import os
from dataclasses import asdict
from pathlib import Path

from .core import (
    DockingParams,
    copy_inputs_to_scratch,
    copy_run_files,
    evaluate_poses,
    export_poses,
    make_scratch,
    parse_box_config,
    pose_results_to_rows,
    prepare_ligand,
    prepare_receptor,
    run_vina,
    summarize_pose_results,
    write_csv,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run one CPU AutoDock Vina docking job and evaluate ligand RMSD.")
    parser.add_argument("--data-dir", default=".", help="Directory containing protein, ligand and box files. Default: current directory.")
    parser.add_argument("--protein", help="Protein filename/path. Defaults to protein_with_h.pdb, then protein.pdb.")
    parser.add_argument("--ligand", help="Ligand SDF filename/path. Defaults to ligand.sdf.")
    parser.add_argument("--box-config", help="Box config filename/path. Defaults to box_config.txt.")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--num-modes", type=int, default=10)
    parser.add_argument("--energy-range", type=int, default=4)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--default-altloc", default="A")
    parser.add_argument("--no-allow-bad-res", action="store_true")
    parser.add_argument("--output-name", default="single_docking")
    parser.add_argument("--work-root", help="ASCII scratch directory. Defaults to the system temp directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    data_dir = Path(args.data_dir).resolve()
    output_dir = data_dir / "vina_runs" / args.output_name
    scratch = make_scratch(args.work_root, "vina_single_docking")

    protein, ligand_sdf, box_file = copy_inputs_to_scratch(
        data_dir,
        scratch,
        protein=args.protein,
        ligand=args.ligand,
        box_config=args.box_config,
    )
    box = parse_box_config(box_file)
    receptor = prepare_receptor(
        protein,
        scratch,
        default_altloc=args.default_altloc,
        allow_bad_res=not args.no_allow_bad_res,
    )
    ligand_pdbqt = prepare_ligand(ligand_sdf, scratch)

    params = DockingParams(args.exhaustiveness, args.num_modes, args.energy_range, args.seed, args.cpu)
    run_dir = scratch / "run"
    docked_pdbqt, _ = run_vina(receptor, ligand_pdbqt, box, run_dir, params)
    docked_sdf = export_poses(docked_pdbqt, run_dir)
    results = evaluate_poses(ligand_sdf, docked_sdf, docked_pdbqt)

    copy_run_files(run_dir, output_dir, include_prepared=(receptor, ligand_pdbqt))
    write_csv(output_dir / "pose_rmsd.csv", pose_results_to_rows(results))
    summary = summarize_pose_results(results)
    metadata = {
        "protein": protein.name,
        "ligand": "ligand.sdf",
        "box": asdict(box),
        "params": asdict(params),
        "summary": summary,
    }
    write_json(output_dir / "run_metadata.json", metadata)

    print("\nDone")
    print(f"  output: {output_dir}")
    print(f"  best affinity: {summary['best_affinity_kcal_mol']} kcal/mol")
    print(f"  best RMSD: {summary['best_rmsd_A']:.3f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
