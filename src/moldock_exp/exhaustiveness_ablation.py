from __future__ import annotations

import argparse
import os
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
)


def parse_int_list(text: str) -> list[int]:
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run exhaustiveness-only ablation with a fixed random seed.")
    parser.add_argument("--data-dir", default=".", help="Directory containing protein, ligand and box files. Default: current directory.")
    parser.add_argument("--protein", help="Protein filename/path. Defaults to protein_with_h.pdb, then protein.pdb.")
    parser.add_argument("--ligand", help="Ligand SDF filename/path. Defaults to ligand.sdf.")
    parser.add_argument("--box-config", help="Box config filename/path. Defaults to box_config.txt.")
    parser.add_argument("--exhaustiveness-values", default="4,8,12,16,24,32,48,64")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--num-modes", type=int, default=10)
    parser.add_argument("--energy-range", type=int, default=4)
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--default-altloc", default="A")
    parser.add_argument("--no-allow-bad-res", action="store_true")
    parser.add_argument("--output-name", default="exhaustiveness_ablation")
    parser.add_argument("--work-root", help="ASCII scratch directory. Defaults to the system temp directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    data_dir = Path(args.data_dir).resolve()
    output_dir = data_dir / "vina_runs" / args.output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch = make_scratch(args.work_root, "vina_exhaustiveness_ablation")

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

    summary_rows: list[dict[str, object]] = []
    pose_rows: list[dict[str, object]] = []
    for exhaustiveness in parse_int_list(args.exhaustiveness_values):
        params = DockingParams(exhaustiveness, args.num_modes, args.energy_range, args.seed, args.cpu)
        label = f"exh_{exhaustiveness}"
        run_dir = scratch / label
        docked_pdbqt, _ = run_vina(receptor, ligand_pdbqt, box, run_dir, params)
        docked_sdf = export_poses(docked_pdbqt, run_dir)
        results = evaluate_poses(ligand_sdf, docked_sdf, docked_pdbqt)

        extra = {"exhaustiveness": exhaustiveness, "seed": args.seed}
        summary_rows.append(summarize_pose_results(results, extra))
        pose_rows.extend(pose_results_to_rows(results, extra))
        copy_run_files(run_dir, output_dir / label)

    write_csv(output_dir / "summary_by_exhaustiveness.csv", summary_rows)
    write_csv(output_dir / "pose_rmsd_by_mode.csv", pose_rows)
    print("\nDone")
    print(f"  output: {output_dir}")
    print(f"  summary: {output_dir / 'summary_by_exhaustiveness.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
