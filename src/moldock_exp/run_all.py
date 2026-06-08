from __future__ import annotations

import argparse
import os
from pathlib import Path

from . import docking_main, exhaustiveness_ablation, num_modes_ablation, seed_ablation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run all molecular docking experiments in sequence.")
    parser.add_argument("--data-dir", default=".", help="Directory containing protein, ligand and box files. Default: current directory.")
    parser.add_argument("--protein", help="Protein filename/path. Defaults to protein_with_h.pdb, then protein.pdb.")
    parser.add_argument("--ligand", help="Ligand SDF filename/path. Defaults to ligand.sdf.")
    parser.add_argument("--box-config", help="Box config filename/path. Defaults to box_config.txt.")
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--default-altloc", default="A")
    parser.add_argument("--no-allow-bad-res", action="store_true")
    parser.add_argument("--output-root", default="results", help="Folder name under data/vina_runs.")
    parser.add_argument("--work-root", help="Optional ASCII scratch root. Defaults to the system temp directory.")
    return parser


def optional_args(args: argparse.Namespace) -> list[str]:
    values: list[str] = []
    for option, value in (
        ("--protein", args.protein),
        ("--ligand", args.ligand),
        ("--box-config", args.box_config),
        ("--work-root", args.work_root),
    ):
        if value:
            values.extend([option, str(value)])
    if args.no_allow_bad_res:
        values.append("--no-allow-bad-res")
    return values


def run_step(name: str, main_func, argv: list[str]) -> None:
    print(f"\n{'=' * 72}")
    print(name)
    print(f"{'=' * 72}")
    exit_code = main_func(argv)
    if exit_code:
        raise SystemExit(exit_code)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    base = [
        "--data-dir",
        str(Path(args.data_dir).resolve()),
        "--cpu",
        str(args.cpu),
        "--default-altloc",
        args.default_altloc,
        *optional_args(args),
    ]

    run_step(
        "1/4 Single docking",
        docking_main.main,
        [*base, "--output-name", f"{args.output_root}/single_docking"],
    )
    run_step(
        "2/4 Exhaustiveness ablation",
        exhaustiveness_ablation.main,
        [*base, "--output-name", f"{args.output_root}/exhaustiveness_ablation"],
    )
    run_step(
        "3/4 Seed ablation",
        seed_ablation.main,
        [*base, "--output-name", f"{args.output_root}/seed_ablation"],
    )
    run_step(
        "4/4 Num modes ablation",
        num_modes_ablation.main,
        [*base, "--output-name", f"{args.output_root}/num_modes_ablation"],
    )

    results_dir = Path(args.data_dir).resolve() / "vina_runs" / args.output_root
    print("\nAll experiments finished")
    print(f"  results: {results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
