# 药物设计学-分子对接实验

本项目整理了 AutoDock Vina 分子对接实验代码，用于在已知蛋白口袋坐标的条件下完成 CPU 对接，并用实验配体坐标评估对接 pose 的重原子 RMSD。

## 功能

- 单次 CPU 对接：自动准备 receptor/ligand，运行 Vina，输出 pose、打分和 RMSD。
- `exhaustiveness × seed` 稳定性验证：比较不同搜索强度和随机种子下结果是否稳定。
- `num_modes` 消融：比较输出构象数量对结果的影响。
- 自动规避 Windows 中文路径兼容问题：运行时会把输入复制到英文临时目录，对接完成后把结果写回 `data/vina_runs`。

## 依赖项

推荐使用已有 Conda 环境：

```powershell
conda activate docking-vina
```

环境中需要：

- Python >= 3.10
- AutoDock Vina，可在命令行调用 `vina`
- Meeko，提供 `mk_prepare_receptor.py`、`mk_prepare_ligand.py`、`mk_export.py`
- RDKit

本项目本身安装方式：

```powershell
cd code
pip install -e .
```

安装后会提供三个命令：

```text
dock-run
dock-exhaustiveness-ablation
dock-num-modes-ablation
```

## 数据目录

默认输入文件放在同一个 `data` 目录下：

```text
data/
  protein_with_h.pdb   推荐，带氢蛋白结构
  protein.pdb          如果没有 protein_with_h.pdb，则使用这个
  ligand.sdf           实验配体坐标，用于 RMSD 评估
  box_config.txt       Vina 口袋坐标和盒子大小
  code/                本项目代码
```

`box_config.txt` 格式示例：

```text
--center_x 13.3 --center_y 3.0 --center_z 0.2 --size_x 18.4 --size_y 16.4 --size_z 19.0
```

## 使用方法

在 `data` 目录运行命令最方便：

```powershell
cd ..
dock-run
```

也可以从任意目录指定数据目录：

```powershell
dock-run --data-dir "D:\path\to\data"
```

常用单次对接参数：

```powershell
dock-run --exhaustiveness 16 --num-modes 10 --energy-range 4 --seed 12345 --cpu 20
```

输出目录：

```text
data/vina_runs/single_docking/
  receptor.pdbqt
  ligand.pdbqt
  docked_poses.pdbqt
  docked_poses.sdf
  vina.log
  pose_rmsd.csv
  run_metadata.json
```

## Exhaustiveness 消融 + Seed 重复验证

推荐用于检查结果稳定性：

```powershell
dock-exhaustiveness-ablation --exhaustiveness-values 8,16,32,64 --seeds 1,2,3,4,5 --num-modes 10 --energy-range 4 --cpu 20
```

输出目录：

```text
data/vina_runs/exhaustiveness_seed_ablation/
  summary_by_run.csv
  pose_rmsd_by_mode.csv
  exh_8_seed_1/
  exh_8_seed_2/
  ...
```

重点查看：

- `rmsd_of_best_affinity_A`：Vina 最佳打分 pose 与实验配体的 RMSD
- `best_rmsd_A`：所有输出 pose 中最接近实验构象的 RMSD
- `success_2A_count`：RMSD <= 2 A 的 pose 数量

## Num Modes 消融

用于比较输出 pose 数量的影响：

```powershell
dock-num-modes-ablation --num-modes-values 10,20,30,40 --exhaustiveness 16 --energy-range 4 --seed 12345 --cpu 20
```

如果 `num_modes` 增大但实际输出 pose 数没有增加，通常是因为 `energy_range` 限制了输出能量窗口。

## 结果判读

RMSD 采用重原子、同一受体坐标系下、不对配体重新叠合的计算方式，并考虑对称原子匹配。

一般判读：

- RMSD <= 2 A：对接成功
- RMSD 2-3 A：基本可接受，需要结合相互作用检查
- RMSD > 3 A：pose 偏差较大

本实验中 `ligand.sdf` 是实验配体坐标，因此 RMSD 可以直接用于评估对接构象复现能力。
