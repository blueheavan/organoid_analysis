# Organoid Pipeline — 3D 类器官分析

> **Current local release assessment (2026-09-10): NOT READY FOR THE SPECIFIED RESEARCH USE.** The full audit found a failed surface-accuracy specification and missing project-specific segmentation, threshold and inferential validation. Engineering and controlled-test results are bounded separately in [the validation report](docs/VALIDATION_REPORT.md).

用于 3D 显微镜 Z-stack 的本地类器官分析工具，提供 Cellpose 3D 分割、可审计的多层级测量、结果浏览与表格统计。项目面向 macOS Apple Silicon，使用 Pixi 管理运行环境。

```text
Organoid instance mask
└── Cell instance mask
    └── Nucleus instance mask
```

所有三维定量均基于 `(Z, Y, X)` instance-label mask 及真实 voxel spacing；MIP、单张切片和密度投影仅可用于显示或质控，不能替代三维测量。

> 仅供研究分析与方法开发使用，不用于临床诊断、治疗选择或生物学结论的自动判定。

## 功能

| 模块 | 功能 |
|---|---|
| Cellpose 3D 分割 | 对 nucleus Z-stack 及可选 cell/cytoplasm Z-stack 进行 3D instance segmentation |
| 多层级 3D 分析 | 基于 organoid、cell、nucleus 三套已配准标签，计算层级、形态、空间组织与 QC |
| 经典分析 CLI | 进行类器官形态、marker readout、viability 与统计汇总 |
| Web 应用 | 上传和预览图像、运行或恢复分割、浏览真实 3D 几何与导出结果 |
| 数据导出 | 生成 Parquet feature tables、OME-TIFF label masks、QC 表与 JSON 汇总 |

## 安装

安装 [Pixi](https://pixi.sh/)，然后在项目根目录执行：

```bash
brew install pixi    # 尚未安装 Pixi 时
pixi install
```

建议使用 Pixi 提供的固定环境，不必另外创建 `venv`、`requirements.txt` 或 Conda 环境。首次运行 Cellpose 可能需要联网下载预训练权重；权重会由 Cellpose 缓存。Apple Silicon 会优先使用 MPS，个别不支持的 PyTorch 操作会回退到 CPU。

可先检查环境：

```bash
pixi run check
pixi run analysis-test
pixi run test
```

## 快速开始

### 网页应用

```bash
pixi run web
```

终端会显示本地访问地址（通常为 `http://localhost:8501`）。应用包含以下工作区：

1. **Upload & preview**：上传并检查 nucleus Z-stack；如需 cell mask，也上传已配准的 cell/cytoplasm Z-stack。
2. **Segmentation results**：执行 Cellpose 3D 分割，查看 mask overlay 和对象结果。
3. **Object features**：查看分割 mask 的对象级描述性特征。
4. **3D Analysis results**：恢复保存的 segmentation，或读取完整多层级分析结果。
5. **Statistical analysis (Excel)**：使用现有表格进行探索性统计。

分割结果保存在：

```text
results/segmentation_output/run-*/
```

刷新页面或关闭浏览器后，可在 **3D Analysis results** 使用 **Restore saved Cellpose segmentation** 恢复 nuclei/cell masks。没有 cell mask 的历史分割不能用于多层级分析，因为不能由 nucleus、MIP 或投影图推断 cell labels。

### 多层级 3D 分析

该工作流要求三套已注册、位于同一物理 grid 的 3D instance-label TIFF/OME-TIFF：

| 输入 | 含义 | 要求 |
|---|---|---|
| `organoid_labels` | whole-organoid instances | `ZYX`、背景为 0、正整数 instance ID |
| `cell_labels` | cell/cytoplasm instances | 与 organoid labels 的 shape 和 grid 一致 |
| `nucleus_labels` | nucleus instances | 与 organoid/cell labels 的 shape 和 grid 一致 |
| `nucleus_intensity`（可选） | 原始核荧光通道 | 同一 grid；仅用于强度 readout |

若 TIFF 含有一致且完整的 OME physical-spacing metadata：

```bash
pixi run python -m organoid_analysis analyze-3d \
  --organoid-labels /absolute/path/organoid_labels.ome.tif \
  --cell-labels /absolute/path/cell_labels.ome.tif \
  --nucleus-labels /absolute/path/nucleus_labels.ome.tif \
  --sample-id sample_001 --well-id A01 --field-id 1 \
  --out /absolute/path/Result
```

没有完整 OME metadata 时，必须显式提供真实采集 spacing（单位为 µm）：

```bash
pixi run python -m organoid_analysis analyze-3d \
  --organoid-labels /absolute/path/organoid_labels.tif \
  --cell-labels /absolute/path/cell_labels.tif \
  --nucleus-labels /absolute/path/nucleus_labels.tif \
  --nucleus-intensity /absolute/path/hoechst_raw.tif \
  --spacing-z-um 2.0 --spacing-y-um 0.65 --spacing-x-um 0.65 \
  --sample-id sample_001 --well-id A01 --field-id 1 \
  --batch-id batch_01 --condition control --timepoint day_7 \
  --out /absolute/path/Result
```

运行规则：

- 三个 spacing 参数要么全部提供，要么全部省略。
- 显式 spacing 与 OME metadata 不一致时，命令会报错。
- `--out` 必须是不存在或空目录。
- TIFF axes metadata 缺失或错误、但文件已人工确认时，可使用 `--axes ZYX`。

查看全部选项：

```bash
pixi run python -m organoid_analysis analyze-3d --help
```

### 查看 Result

在 **3D Analysis results** 的 **Result directory** 中填写 `Result` 的绝对路径。页面只读取已有导出，支持：

- Organoid、Cell、Nucleus、Topology 与 QC 表的筛选、排序和下载；
- 三个标签层级的真实 3D surface geometry；
- Z-slice 质控，包括 Mask、Boundary、Object ID 和 QC flags；
- 对象的 volume、surface area、sphericity、parent ID 与 QC status。

## 测量与 QC

### 层级分配

Parent assignment 使用 maximum voxel overlap：

```text
cell_id     → organoid_id
nucleus_id  → cell_id → organoid_id
```

Nucleus 的 organoid membership 从其已分配的 cell 继承。系统同时保留直接 nucleus-to-organoid overlap 作为审计字段；两者不一致时写入 `direct_organoid_parent_mismatch` QC。

### 三维形态与空间组织

每个 organoid、cell 和 nucleus 可输出 voxel count、体积、物理坐标质心、surface area、sphericity、主轴、elongation、prolate/oblate ratios 和 surface-to-volume ratio。

- Volume、surface area 与 sphericity 均来自原始 label mask。
- Surface area 使用 native physical spacing 的 marching cubes 估计。
- Sphericity 为 `π^(1/3) × (6V)^(2/3) / A`。
- Cell topology 仅把共享 voxel face 的两个 cell 视为直接接触；contact area 会按各向异性 voxel face area 计算。
- 距离与 radial position 使用 physical spacing；`normalized_radial_position_equivalent_radius` 是等效球半径归一化值，不是 local-radius estimate。

Cell 特征包含 nucleus count、anucleate/multinucleated 状态、nucleus-to-cell volume ratio 等。只有提供 `--nucleus-intensity` 时才导出核强度 readout；缺少原始强度通道时不会生成推测值。

QC 是正式输出的一部分，涵盖 border touch、fragmentation、对象过小、parent assignment、低 overlap、跨层级归属、anucleate/multinucleated 与 volume outlier 等情形。Volume outlier 默认采用 robust MAD。

## 输出目录

```text
Result/
├── features/
│   ├── organoid_features.parquet
│   ├── cell_features.parquet
│   ├── nucleus_features.parquet
│   └── cell_topology_edges.parquet
├── qc/
│   └── qc_flags.parquet
├── summary/
│   └── analysis_summary.json
└── masks/
    ├── organoid_labels.ome.tif
    ├── cell_labels.ome.tif
    └── nucleus_labels.ome.tif
```

Parquet 是主要 feature database，可用 pandas、Polars、DuckDB 或 R/Arrow 读取。命令行提供的 `sample_id`、`well_id`、`field_id`、`batch_id`、`condition`、`treatment`、`dose` 和 `timepoint` 会保留在 feature/QC 表中。

## 其他命令行工作流

### 类器官形态、marker 与 viability 分析

```bash
pixi run analysis-analyze \
  --manifest /absolute/path/samples.csv \
  --config configs/structural_fluorescence.yaml \
  --out /absolute/path/organoid_results
```

可选配置：

| 配置 | 输入与用途 |
|---|---|
| `structural_fluorescence.yaml` | structural fluorescence + Calcein/PI |
| `brightfield_probability.yaml` | brightfield + 已配准 foreground probability map |
| `imported_instances.yaml` | 导入已验证的 3D instance labels |
| `brightfield_exploratory.yaml` | 明场 threshold/watershed 的探索性起点 |

### 已注册 cell–nucleus 配对

```bash
pixi run analysis-cells \
  --cell-labels /absolute/path/cell_masks.tif \
  --nucleus-labels /absolute/path/nucleus_masks.tif \
  --spacing-z-um 2.0 --spacing-y-um 0.5 --spacing-x-um 0.5 \
  --out /absolute/path/cell_results
```

此工作流提供 cell/nucleus 配对特征；若需要完整 organoid → cell → nucleus hierarchy、3D topology 或空间组织，请使用 `analyze-3d`。

### 合成示例

```bash
pixi run analysis-demo
open data/synthetic_demo/results/report.html
```

合成数据用于软件测试和演示，不能替代真实生物学验证。

## 常见问题

**“3D Analysis results 是空的”**

该页面读取已完成的 `Result/` 目录。运行 `analyze-3d`，或在网页恢复 segmentation 后提供已配准的 organoid label 并执行多层级分析，再填写生成的目录。

**“Voxel spacing is missing”**

TIFF 缺少完整物理 metadata。请提供：

```bash
--spacing-z-um Z --spacing-y-um Y --spacing-x-um X
```

这些值应来自显微镜采集记录或可信 metadata；不要以 `1.0` 作为最终定量的占位值。

**“Input TIFF voxel spacings differ” 或 shape 不一致**

所有 labels 和可选 intensity 必须位于同一 `ZYX` grid。请先完成上游 registration/resampling；分析不会静默裁剪、拉伸或重采样 instance labels。

## 项目结构与测试

```text
src/organoid_analysis/
├── command_line/           # CLI 入口（organoid_commands.py 等）
├── config.py                # 顶层运行配置
├── microscopy_io/          # TIFF/OME 读取、体素间距 (spacing) 解析与统一契约
├── segmentation/            # Cellpose 3D 与经典 (watershed) 分割
├── quantification/          # per-object mask 特征、multilevel 关系
├── phenotyping/             # 细胞/器官表型分类
├── statistics/               # 混合效应/聚类稳健推断
├── validation/               # 科学验证辅助
├── workflows/                # Web/CLI 共用的读取→分割→测量→导出编排
├── result_export/            # 结果导出 (CSV/OME-TIFF 等)
├── visualization/             # 3D 体渲染查看器
└── web_interface/            # Streamlit 应用
configs/                   # YAML 分析配置
data/                      # 示例与本地数据
docs/                      # 方法、参数与验证说明
results/                   # 本地分割和分析输出
tests/                     # 单元、集成与 UI 测试，镜像上述包结构
```

```bash
pixi run lint
pixi run test
pixi run test-render
pixi run check
pixi run typecheck
pixi run python -m organoid_analysis --help
```

方法、参数、输入域、验证范围与已知限制见 [docs](docs/) 下的科学文档。真实数据上的生物学精度与参数敏感性应针对实际样本、成像条件和研究问题独立评估。

## License

本项目采用 [MIT License](LICENSE)。本地化基础工作参考 [bnsreenu/3D-Organoid-Analysis](https://github.com/bnsreenu/3D-Organoid-Analysis)。
