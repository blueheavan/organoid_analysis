# 当前代码科学性审核 — 2026-09-11

**总体结论：NOT READY FOR THE SPECIFIED RESEARCH USE。**

适用范围：当前版本从三维显微图像/实例标签生成物理测量、层级关系、活性信号代理及探索性统计结果的流程。当前证据不支持把这些输出整体作为已验证的定量研究结果。发现两个 P1 问题及四个 P2 问题/默认配置风险；表面积精度失败和元数据冲突漏检阻止正面的研究就绪结论。这里不声称每个输入或每个测量均错误。

本次只审核：未修改生产代码、测试套件、科学算法、默认阈值或验收标准，未 commit/push。新增本报告及 [本次证据目录](evidence/2026-09-11-code-science-audit/)。既有报告仅用于理解声明与验收要求；以下缺陷来自当前源码和新执行的反例。不是把旧 PASS/FAIL 当作本次证明。

## 1. 审核对象、预设标准与证据边界

- 基线提交：`3639bec910337ec623782b657fe62cf1dafc3800`，初始 `main` 工作区干净。
- 环境：锁定 Pixi；Python 3.12.14、NumPy 2.5.2、SciPy 1.17.1、pandas 2.3.3、scikit-image 0.26.0、scikit-learn 1.9.0、statsmodels 0.14.6、Cellpose 4.2.1.1、tifffile 2026.3.3；macOS ARM64。详见 [snapshot.json](evidence/2026-09-11-code-science-audit/snapshot.json) 与 [probe.json](evidence/2026-09-11-code-science-audit/probe.json)。
- 审核路径：Audit-only，最高影响 S3。`config.py` 默认启用统计；`statistics/inference.py` 产生模型检验、置信区间和 FDR 判定。I/O 尺度与定量测量为 S2；探索性分类/聚类也是 S3。无临床决策用途，S4 评估 NOT APPLICABLE。
- 科学单位：一个 organoid/cell/nucleus 是测量单位；供推断的独立单位必须由实际培养、供者、随机化与处理分配决定。代码中的字符串标识不能证明实验独立性。
- [预设检查计划](evidence/2026-09-11-code-science-audit/PLAN.md) 先于反例执行冻结；C8 饱和检查另行记录于执行之前。体积 `<1%`、面积 `<5%` 沿用 `SCIENTIFIC_SPEC.md` 第 9 节，未放宽。
- 分离级别：Tier D，**LIMITED INDEPENDENCE**。同一审核上下文执行源码重建与测试；解析几何、精确整数、单位变换是外部于被测实现的判定依据，但不使本次审核成为独立第三方验证。

## 2. 已复现发现

### F1 — P1：缺一个轴的标定会丢弃其他轴的已知信息，放过冲突

**位置：** `microscopy_io/tiff_contract.py:166–167, 204–209, 224–226`；相关 `microscopy_io/metadata.py:114–115`。

`ome_spacing()` 只要一个轴缺失便返回整个 `None`；`load_sample()` 仅在完整 `metadata_spacing` 存在时比较显式标定。并行的 `compare_registered_grid()` 也在任何轴缺失时直接返回 `partial`，没有先比较可比较的轴。

**反例 C2：** OME 明确记录 X=Y=1 µm、缺 Z；manifest 提供 Z=2、X=Y=2 µm。加载成功，返回 `(2,2,2)`，未报告已知 XY 冲突。另一对 `(1,1,None)` 与 `(2,2,2)` 的网格比较仅返回 `partial`。若 Z=2 正确，实际每体素体积应为 2 µm³，却会采用 8 µm³，差 4 倍；尺寸、分割物理阈值和体积均受影响。

**最小修复方向：** 保留逐轴可选值；先拒绝任一已知轴上的冲突，再判定整体是否完整；把同一约束应用于主图、独立通道、标签及真值读取。回归应覆盖每个轴缺失、其余轴匹配/冲突及真正全未知情况。只补缺失值，不隐式覆盖已知值。

**状态：FAIL，未修复。** [新反例脚本](evidence/2026-09-11-code-science-audit/probe.py)，结果 `C2_partial_metadata`。

### F2 — P1：表面积估计不满足既定精度，球形度随之偏低

**位置：** `quantification/features.py:101–107`，`surface_mesh()` 使用二值、无平滑、0.5 等值面的 marching cubes；各测量流程复用该实现。

**反例 C1：** 使用解析球体而非旧报告数值作判定，半径 10/20/40 µm，各取 `(1,1,1)` 与 `(2,1,1)` µm 的 ZYX 体素间距。六个面积相对误差均不满足 `<5%`。

| 半径 µm | Z 间距 µm，XY=1 | 体积相对误差 | 面积相对误差 | 球形度 |
|---|---|---|---|---|
| 10 | 1 | −0.472% | +9.184% | 0.9130 |
| 10 | 2 | −2.263% | +11.355% | 0.8844 |
| 20 | 1 | −0.326% | +8.483% | 0.9198 |
| 20 | 2 | −0.657% | +11.746% | 0.8910 |
| 40 | 1 | −0.120% | +8.676% | 0.9194 |
| 40 | 2 | −0.227% | +11.958% | 0.8918 |

球形度公式本身正确；偏差来自离散支持与表面估计。一个粗采样球体还不满足体积 `<1%` 要求。这不是说体素计数乘体素体积的算术有误，而是连续对象的离散近似未达到声明的精度。将所有间距乘 3 时，体积乘 27、面积约乘 9，说明单位缩放正确仍不足以证明准确性。

**最小修复方向：** 按当前标准保留失败状态；任何表面估计器变更须有明确测量定义、开发/保留验证分离、分辨率/形状/各向异性分层和方法版本。不能用放宽到现有回归测试的 15% 或截断球形度来消除失败。支持范围的改变也需要独立证据，不能据此次结果事后挑选范围。

**状态：FAIL，当前已知问题经新反例确认。** [probe.json 的 C1_geometry](evidence/2026-09-11-code-science-audit/probe.json)。[scikit-image 官方接口](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes)说明方法/spacing 语义，不提供本项目的精度保证。


**状态更新 2026-09-12：已关闭（限定适用域）。**
生产估计量已替换为 `crofton_minimax_sym_v3`（权重来自平面朝向上的 minimax
线性规划，取值前即给出先验误差预算）。在冻结后、确认集未被触碰的前提下，
96 个域内确认案例最差面积误差 0.951 %；项目自身冻结网格的 40 个域内案例最差
0.790 %（原估计量为 17.85 %）。适用域为声明并冻结的 ρ_in ≥ 10、anisotropy ≤ 4、
光滑闭合曲面；带二面角折痕的曲面在任何分辨率下均 `NOT QUALIFIED`。域内判定逐
对象导出并对域外对象加 review 标记，不拒绝也不静默通过；光滑性无法从 mask 机器
检验，由调用方负责。

未一并解决、必须照实记录的部分：SG-1 为无限定项，冻结网格含折痕与欠分辨形状，
因此仍为 FAIL（无限定最差 25.81 %，出现在 3 体素圆柱上，比原估计量的 18.74 %
更差——该估计量在其声明域内远更准确，在域外退化更快）。验收判据、phantom 与
gate 规则一概未改，域内主张由独立测试承载。体积仍为体素计数、未改动，SG-2 与
上一份记录逐位一致，并且在域内体积而非面积成为约束项（F2 中"球形度偏低"的
部分随之改变方向：准确面积不再抵消体素体积的正偏差，小对象球形度略高于 1，
按既有策略标记而不裁剪）。证据：`docs/evidence/2026-09-12-surface-crofton-v3`
（canonical），说明见 `docs/VALIDATION_UPDATE_2026-09-12.md`。

### F3 — P2：聚类推荐、拟合、评分用了不同的特征尺度

**位置：** `web_interface/analysis_ui.py:408–416`；`statistics/exploration.py:538–575`。

推荐 k 的调用传入原始特征，实际 `perform_kmeans_clustering()` 做 StandardScaler，而 `silhouette_for()` 又用原始特征计算分数。它们评价的几何空间不同。

**反例 C3：** 固定 180 行合成数据、随机种子与 k=3。仅将第二个特征单位乘 1000，标准化后的拟合标签完全相同，正确的拟合空间 silhouette 均约 0.409983；UI 对应评分从 0.389206 变为 0.178940，推荐 k 从 3 变成 8。标准化空间搜索在两种单位下均推荐 8；此处不声称 8 有生物学真实性。

**最小修复方向：** 推荐、最终拟合和评分共享同一填补、标准化及距离定义；为单位变换、缺失值处理和固定标签的评分一致性补充回归。silhouette 是输入距离空间的分离度，见 [官方定义](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html)。

**状态：FAIL，未修复。** 结果 `C3_clustering`。

### F4 — P2：含未分配核时，大整数父 ID 经浮点中间值损坏

**位置：** `workflows/multilevel_measurement_workflow.py:203`。

`Series.map(...).fillna(0).astype(int)` 的缺失匹配让映射结果经历 float64。当前输入契约允许 int64 ID；超过 2^53 的整数不一定可由 float64 精确表示。

**反例 C4：** organoid ID=`9007199254740993`；一个 cell 及一个核正确位于内部，另有一个无父核。cell 父 ID 保持正确，已分配核的父 ID 却变成 `9007199254740992`，真实父对象的 `nucleus_count` 从应有的 1 变成 0。附加 mismatch QC 能提示异常，但未防止错误聚合。

**最小修复方向：** 全程使用整数或 pandas nullable integer，不以浮点存储身份。回归覆盖无父核与大 ID 共存、相邻大 ID、精确计数及导出往返。该反例针对直接 multilevel API 的 int64 契约；不是声称正常 uint32 标签也必然出现这个损坏。

**状态：FAIL，未修复。** 结果 `C4_parent_identity`。

### F5 — P2：多层测量的几何异常未进入 QC，仍输出 pass

**位置：** `quantification/multilevel_relationships/qc.py:44–72`；测量来自 `multilevel_relationships/morphology.py`。

**反例 C5：** 图像内部的 2×2×2 标签满足默认 `minimum_voxels=5`，连续几何不可能的球形度为 1.192792，但 organoid 的 `qc_status="pass"` 且 flags 为空。同一标签在 Web mask features 路径被标记 `sphericity_above_geometric_range`。

这是科学/数值鲁棒性缺口，并非断言该模块违反了自己已实现的 flags 清单。`PARAMETERS.md` 已明确该阈值不是 multilevel 规则；文档披露不能使这个 `pass` 具备几何质量保证。

**最小修复方向：** 明确 `pass` 的有限含义，并对物理不可能/不可估计值给出机器可读状态；可评估共用现有 1.05 review 阈值，但不可把它解释为生物学验证边界。保留原始值，避免裁剪或默认删除。

**状态：FAIL（几何异常识别），未修复。** 结果 `C5_qc`。

### F6 — P2：默认整数饱和判定把存储位宽当作探测器上限

**位置：** `quantification/features.py:177–197`，默认 `quality.*_saturation_value=null`。

**反例 C8：** 已知 12-bit 上限 4095 的 Calcein ROI，100% 像素饱和，以 uint16 保存。默认逻辑以 65535 为上限，报告饱和比例 0、`viability_measurement_eligible=true`。显式提供已知上限 4095 后，比例为 1、eligibility=false。

**范围与现有缓解：** `configs/structural_fluorescence.yaml` 已提醒 12-bit 用户设置 4095，已有测试也验证显式设置有效。因此这是默认配置/采集契约风险；不是显式正确配置路径的新回归，也不证明真实样本的信号已饱和。不能从 dtype、图像最大值或名称推导真实探测器上限。

**最小修复方向：** 饱和上限必须来自可追溯采集元数据或用户明确声明；未知时保持 unknown/ineligible，而不是给出“已排除饱和”的结果。回归覆盖低有效位宽整型、显式上限与浮点强度。不能改变活性 gates 来掩盖饱和。

**状态：FAIL（该默认条件），显式上限路径 PASS。** [saturation_probe.py](evidence/2026-09-11-code-science-audit/saturation_probe.py) 与 [结果](evidence/2026-09-11-code-science-audit/saturation_probe.json)。

## 3. 算法科学依据与适用性

| 环节 | 当前实现与科学假设 | 依据强度 / 适用性 / 本次证据 |
|---|---|---|
| 图像读入 | ZYX/CZYX，显式时间点、物理单位、均匀 Z；形状/spacing 一致不足以证明配准 | 单位定义 ESTABLISHED；部分元数据完整性 FAIL（F1）；采集与配准资格 INSUFFICIENT EVIDENCE |
| 经典分割 | 强度阈值/Otsu、物理尺度平滑、闭运算、填洞、距离变换 watershed；要求结构信号对应目标边界 | 方法 CONTEXT_DEPENDENT；尚无目标图像精度支持，适用性 INSUFFICIENT EVIDENCE。核染色边界不能自动替代整类器官边界 |
| Cellpose | 4.2.1.1；默认 cpdino-vitb，do_3D，正交平面网络、diameter/anisotropy 与下采样尺度调整 | 模型数据适配 CONTEXT_DEPENDENT；接口回归 PARTIAL，真实精度 INSUFFICIENT EVIDENCE。未新运行网络或下载权重 |
| 体积、矩与表面 | 体素计数、含体素内二阶矩的等效椭球轴；表面自 2026-09-12 起为 minimax 权重的格点截线计数（marching cubes 保留为 legacy）；经典流程填洞，multilevel/Web 默认原始标签 | 公式 ESTABLISHED；支持定义不同必须保留；面积精度在声明域内 PASS、域外无主张、无限定 SG-1 仍 FAIL；体积精度 ρ_in < 12 仍 FAIL（F2 更新见下） |
| 层级/接触/位置 | 最大重叠分配，完整子对象体积作分母；6 邻接体素面面积；EDT 取舍入后质心处 | 离散定义 ESTABLISHED/CONTEXT_DEPENDENT；身份完整性 FAIL（F4）。体素接触是离散代理；EDT 是到背景体素中心，不等同连续膜面距离 |
| 活性代理 | 局部背景中位数扣除、批内独立对照汇总、双标记规则；保留负扣背景值与超对照范围标记 | HEURISTIC；只能称信号状态。饱和条件 FAIL（F6），独立 assay 适用性 INSUFFICIENT EVIDENCE |
| 推断与探索 | 生物重复随机截距 LMM，失败/退化时 clustered OLS；对象级分类、KMeans | CONTEXT_DEPENDENT/HEURISTIC；研究设计与统计校准 INSUFFICIENT EVIDENCE；聚类流程一致性 FAIL（F3） |

已读取安装版 Cellpose 的 `models.py`/`dynamics.py`：`flow_threshold` 在 do_3D 时不使用，当前 Web 已禁用并说明，未将其误列为缺陷；与 [官方 API](https://cellpose.readthedocs.io/en/latest/api.html)一致。默认 `min_size=15` 工作体素、`max_size_fraction=0.4`、normalize、迭代等仍影响结果。标量 diameter/anisotropy 修正不能证明不同下采样倍率的分割等价。

## 4. 参数/阈值依据与敏感性

从 `config.py`、`SegmentationConfig`、multilevel config 和当前安装库重建的重要选择包括：

| 参数组 | 当前值举例 | 来源与依据 | 校准/敏感性结论 |
|---|---|---|---|
| 分割 | Gaussian 1 µm、closing 1 µm、最小体积 2000 µm³、h=2 µm、seed distance=15 µm、probability=0.5 | 项目 HEURISTIC；具体阈值不是由方法文献自动确定 | 未找到合格独立标注支持这些精确值，INSUFFICIENT EVIDENCE |
| 几何采样 | 原始 spacing、level=0.5、padding=1、step=1、无平滑 | 实现选择 CONTEXT_DEPENDENT | C1 已检验分辨率/各向异性变化，六例面积 FAIL；不能据此宣称所有其他形状误差大小 |
| Cellpose | diameter=30/50 px、anisotropy=2.9、XY=0.414 µm、flow smoothing=1；库内 min_size/max_size_fraction | SOFTWARE_OR_MODEL_DEFAULT 及项目 HEURISTIC | 默认物理尺寸必须按采集资格验证；实际目标域基准、保留验证与稳定性 INSUFFICIENT EVIDENCE |
| 活性与 QC | gates=0.30/0.60、对照重复≥2、separation SNR≥3、背景壳2–7 µm、饱和比例1% | gates/阈值 HEURISTIC；运行时控制端点 CONTROL_DERIVED | 控制端点是校准，不是独立验证。未做跨批保留 assay 验证，精度与判别不确定性 NOT ASSESSED |
| 多层关系 | overlap=0.5、minimum_voxels=5、MAD z=3.5、core/peripheral=0.5/0.8等效半径 | HEURISTIC | 可追溯的操作定义；不等于组织区域或异常细胞的已验证边界 |
| 统计 | 每条件≥3重复、volume log10、LMM退化比1e-6、t(G−1)、BH按feature、bootstrap=2000 | 公式/成熟方法与项目 HEURISTIC 组合 | 未建立最小 N 的精度依据、I 类错误或 CI coverage；仅有可执行性/回归证据 |

未对科学阈值作事后调参。图像/标签中的有效对象数量也不能充当生物重复 N。未来验证须先定义 assay、对象/样本独立性、参考标准资格、所需精度和分层，再规定校准与保留验证方案。

## 5. 统计有效性、QC 传播与解释边界

`aggregation.py` 在完整 acquisition unit 上按 well 汇总，再按 condition/biological_replicate 汇总；重复跨 batch 不自动增加独立 N。`organoid_measurement_workflow.py` 把不完整 unit 排除出下游推断；这些设计比直接按对象行数作独立 N 更审慎，但不能补足缺失的真实实验设计。

`inference.py:61–71, 135–154` 的模型为 `feature ~ condition`，随机/聚类键是 condition 与 biological_replicate 的组合。未建立 batch 效应、额外 well 嵌套效应或跨条件配对结构；文档已将配对/跨条件重复测量列为不支持，故本次不把缺少这些模型项本身判作支持域内的新缺陷。它仍不能支持这些设计上的确认性处理效应。

LMM 的 t(G−1) 近似不是 Satterthwaite/Kenward–Roger；仅收敛不足以验证误差率。BH 只覆盖每个 feature 内的 contrasts，未覆盖所有 feature×contrast/omnibus。描述性等权 well 中位数与对象级模型也有不同 estimand。需要实际数据设计、模型诊断、效应/不确定性及相应模拟/独立验证。

`exploration.py` 虽已把填补/标准化放入训练折，仍按对象随机拆分，未隔离同一供者/井/图像。因而 Test Accuracy/AUC 不能作为生物样本外推性能；[scikit-learn 分组验证说明](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation-iterators-for-grouped-data)支持这一限制。UI 已有相应提示，故将生物泛化评估列为 INSUFFICIENT EVIDENCE。在构成聚类的同一特征上做 ANOVA 是数据选择后的描述，不能当作独立发现群体的显著性证据。

经典分割方法不一致警告保留于 sample/segmentation_qc，但不会自动让所有对象 `morphology_eligible=false`。multilevel 的低重叠、无核、多核等 flags 也不自动排除聚合。当前文档/UI披露了保留 flagged 对象的策略，因此 C6 为 PARTIAL，不凭主观偏好把“未自动删除”判作算法错误；使用者仍须预设审核与纳入规则，不能把 complete/pass 解释为科学有效。

## 6. 真实数据与独立标准

实际遍历包含 Git 忽略文件的 `data/`，发现五个真实 TIFF 和一个二维 JPEG，详见 [实时清单](evidence/2026-09-11-code-science-audit/real_data_inventory.json)。五个 TIFF 都是 ImageJ ZYX；HCT116-C1、PDAC-C1、PDAC-C2 缺可解析 Z spacing，另两个文件含 Z spacing 1.5/2.0 µm。文件名只作身份记录，不用作已确认的染色/组织来源。

合成 demo 的图像、truth、标签与输出仍存在；它们可作受控工程证据。当前目录中未找到与这些真实图像匹配、具备标注协议/盲法/一致性/裁决记录的独立三维 organoid/cell/nucleus 标注集，也未找到正交活性 assay 或足以确认供者、培养、well、batch 独立性的实验设计。文件 hash 不证明这些资格。

因此：真实分割准确性、罕见/困难场景检测、跨批活性代理有效性和实际研究推断均为 INSUFFICIENT EVIDENCE 或 NOT ASSESSED。此次只重新核实文件与元数据，没有重新跑真实图像网络；即使跑通也不能补足上述标准。未替这些缺失研究臆造最小 N 或允许误差。

## 7. 工程、数值与科学鲁棒性分别评价

| 项目 | 本次观测 | 状态 |
|---|---|---|
| 工程回归 | `pixi run --locked python -m pytest -q --junitxml=...`：436 passed，8 skipped，0 failure/error，54.527 s；8个浏览器布局/渲染测试因缺 Playwright 跳过；render marker 的原生渲染未运行 | PASS（已执行范围） |
| Lint | `pixi run --locked lint` exit 0 | PASS |
| 类型检查 | `typecheck-ratchet` exit 0：baseline/current 110，新增0。内部类型34、缺第三方stubs44、科学库边界32 | PARTIAL；不能说 mypy 无错误 |
| Notebook | `check` exit 0，明确没有 notebook | NOT APPLICABLE |
| 独立反例 | 两个脚本均 exit 0，成功生成严格 JSON；脚本成功是复现成功，不代表被测性质通过 | F1–F6 见上文 |
| 科学/数值鲁棒性 | 物理缩放符合关系；部分元数据、面积精度、特征空间一致性、大ID和QC等失败 | FAIL |
| 工程鲁棒性总体 | 现有失败路径/导出/并发相关回归已随测试执行；未新做负载压力、实际 GUI/native、干净安装发行包等完整评估 | PARTIAL |
| 科学门禁 | `science-gate` exit 1，0/6 PASS；它核查已有封存记录，没有重新计算188个案例 | FAIL（门禁结果）；新科学证明来自C1等反例 |

所有命令状态在 [command_exit_codes.json](evidence/2026-09-11-code-science-audit/command_exit_codes.json)，原始日志、pytest XML 与反例 JSON 在同目录。工程绿灯与科学有效性独立报告。

## 8. 追溯矩阵与后续闭环

| 要求 | 源码/受影响输出 | 本次独立判定依据 | 观测/状态 | 必要闭环 |
|---|---|---|---|---|
| C1 几何精度 | features.geometry；volume、area、sphericity | 6个解析球体、3倍尺度关系 | area 6/6失败，volume 1/6失败；FAIL | 保持标准，方法验证或有证据的适用域资格 |
| C2 已知单位冲突 | tiff_contract / metadata；所有物理测量 | 手工构造的部分OME和不一致manifest | 接受冲突；FAIL | 逐轴校验及跨入口回归 |
| C3 聚类几何一致 | analysis_ui / exploration；k推荐、silhouette | 单位变换与拟合空间独立评分 | 标签不变而评分/推荐变；FAIL | 共享预处理与距离空间 |
| C4 父子身份 | multilevel workflow；核计数/密度 | 精确int64 ID和手工计数 | 父ID损坏、漏计；FAIL | 无浮点身份映射及端到端回归 |
| C5 几何异常状态 | multilevel QC | 8体素标签、现有物理范围 | 球形度1.1928且pass；FAIL | 机器可读异常、有限pass语义 |
| C6 QC到推断 | workflow、aggregation、report | 源码追踪/当前回归 | flags保留但不总排除；PARTIAL | 按实际研究预设纳入规则 |
| C7 工程检查 | 当前源码/测试/锁文件 | 实际命令与退出值 | 测试/lint通过；类型债110；PARTIAL总体 | 改动后受影响回归，发布另验发行包 |
| C8 采集饱和 | marker_measurements；eligibility | 已知12-bit上限与显式配置对照 | 默认误报未饱和；FAIL | 上限资格/未知状态契约 |
| 生物分割/活性/推断 | 全流程 | 合格独立标注、正交assay、独立实验单位 | INSUFFICIENT EVIDENCE | 冻结科学验证计划后获取并评估独立数据 |

本次 Gate 1/2 完成当前仓库与影响重建；Gate 3/4 核对现有规范与方法参数（没有实施新设计）；Gate 5 冻结反例标准；Gate 6/7/8 审阅实现并执行相应检查；Gate 9 明示审核分离限制；Gate 10 完成分级与最小修复建议，按 audit-only 范围不修改生产逻辑；Gate 11 记录源码/锁文件/数据身份；Gate 12 以本报告及矩阵给出受范围约束结论。无新实现的设计审批、无发布的发行物/部署评估均不在本次范围。

F2 已于 2026-09-12 在限定适用域内关闭（无限定 SG-1 仍 FAIL，体积未改动）；优先关闭 F1；F3/F4/F5 的修复可保持既有科学模型与阈值，F6 需要明确采集上限契约。修复工程问题后仍须分别建立真实图像分割、活性代理和实验设计的科学证据。**当前总体结论仍是 NOT READY FOR THE SPECIFIED RESEARCH USE。**
