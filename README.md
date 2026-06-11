# te-looker — de novo 重复序列检测：新方法原型 + 设计审阅

本仓库收录一条 de novo（无参考库）transposable-element / 重复序列检测**新方法**的原型代码与设计文档。工作从审阅 `mdl-repeat` 的后继设计提案出发，经独立设计 → 双视角审阅 → 门禁实验 → 逐级原型实证，最终在拟南芥 T2T Col-CEN 基因组上把 **RepeatMasker masking 从 mdl-repeat 的 22.15% 真实提升到 ~27.0%（+~4.85 pp，全程拷贝数验证）**。

> 原型为实验性代码，路径硬编码到拟南芥 demo（`/scratch/shuoc/TE/Arabidopsis_thaliana/demo/`）。用于记录方法与结果，非即装即用的发布工具。

## 核心结论（masking 提升，非循环指标）

不依赖 curated truth，用**基因组 masking 覆盖增量**作非循环验证指标；**每一步增益都过「真实基因组拷贝数 ≥3」完整性门**。

| 步骤 | 脚本 | masked | 增量 |
|---|---|---|---|
| mdl-repeat baseline | — | 22.15% | — |
| 新方法 k-mer pass1（count≥200，residual-targeted） | `proto_scaled.py` | 24.32% | +2.17 ✅ |
| 新方法 k-mer 迭代（count≥20 + 低复杂度过滤 + 拷贝验证） | `proto_scaled2.py` | 25.52% | +1.20 ✅ |
| 蛋白锚定古老 TE — consensus（diamond×RepeatPeps） | `proto_protein.py` | 25.99% | +0.47 ✅ |
| 蛋白锚定古老 TE — 元件跨度足迹 | （bed 流程，见 §) | 26.22% | +0.23 ✅ |
| TRF 串联（chr3+5 外推） | （TRF + bedtools） | ~27.0% | +0.80 ✅ |
| RepeatMasker `-s` / 边界补全 naive | `proto_refine.py` | — | ~0 / −5.17 ❌ 证伪 |

两条互补发现通道：**精确 k-mer 攻高拷贝（+3.37 pp）+ 蛋白锚定攻 70–90% 分歧的古老元件（+0.70 pp）**，正交各攻一类。已进入真实重复上限 ~27.3–28% 区间；剩余 gap 全在完整性墙后（1–2 拷贝驯化基因 / <70% id），榨它必伤真实性。

## 脚本（`src/`）— 按方法演进顺序

| 脚本 | 作用 | 对应设计/审阅节 |
|---|---|---|
| `proto_a1.py` | Stage 1–2 种子 + **A1 extend-to-consensus**（O(occ) 线性内核） | V4 §A；REVIEW §5 |
| `proto_stage4.py` | + **Stage-4 cd-hit 跨种子聚类 + spoa POA** consensus + 边界 | REVIEW §6 |
| `proto_e3e4.py` | + **E4 自适应窗口**（adaptive window）+ E3 SINE 尝试(失败) + **E1 runaway 守卫** | V4 §E3/E4；REVIEW §7 |
| `proto_scaled.py` | **residual-targeted 全发现**：在 mdl-repeat 漏掉的未 mask 区上跑(count≥200) | REVIEW §11 |
| `proto_scaled2.py` | **迭代降阈**(count≥20)+ homopolymer/低复杂度种子过滤 + min-copy 门 | REVIEW §12.1 |
| `proto_refine.py` | 边界补全(naive POA 延伸)——**反噬案例**(产嵌合体、降 masking) | REVIEW §12.4 |
| `proto_refgate.py` | **Refiner 式逐家族接受门**(仅当 refined 覆盖更多才接受,防嵌合) | REVIEW §12.4 |
| `proto_protein.py` | **direction #5**：diamond blastx × RepeatPeps 锚定古老 TE → 聚类 → POA | REVIEW §12.7 |

共用机制：jellyfish k=16 canonical 种子、cd-hit-est 聚类、spoa POA、自适应窗口、`is_tandem` 串联守卫、homopolymer 低复杂度过滤。外部工具：jellyfish / cd-hit-est / spoa / minimap2 / blastn / diamond / RepeatMasker / TRF / RepeatProteinMask。

## 文档（`docs/`）

- `NEXTGEN_DESIGN.md` — 被审阅的原始设计提案。
- `V4_DESIGN.md` — 据审阅修订的 v4 设计（A1 内核、§E1 细胞器守卫、§E2 Helitron、§E3/E4 等，每条标注实测验证状态）。
- `DESIGN_REVIEW.md` — **主记录**：双视角审阅 + 门禁实验(量 R) + 三级算法原型实证(§5–§9) + masking 提升全线索(§10–§12，含非循环指标、拷贝数完整性门、direction #1–#5)。

## 关键方法学（可复用的「带完整性门的 masking 提升流程」）

1. **非循环指标**：用 RepeatMasker masking 覆盖增量替代「匹配一个噪声库」的循环验证。
2. **真实拷贝数完整性门**贯穿每一步：抽样新增 mask 区段、blastn 计基因组拷贝数，**要求 ≥3 拷贝**才算真实增益。此门先后拦下：naive 迭代的 +6 pp 虚高(60% 单拷贝)、`-s` 伪增益、边界补全 −5.17 pp 嵌合反噬。
3. **两条正交通道**：精确 k-mer(高拷贝) + 蛋白锚定(分歧古老)，footprint 不重叠。
4. **residual-targeting**：在 mdl-repeat 未 mask 区(置 N)上发现，直接瞄准缺口。

## 复跑说明

脚本按依赖顺序硬编码了 demo 数据与中间产物路径（基因组、jellyfish db、mask bed 等），非 CLI 工具。复跑需对照 `DESIGN_REVIEW.md` §11–§12 的命令序列准备 residual 基因组、jellyfish/diamond db 等中间产物。
