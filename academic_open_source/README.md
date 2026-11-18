# academic_open_source —— 学术交流开源版

本目录为 **学术交流开源版**，如需 MCP 及复合物穿模商业包，请联系作者。

---

## 本目录包含什么

| 文件 | 作用 |
|---|---|
| `probe_backbone_hallucination.py` | 主链拓扑连续性校验（Cα–Cα 间距、肽键 C–N 键长、连续越界熔断、残基跳号过滤） |
| `probe_chirality_violation.py` | 立体化学合规性校验（手性标量三重积，识别 D 型氨基酸） |
| `probe_vdw_overlap.py` | 空间位阻防碰撞（非键合重原子范德华重叠检测，二硫键豁免） |
| `op_full_atom_grid.py` | 三维空间哈希网格（把 O(N²) 全原子比对降到近 O(N)），上面三个探针的共享底座 |

全部为**可读的 `.py` 纯源码**，无 `.pyc`、无混淆、无外部依赖（仅 Python 标准库）。

## 它不是什么

本目录**不含**以下商业能力，这些不在开源范围内：

- 蛋白–配体复合物空间穿模审计（含共价键模式白名单、链邻接扩展、极性原子对豁免）
- 时间维度拓扑寿命审计（VR 持久同调、空腔体积突降检测）
- 轨迹与能量日志的法医级真伪判别（隐藏时间线、复制帧、能垒曲率）
- MCP Server 封装、输出脱敏层、字节码发布打包

## 最小可运行示例

```python
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import probe_backbone_hallucination as bb

report = bb.scan_backbone_hallucination("1UBQ.pdb")
print(report["status"])        # [OK_BACKBONE_PHYSICAL]
print(report["n_out_of_bounds"])
```

可在本机跑通的自检 / 命令行入口（均在**本目录内**执行）：

```powershell
# 主链探针：8 项内置自检（理想链 / 坍塌 / 孤立断链 / 肽键断裂 / 空文件 / 双跑确定性 / 白盒纯度）
python probe_backbone_hallucination.py --selftest

# 立体化学与空间位阻探针：传入任意真实 PDB 路径作为参数
python probe_chirality_violation.py .\your_structure.pdb
python probe_vdw_overlap.py .\your_structure.pdb
```

> 说明：`probe_backbone_hallucination.py` 的 CIF 分支依赖未随本目录分发的 `op_cif_parser`，
> 缺失时自动降级为仅支持 `.pdb/.ent`，不影响 PDB 路径下的任何判定结果。

## 引用与联系

本目录代码可用于教学、复现与二次研究；商用（含 SaaS、内部生产流水线）请先联系作者获取授权。
误报率基准与复现方法见仓库根目录 `REPRODUCE_FPR.md`。

---

*本目录为学术交流开源版，如需 MCP 及复合物穿模商业包，请联系作者。*
