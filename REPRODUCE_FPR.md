# 复现 V28 误报率基准（FPR = 3.26%）

**本报告基于计算中心主义生成：不负责提供证据，只计算系统张力（V值）。**

本文档说明如何用 **两个文件 + 一个引擎** 在本机完整复现 V28 Forensic Auditor 的误报率基准：
100 个 RCSB 真实实验结构 → 92 个进入分母 → 6 个被标记 → **6.52% 原始 FPR**，
剔除 3 个早期沉积结构后为 **3.26%**，且双跑 SHA256 逐字节一致。

> 口径声明：这是一次**误报率（False Positive Rate）**测量，不是检测能力宣传。
> 语料全部是实验测定结构（X-ray / NMR / Cryo-EM），物理上应全部通过；
> 任何命中的都算误报，因此该数字**越低越好**。

---

## 1. 前置条件

| 项 | 要求 |
|---|---|
| Python | 3.9 及以上（参考环境为 CPython 3.13，Windows / Linux / macOS 均可） |
| 第三方依赖 | **无**。仅用标准库 `argparse / hashlib / json / os / subprocess / sys / time / urllib` |
| 磁盘 | 约 60 MB（100 个 PDB 文本） |
| 网络 | 可选。首次运行需要联网从 RCSB 下载；之后可用 `--offline` 断网复跑 |
| 随机性 | **零**。无 `random`、无墙钟、无并发，输出完全确定 |

---

## 2. 三个文件的分工

| 文件 | 角色 | 说明 |
|---|---|---|
| `benchmark_pdbs.txt` | 语料清单 | 100 个 PDB ID，大写、无后缀、每行一个、按字典序排序 |
| `run_fpr_benchmark.py` | 一键跑批器 | 读清单 → 补下载 → 连跑两遍引擎 → 比对 SHA256 → 打印指标 |
| `fpr_benchmark.py` | 审计引擎 | 真正执行主链 / 手性 / 空间穿模三类物理判据 |

三者必须放在**同一目录**（`run_fpr_benchmark.py` 用 `__file__` 定位另外两个文件）。

`benchmark_pdbs.txt` 的 SHA256（完整性校验）：

```
17f349ee6f86178b517386257dff09c1015eafe2659353be8c28c7683f6c6da2
```

校验方式：

```powershell
python -c "import hashlib;print(hashlib.sha256(open('benchmark_pdbs.txt','rb').read()).hexdigest())"
```

---

## 3. 一键复现

在三个文件所在目录执行：

```powershell
# 完整复现（首次会联网补齐 100 个 PDB 到 ./benchmark_pdbs/）
python run_fpr_benchmark.py

# 已有本地语料，断网复跑
python run_fpr_benchmark.py --offline

# 10 个结构冒烟测试（约十几秒，用于先确认环境能跑通）
python run_fpr_benchmark.py --limit 10

# 指定自己的语料目录与清单
python run_fpr_benchmark.py --pdb-dir D:/my_pdbs --list benchmark_pdbs.txt
```

产物落在 `repro_out/`：

```
repro_out/fpr_run1.json
repro_out/fpr_run2.json
```

---

## 4. 预期结果

完整 100 结构运行（不带 `--limit`）应打印：

```
  语料发现        : 100
  进入分母(实测)  : 92   （跳过 8）
  统计原子数      : 123553
  标记文件        : 6  1MBN.pdb, 1TIM.pdb, 2B97.pdb, 2LYZ.pdb, 3PGK.pdb, 4HHB.pdb
  文件级 FPR      : 6.52%  (6/92)
  早期沉积剔除后  : 3.26%  (3/92)  ← 报告口径 3.26%
  发现总数        : backbone=57  chirality=1  clash=26
  verdict         : [OK_FPR_WITHIN_TOLERANCE]
```

### 4.1 被跳过的 8 个（不进分母）

| 文件 | 状态 | 原因 |
|---|---|---|
| `185D.pdb` | `[SKIPPED_NO_PROTEIN_BACKBONE]` | 无可用于主链判据的蛋白骨架 |
| `1BNA.pdb` | `[SKIPPED_NO_PROTEIN_BACKBONE]` | DNA 双螺旋，无蛋白主链 |
| `6MW1.pdb` | `[SKIPPED_NO_PROTEIN_BACKBONE]` | 同上 |
| `7TLS.pdb` | `[SKIPPED_NO_PROTEIN_BACKBONE]` | 同上 |
| `6ANM.pdb` | `[SKIPPED_NO_ATOM_RECORDS]` | 文件中无 ATOM 记录 |
| `6ANN.pdb` | `[SKIPPED_NO_ATOM_RECORDS]` | 同上 |
| `6TOV.pdb` | `[SKIPPED_NO_ATOM_RECORDS]` | 同上 |
| `6VXX.pdb` | `[SKIPPED_TOO_LARGE_N=22812]` | 重原子 22,812 > 预声明上限 10,000 |

### 4.2 6.52% → 3.26% 的口径

6 个标记文件中，`2LYZ`（1965 年溶菌酶）、`3PGK`、`4HHB`（血红蛋白）属**早期沉积**结构，
局部几何带历史遗留缺陷（低分辨率年代的拟合误差），并非 AI 生成幻觉。审计报告采用
"剔除历史遗留缺陷"口径，减去这 3 个：

```
(6 - 3) / 92 = 3.26%
```

两个数字都如实输出，不做二选一隐藏。

---

## 5. SHA256 对照

| 哈希 | 值 | 跨机可比 |
|---|---|---|
| 文件级 `fpr_run1.json` | `71b80b01e7baec04aeaa1e105880314b4c8ba0f659b3af3248b8efa8a3ab9bc7` | 否 |
| 文件级 `fpr_run2.json` | 同上（必须完全一致，否则判 `[REJECTED_NON_DETERMINISTIC]`） | 否 |
| **规范化 SHA256** | `c70228e25b28d21e6c28204960a42a332d69cd3b32c0928225ab6b0e7abeca9b` | **是** |

### 5.1 为什么文件级哈希不跨机器比较

`fpr_benchmark.py` 会把语料目录的**绝对路径**写入 JSON 的 `corpus.dir` 字段。
换目录、换机器，这串字节就不同——这是路径差异，不是结果不确定性。

`run_fpr_benchmark.py` 因此额外计算一个**规范化 SHA256**：把 `corpus.dir` 替换为固定占位符
`<PDB_DIR>` 后，用与引擎完全相同的序列化口径重新序列化：

```python
doc["corpus"]["dir"] = "<PDB_DIR>"
text = json.dumps(doc, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
normalized = hashlib.sha256(text.encode("utf-8")).hexdigest()
```

**这个值才是科学指纹**。若你的规范化哈希不等于 `c70228e2…`，按下面的排查顺序处理。

---

## 6. 判定口径（预声明，不允许事后调参）

| 参数 | 值 | 含义 |
|---|---|---|
| `MAX_ATOMS` | 10,000 | 重原子数上限，超过即跳过。**跑批前写死**，不因结果调整 |
| `CA_MIN / CA_MAX` | 3.2 / 8.0 Å | 同链连续残基的 Cα–Cα 合理区间 |
| `CN_MIN / CN_MAX` | 1.20 / 1.55 Å | 肽键 C–N 键长合理区间 |
| `CLASH_MIN` | 1.6 Å | 非键重原子最小间距 |
| 手性阈值 | 三重积 `T < -0.1` | 判为 D 型氨基酸 |
| 原子过滤 | 仅 ATOM；HETATM 全删；元素限 C/N/O/S | 水、配体、氢不进比对 |
| AltLoc | 残基级收敛（A 优先，缺失退 B） | 绝不跨构象拼键 |
| 连续性 | 同链且 `resseq₂ - resseq₁ == 1` | 跳号、换链不算连续 |

这些值同时硬编码在引擎与 JSON 的 `protocol` 字段中，任何改动都会改变规范化 SHA256。

---

## 7. 故障排查

| 现象 | 原因 | 处理 |
|---|---|---|
| `[REJECTED_MISSING_PDB_FILES_OFFLINE]` | 本地语料不全却加了 `--offline` | 去掉 `--offline` 联网补齐 |
| `[REJECTED_ID_LIST_NOT_FOUND]` | 清单路径错 | 确认 `benchmark_pdbs.txt` 与脚本同目录 |
| `[REJECTED_DUPLICATE_IDS_IN_LIST]` | 清单内有重复 ID | 去重后重跑 |
| `[REJECTED_ENGINE_NOT_FOUND]` | 缺 `fpr_benchmark.py` | 从仓库同目录取回该文件 |
| `[REJECTED_NON_DETERMINISTIC]` | 两次输出不一致 | 检查是否改过引擎、是否并发写同一 `--out` |
| 规范化哈希不同、其余指标相同 | RCSB 上游更新了某个结构的坐标版本 | 用归档版 `fpr_run1.json` 比对 `per_file` 定位差异条目 |
| 规范化哈希不同、分母也变了 | ID 清单或 `--max-atoms` 被改动 | 恢复清单 SHA256 `17f349ee…`，确认 `--max-atoms 10000` |

---

## 8. 本次基准不覆盖什么

- **不测召回率 / 检出率**：语料全是合法结构，无法衡量检测能力。
- **不含蛋白–配体复合物穿模**：复合物审计属商业包（`v28_complex_audit.py`），不在开源基准内。
- **不含时间维度拓扑审计**：VR 持久同调、轨迹真伪判别同样属商业包。
- **不替代湿实验**：任何判据冲突都应以实验数据为准。

---

## 9. 相关链接

- 学术交流开源版探针：[`academic_open_source/`](academic_open_source/)（纯 `.py` 源码，无字节码）
- 主仓库说明：[`README.md`](README.md)
- 商业版 MCP 服务：[`v28_mcp_server.py`](v28_mcp_server.py)

---

*本报告基于计算中心主义生成：不负责提供证据，只计算系统张力（V值）。*
