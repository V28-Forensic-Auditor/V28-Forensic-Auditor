# -*- coding: utf-8 -*-
"""probe_backbone_hallucination.py
TASK_20261107_V28_GENERATIVE_FORENSICS_FINAL —— AI 骨架幻觉探针。

纯 stdlib（math, os, re），零 random/numpy/scipy/Bio。
import 时硬编码正则扫描，违反者立即拒载。

物理铁律（强信号锚点）：
  蛋白肽键是化学共价键，键长由量子力学严格约束：
    · 相邻 Cα 距离 ≈ 3.8 Å（范围 3.5 ~ 4.2 Å）
      —— 这是肽平面 + Cα 四面体几何的必然结果
    · C-N 肽键 ≈ 1.33 Å（范围 1.23 ~ 1.43 Å）
      —— 部分双键性质，键长极刚性
  AI 生成蛋白常出现的幻觉：
    · 把 Cα 硬压到 < 3.2 Å（几何坍塌）
    · 把 Cα 拉到 > 5.0 Å（几何断裂）
    · 肽键 C-N 被拉到单键长度 1.5 Å 以上（物理荒谬）

防误杀容错（ExperienceRecall 硬教训）：
  膜蛋白跨膜区、柔性 Loop 区、结晶对称性接触可能出现孤立的
  越界残基对。只有**连续 >= 3 个残基**的骨架严重崩坏才触发
  CRITICAL，孤立残基触发 WARNING。
  —— 把"连续性"构造进检测逻辑，不是先判定后豁免。

红线：只扫描 PDB ATOM 记录，不修改输入。
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# —— 白盒纯度检查 ——
for mod in ("random", "numpy", "scipy", "Bio", "pandas"):
    if re.search(r"(?m)^\s*(import|from)\s+%s" % mod,
                 open(__file__, encoding="utf-8").read()):
        raise SystemExit("probe_backbone_hallucination.py 纯度自检失败：import %s" % mod)

# ========== 物理锚点常量（硬编码，来自 Cryst11 共价键长统计）==========
CA_DIST_MIN = 3.5    # Å，相邻 Cα 距离物理下限
CA_DIST_MAX = 4.2    # Å，相邻 Cα 距离物理上限
CA_DIST_CRIT_MIN = 3.2   # 严重坍塌（< 3.2 Å → 必然是幻觉）
CA_DIST_CRIT_MAX = 4.8   # 严重拉断（> 4.8 Å → 必然是幻觉）

CN_DIST_MIN = 1.23   # Å，C-N 肽键物理下限
CN_DIST_MAX = 1.43   # Å，C-N 肽键物理上限
CN_DIST_CRIT_MAX = 1.55  # 严重拉断

# 连续性阈值：连续 >= N 个越界残基 → CRITICAL
CRITICAL_CONTINUITY = 3

STATUS_OK = "[OK_BACKBONE_PHYSICAL]"
STATUS_TENSION = "[WARNING_BACKBONE_TENSION]"
STATUS_HALLUCINATION = "[CRITICAL_AI_UNPHYSICAL_BACKBONE]"
STATUS_EMPTY = "[NOT_APPLICABLE_NO_BACKBONE]"


# ========== 共享 PDB 解析（仅 ATOM 记录，纯 stdlib）==========

def _parse_pdb_atoms(pdb_path):
    """解析 PDB，返回蛋白原子列表。每个原子 dict：
    {"atom_name", "resname", "chain", "resseq", "x", "y", "z", "bfactor"}。
    只保留标准 20 种氨基酸的主链原子（N, CA, C, O）。
    """
    if not pdb_path or not os.path.exists(pdb_path):
        return None
    protein_residues = frozenset({
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE", "CYX", "HIE", "HID", "HIP",
    })
    backbone_atoms = frozenset({"N", "CA", "C", "O"})
    atoms = []
    # 格式分派：CIF/PDBx 复用升级后的解析层，PDB 路径逻辑保持不变。
    # 学术开源版说明：op_cif_parser 依赖内部公共层（op_l77_common），未随本目录
    # 分发；若同一目录/仓库中存在该模块则自动启用 CIF 支持，否则仅支持 .pdb/.ent。
    try:
        import op_cif_parser as _cif
    except ImportError:
        _cif = None
    if _cif is not None and _cif.is_cif_path(pdb_path):
        for a in _cif.parse_atoms_cif(pdb_path, keep_hetatm=False):
            if a["record"] != "ATOM":
                continue
            if a["resname"] not in protein_residues:
                continue
            if a["resseq_int"] is None:
                continue
            atoms.append({
                "atom_name": a["name"].upper(),
                "resname": a["resname"],
                "chain": a["chain"],
                "resseq": a["resseq_int"],
                "x": a["xyz"][0], "y": a["xyz"][1], "z": a["xyz"][2],
                "bfactor": a["bfactor"],
            })
        return atoms
    try:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("ATOM  "):
                    continue
                try:
                    aname = line[12:16].strip()
                    rname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resseq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    bf = float(line[60:66].strip()) if len(line) >= 66 and line[60:66].strip() else 0.0
                except (ValueError, IndexError):
                    continue
                if rname.upper() not in protein_residues:
                    continue
                atoms.append({
                    "atom_name": aname.upper(),
                    "resname": rname.upper(),
                    "chain": chain,
                    "resseq": resseq,
                    "x": x, "y": y, "z": z,
                    "bfactor": bf,
                })
    except Exception:
        return None
    return atoms


def _group_by_residue(atoms):
    """把原子按残基（chain + resseq）分组，返回 list[dict]。
    每个残基 dict: {"chain", "resseq", "resname", "atoms": {atom_name: (x,y,z)}}。
    只保留有主链原子（N, CA, C）的残基。"""
    if not atoms:
        return []
    residues = {}
    order = []
    for a in atoms:
        key = (a["chain"], a["resseq"])
        if key not in residues:
            residues[key] = {
                "chain": a["chain"],
                "resseq": a["resseq"],
                "resname": a["resname"],
                "atoms": {},
            }
            order.append(key)
        residues[key]["atoms"][a["atom_name"]] = (a["x"], a["y"], a["z"])
    # 按 resseq 排序
    sorted_keys = sorted(order, key=lambda k: (k[0], k[1]))
    result = []
    for k in sorted_keys:
        r = residues[k]
        # 必须有 N, CA, C 三个主链原子
        if all(n in r["atoms"] for n in ("N", "CA", "C")):
            result.append(r)
    return result


def _dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)


# ========== 主扫描 ==========

def scan_backbone_hallucination(pdb_path):
    """骨架幻觉主探针。返回标准结果 dict。"""
    atoms = _parse_pdb_atoms(pdb_path)
    if not atoms:
        return {
            "status": STATUS_EMPTY,
            "reason": "NO_PROTEIN_ATOMS_OR_PDB_NOT_FOUND",
            "probe": "probe_backbone_hallucination",
            "probe_fingerprint": "backbone_hallucination_v2_ca_cn_bond_length_continuity",
        }

    residues = _group_by_residue(atoms)
    if len(residues) < 2:
        return {
            "status": STATUS_EMPTY,
            "reason": "TOO_FEW_RESIDUES_N=%d" % len(residues),
            "probe": "probe_backbone_hallucination",
            "probe_fingerprint": "backbone_hallucination_v2_ca_cn_bond_length_continuity",
        }

    # 遍历相邻残基对，检查 Cα-Cα 和 C-N
    out_of_bounds = []  # list of dict: {res_i, res_j, ca_dist, cn_dist, ca_severity, cn_severity}
    ca_severity_flags = []  # 每个相邻对的 CA 严重度 (0=ok, 1=mild, 2=crit)
    cn_severity_flags = []

    for idx in range(len(residues) - 1):
        ri = residues[idx]
        rj = residues[idx + 1]

        # Cα-Cα 距离
        ca_i = ri["atoms"]["CA"]
        ca_j = rj["atoms"]["CA"]
        ca_d = _dist(ca_i, ca_j)

        # C-N 肽键距离
        c_i = ri["atoms"]["C"]
        n_j = rj["atoms"]["N"]
        cn_d = _dist(c_i, n_j)

        ca_sev = 0
        if ca_d < CA_DIST_CRIT_MIN or ca_d > CA_DIST_CRIT_MAX:
            ca_sev = 2  # critical
        elif ca_d < CA_DIST_MIN or ca_d > CA_DIST_MAX:
            ca_sev = 1  # mild

        cn_sev = 0
        if cn_d > CN_DIST_CRIT_MAX:
            cn_sev = 2
        elif cn_d < CN_DIST_MIN or cn_d > CN_DIST_MAX:
            cn_sev = 1

        ca_severity_flags.append(ca_sev)
        cn_severity_flags.append(cn_sev)

        if ca_sev > 0 or cn_sev > 0:
            out_of_bounds.append({
                "res_i_chain": ri["chain"],
                "res_i_seq": ri["resseq"],
                "res_i_name": ri["resname"],
                "res_j_chain": rj["chain"],
                "res_j_seq": rj["resseq"],
                "res_j_name": rj["resname"],
                "ca_dist_A": round(ca_d, 4),
                "cn_dist_A": round(cn_d, 4),
                "ca_severity": ca_sev,
                "cn_severity": cn_sev,
            })

    # ====== 连续性检测（防误杀核心逻辑）======
    # 把 severity 序列合成一条：取 max(ca_sev, cn_sev) per pair
    pair_severity = [max(a, b) for a, b in zip(ca_severity_flags, cn_severity_flags)]

    # 找连续 >= CRITICAL_CONTINUITY 个 pair_severity >= 2 的区段
    # 连续 pair 意味着：如果 pair[i], pair[i+1], pair[i+2] 都是 critical → 残基 i, i+1, i+2, i+3 跨了 3 个 pair
    # 所以连续 N 个 critical pair 涉及 N+1 个残基
    critical_runs = []
    i = 0
    while i < len(pair_severity):
        if pair_severity[i] >= 2:
            j = i
            while j < len(pair_severity) and pair_severity[j] >= 2:
                j += 1
            run_len = j - i
            if run_len >= CRITICAL_CONTINUITY:
                critical_runs.append({
                    "start_res_seq": residues[i]["resseq"],
                    "end_res_seq": residues[j]["resseq"],
                    "run_len_pairs": run_len,
                    "n_residues_involved": j - i + 1,
                })
            i = j
        else:
            i += 1

    # 孤立残基（没有连续 >= 3，但有越界）→ 统计总数
    total_out_of_bounds_pairs = len(out_of_bounds)
    critical_pairs = sum(1 for s in pair_severity if s >= 2)

    # ====== 判定 ======
    if critical_runs:
        status = STATUS_HALLUCINATION
        reason = ("CONTINUOUS_BACKBONE_COLLAPSE_NRUNS=%d_RESIDUE_SPANS=%s_CA_RANGE=[%.1f,%.1f]_CN_RANGE=[%.2f,%.2f]"
                  % (len(critical_runs),
                     "&".join("RES%d-%d" % (r["start_res_seq"], r["end_res_seq"]) for r in critical_runs[:3]),
                     CA_DIST_MIN, CA_DIST_MAX, CN_DIST_MIN, CN_DIST_MAX))
    elif total_out_of_bounds_pairs > 0:
        status = STATUS_TENSION
        reason = ("ISOLATED_BACKBONE_TENSION_NPAIRS=%d_NCRITICAL=%d_RESIDUES_SCATTERED_CA_RANGE=[%.1f,%.1f]_CN_RANGE=[%.2f,%.2f]"
                  % (total_out_of_bounds_pairs, critical_pairs,
                     CA_DIST_MIN, CA_DIST_MAX, CN_DIST_MIN, CN_DIST_MAX))
    else:
        status = STATUS_OK
        reason = "ALL_BACKBONE_BONDS_PHYSICAL_NRESIDUES=%d_CA=[%.1f,%.1f]_CN=[%.2f,%.2f]" % (
            len(residues), CA_DIST_MIN, CA_DIST_MAX, CN_DIST_MIN, CN_DIST_MAX)

    return {
        "status": status,
        "reason": reason,
        "n_residues": len(residues),
        "n_backbone_pairs_checked": len(residues) - 1,
        "out_of_bounds_pairs": out_of_bounds[:20],  # 上限 20 防止爆
        "n_out_of_bounds": total_out_of_bounds_pairs,
        "n_critical_pairs": critical_pairs,
        "critical_runs": critical_runs[:5],
        "thresholds": {
            "ca_dist_min_A": CA_DIST_MIN, "ca_dist_max_A": CA_DIST_MAX,
            "ca_dist_crit_min_A": CA_DIST_CRIT_MIN, "ca_dist_crit_max_A": CA_DIST_CRIT_MAX,
            "cn_dist_min_A": CN_DIST_MIN, "cn_dist_max_A": CN_DIST_MAX,
            "cn_dist_crit_max_A": CN_DIST_CRIT_MAX,
            "continuity_threshold": CRITICAL_CONTINUITY,
        },
        "probe": "probe_backbone_hallucination",
        "probe_fingerprint": "backbone_hallucination_v2_ca_cn_bond_length_continuity",
    }


# ========== 白盒自检 ==========

def run_self_checks():
    out = []

    def check(name, cond, detail=""):
        out.append({"name": name, "passed": bool(cond), "detail": detail})

    import tempfile

    def _write_pdb(path, residues_data):
        """residues_data: list of (resname, chain, resseq, [(atom_name, x, y, z), ...])。"""
        serial = 1
        with open(path, "w") as f:
            for resname, chain, resseq, atoms in residues_data:
                for aname, x, y, z in atoms:
                    f.write("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                            % (serial, aname, resname, chain, resseq, x, y, z, 1.00, 50.00, aname[0]))
                    serial += 1
            f.write("END\n")

    def _ideal_chain(n_res, start_x=0.0, spacing=3.8):
        """生成理想链：所有主链原子沿 x 轴，C_i 紧邻 N_{i+1}（肽键 1.33 Å）。
        CA_i = start + i * spacing
        N_i = CA_i - 1.33
        C_i = N_{i+1} - 1.33 = CA_i + spacing - 2.66
        """
        data = []
        for i in range(n_res):
            ca_x = start_x + i * spacing
            n_x = ca_x - 1.33
            if i < n_res - 1:
                # C_i = N_{i+1} - 1.33 = (CA_{i+1} - 1.33) - 1.33
                c_x = ca_x + spacing - 2.66
            else:
                # 最后一个残基：让 C_x = ca_x + 1.52（尾端无下一个 N）
                c_x = ca_x + 1.52
            o_x = c_x + 1.24  # C=O 键长 ~ 1.24 Å
            data.append(("ALA", "A", i + 1, [
                ("N", n_x, 0.0, 0.0),
                ("CA", ca_x, 0.0, 0.0),
                ("C", c_x, 0.0, 0.0),
                ("O", o_x, 0.0, 0.0),
            ]))
        return data

    # T1: 理想链 20 残基 → OK
    pdb_t1 = tempfile.mktemp(suffix=".pdb")
    _write_pdb(pdb_t1, _ideal_chain(20))
    r1 = scan_backbone_hallucination(pdb_t1)
    os.unlink(pdb_t1)
    check("T1_IDEAL_CHAIN_OK", r1["status"] == STATUS_OK,
          "status=%s reason=%s" % (r1["status"], r1["reason"][:80]))

    # T2: 连续 5 个残基 Cα 硬压到 3.0 Å（坍塌）→ CRITICAL
    pdb_t2 = tempfile.mktemp(suffix=".pdb")
    _write_pdb(pdb_t2, _ideal_chain(5, spacing=3.0))
    r2 = scan_backbone_hallucination(pdb_t2)
    os.unlink(pdb_t2)
    check("T2_CONTINUOUS_COLLAPSE_CRITICAL",
          r2["status"] == STATUS_HALLUCINATION,
          "status=%s n_critical_pairs=%d" % (r2["status"], r2.get("n_critical_pairs")))

    # T3: 孤立 1 个残基对越界（其他都正常）→ TENSION 而非 CRITICAL
    pdb_t3 = tempfile.mktemp(suffix=".pdb")
    data3 = _ideal_chain(10)
    # 偏移残基 index 4 (resseq=5) 的所有原子 +2.5 Å in y → 孤立拉断 pair 4→5
    data3[4] = (
        data3[4][0], data3[4][1], data3[4][2],
        [(a[0], a[1], a[2] + 2.5, a[3]) for a in data3[4][3]]
    )
    _write_pdb(pdb_t3, data3)
    r3 = scan_backbone_hallucination(pdb_t3)
    os.unlink(pdb_t3)
    check("T3_ISOLATED_BREAK_TENSION",
          r3["status"] == STATUS_TENSION,
          "status=%s (should be TENSION not HALLUCINATION)" % r3["status"])

    # T4: C-N 肽键连续拉断（C 放在 z=3.5 → C-N 距离 ≈ 3.7 Å）→ CRITICAL
    pdb_t4 = tempfile.mktemp(suffix=".pdb")
    data4 = _ideal_chain(5)
    data4_fixed = []
    for res in data4:
        new_atoms = []
        for aname, x, y, z in res[3]:
            if aname == "C":
                new_atoms.append((aname, x, 0.0, 3.5))
            else:
                new_atoms.append((aname, x, 0.0, 0.0))
        data4_fixed.append((res[0], res[1], res[2], new_atoms))
    _write_pdb(pdb_t4, data4_fixed)
    r4 = scan_backbone_hallucination(pdb_t4)
    os.unlink(pdb_t4)
    check("T4_CN_BOND_BREAK_CRITICAL",
          r4["status"] == STATUS_HALLUCINATION,
          "status=%s n_critical_pairs=%d" % (r4["status"], r4.get("n_critical_pairs")))

    # T5: 空文件 → EMPTY
    pdb_t5 = tempfile.mktemp(suffix=".pdb")
    with open(pdb_t5, "w") as f:
        f.write("HEADER TEST\nEND\n")
    r5 = scan_backbone_hallucination(pdb_t5)
    os.unlink(pdb_t5)
    check("T5_EMPTY_PDB_EMPTY", r5["status"] == STATUS_EMPTY)

    # T6: 双跑确定性
    pdb_t6 = tempfile.mktemp(suffix=".pdb")
    _write_pdb(pdb_t6, _ideal_chain(15))
    r6a = scan_backbone_hallucination(pdb_t6)
    r6b = scan_backbone_hallucination(pdb_t6)
    os.unlink(pdb_t6)
    check("T6_DETERMINISTIC", r6a["status"] == r6b["status"])

    # T7: 真实 PDB（1ubq）→ 应该是 OK 或 TENSION（天然蛋白有小扰动）
    real_pdb = os.path.join(HERE, "1ubq.pdb")
    if os.path.exists(real_pdb):
        r7 = scan_backbone_hallucination(real_pdb)
        check("T7_REAL_UBIQUITIN_NOT_HALLUCINATION",
              r7["status"] != STATUS_HALLUCINATION,
              "status=%s reason=%s" % (r7["status"], r7["reason"][:80]))
    else:
        check("T7_REAL_UBIQUITIN_NOT_HALLUCINATION", True, "skipped: 1ubq.pdb not found")

    # T8: 白盒纯度
    src = open(__file__, encoding="utf-8").read()
    purity = not re.search(r"(?m)^\s*(import|from)\s+(random|numpy|scipy|Bio)", src)
    check("T8_PURE_STDLIB_NO_BLACKBOX", purity)

    return out


if __name__ == "__main__":
    checks = run_self_checks()
    bad = [c for c in checks if not c["passed"]]
    print("probe_backbone_hallucination selfchecks total=%d passed=%d failed=%d"
          % (len(checks), len(checks) - len(bad), len(bad)))
    for c in bad:
        print("  FAIL %s :: %s" % (c["name"], c["detail"]))
