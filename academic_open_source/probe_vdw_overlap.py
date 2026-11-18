# -*- coding: utf-8 -*-
"""probe_vdw_overlap.py
TASK_20261110_V28_PHYSICS_PURIFIER —— 范德华穿模探针（L7.6）。

职责：检测非键合重原子（C/N/O/S）之间的极端空间重叠。
  判据：非键合重原子距离 < 1.6 A -> [REJECTED_AI_VAN_DER_WAALS_OVERLAP]

硬编码前置排除（全部白盒，无外部表）：
  1) 金属原子：按元素列排除 Zn/Mg/Ca/Fe/Mn/Cu/Ni/Co/Na/K 等；
     只保留 C/N/O/S 四类重原子参与比较（Cα 元素为 C，不会被误当钙）。
  2) S-S 二硫键：CYS 的 S 原子对距离 < 2.5 A 视为二硫键，整对排除。
  3) 相邻残基：同链且 |Δresseq| < 2（同残基+相邻残基）视为键合邻域，
     不参与穿模比较（避免肽键/1-2/1-3 合法短距误判）。

关键设计（防止把冲突误判成键）：不得用"距离 < 共价半径和"来建键图——
否则任何 <1.6 A 的冲突对都会被当成共价键而排除，探头永不触发。
因此本探针只用「残基邻接」这一拓扑先验划出键合邻域，非相邻残基对
默认视为非键合，直接做距离比较。

红线：纯 stdlib；只读坐标，绝不修改任何坐标；确定性遍历，双跑一致。
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import op_full_atom_grid as fag

VDW_HARD_MIN = 1.6                      # 极端重叠阈值 (A)
DISULFIDE_MAX = 2.5                     # S-S 二硫键判定上限 (A)
TARGET_ELEMENTS = ("C", "N", "O", "S")  # 参与比较的重原子
METAL_ELEMENTS = {
    "Zn", "Mg", "Ca", "Fe", "Mn", "Cu", "Ni", "Co", "Na", "K",
    "Cd", "Hg", "Sr", "Cs", "Ba", "Mo", "W", "Pt", "Au", "Ag",
}
ADJACENCY_SKIP = 2                      # 同链 |Δresseq| < 2 视为相邻，跳过

STATUS_OVERLAP = "[REJECTED_AI_VAN_DER_WAALS_OVERLAP]"
STATUS_OK = "[OK_NO_VDW_OVERLAP]"
STATUS_EMPTY = "[NOT_APPLICABLE_NO_ATOMS]"


def _resseq_int(raw):
    s = (raw or "").strip()
    digits = ""
    for ch in s:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    try:
        return int(digits)
    except ValueError:
        return None


def _element(line, atom_name):
    raw = line[76:78] if len(line) >= 78 else ""
    return fag.normalize_element(raw, atom_name)


def parse_target_atoms(pdb_path, chain=None):
    """解析并只保留 C/N/O/S 重原子，附带原始 PDB 原子序号。"""
    atoms = []
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[16:17] not in (" ", "", "A"):
                continue
            ch = line[21].strip() if len(line) > 21 else ""
            if chain is not None and ch != chain:
                continue
            atom_name = line[12:16].strip()
            element = _element(line, atom_name)
            if element not in TARGET_ELEMENTS:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            try:
                serial = int(line[6:11])
            except ValueError:
                serial = 0
            atoms.append({
                "index": len(atoms),
                "serial": serial,
                "name": atom_name,
                "element": element,
                "resname": line[17:20].strip(),
                "chain": ch,
                "resseq": _resseq_int(line[22:26]),
                "resseq_raw": line[22:26].strip(),
                "xyz": (round(x, 6), round(y, 6), round(z, 6)),
            })
    return atoms


def find_disulfide_pairs(atoms, coords):
    """CYS 的 S 原子对距离 < DISULFIDE_MAX -> 二硫键配对索引集合。"""
    s_idx = [a["index"] for a in atoms
             if a["element"] == "S" and a["resname"] == "CYS"]
    pairs = set()
    for ii in range(len(s_idx)):
        for jj in range(ii + 1, len(s_idx)):
            i, j = s_idx[ii], s_idx[jj]
            if atoms[i]["chain"] == atoms[j]["chain"] and \
               atoms[i]["resseq"] is not None and \
               atoms[i]["resseq"] == atoms[j]["resseq"]:
                continue
            p, q = coords[i], coords[j]
            d = math.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 +
                          (p[2] - q[2]) ** 2)
            if d < DISULFIDE_MAX:
                pairs.add((min(i, j), max(i, j)))
    return pairs


def _is_adjacent(a, b):
    """同链且 |Δresseq| < ADJACENCY_SKIP 视为键合邻域。"""
    if a["chain"] != b["chain"]:
        return False
    if a["resseq"] is None or b["resseq"] is None:
        return False
    return abs(a["resseq"] - b["resseq"]) < ADJACENCY_SKIP


def scan_vdw_overlap(pdb_path, chain=None):
    """返回穿模检测结果 dict。只读，不写回任何坐标。"""
    atoms = parse_target_atoms(pdb_path, chain=chain)
    if not atoms:
        return {
            "probe": "probe_vdw_overlap",
            "status": STATUS_EMPTY,
            "n_atoms_compared": 0,
            "n_overlaps": 0,
            "overlaps": [],
        }
    coords = [a["xyz"] for a in atoms]
    ss_pairs = find_disulfide_pairs(atoms, coords)

    grid = fag.build_grid(coords)
    overlaps = []
    for i in range(len(coords)):
        ai = atoms[i]
        for (j, d, _d2) in grid.query(coords[i], VDW_HARD_MIN, coords):
            if j <= i:
                continue
            if (i, j) in ss_pairs:
                continue
            aj = atoms[j]
            if _is_adjacent(ai, aj):
                continue
            if d < VDW_HARD_MIN:
                overlaps.append({
                    "atom_i": {"serial": ai["serial"], "name": ai["name"],
                               "resname": ai["resname"],
                               "chain": ai["chain"],
                               "resseq": ai["resseq_raw"]},
                    "atom_j": {"serial": aj["serial"], "name": aj["name"],
                               "resname": aj["resname"],
                               "chain": aj["chain"],
                               "resseq": aj["resseq_raw"]},
                    "element_pair": ai["element"] + "-" + aj["element"],
                    "distance_A": round(d, 4),
                })
    overlaps.sort(key=lambda r: (r["distance_A"], r["atom_i"]["serial"],
                                 r["atom_j"]["serial"]))
    return {
        "probe": "probe_vdw_overlap",
        "status": STATUS_OVERLAP if overlaps else STATUS_OK,
        "threshold_A": VDW_HARD_MIN,
        "n_atoms_compared": len(atoms),
        "n_disulfide_pairs_excluded": len(ss_pairs),
        "n_overlaps": len(overlaps),
        "overlaps": overlaps,
        "exclusions": {
            "metals_by_element": sorted(METAL_ELEMENTS),
            "target_elements": list(TARGET_ELEMENTS),
            "disulfide_max_A": DISULFIDE_MAX,
            "adjacency_skip_delta_res": ADJACENCY_SKIP,
            "bond_graph_used": False,
            "note": "仅用残基邻接划键合邻域，禁止用距离建键以免吞掉冲突对",
        },
    }


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "1ubq.pdb")
    r = scan_vdw_overlap(p)
    print("input:", os.path.basename(p))
    print("status:", r["status"])
    print("atoms compared:", r.get("n_atoms_compared"),
          "overlaps:", r.get("n_overlaps"))
    for o in r.get("overlaps", [])[:10]:
        print("  %s%d-%s  %s%d-%s  %s  %.3f A"
              % (o["atom_i"]["resname"], o["atom_i"]["serial"],
                 o["atom_i"]["name"], o["atom_j"]["resname"],
                 o["atom_j"]["serial"], o["atom_j"]["name"],
                 o["element_pair"], o["distance_A"]))
