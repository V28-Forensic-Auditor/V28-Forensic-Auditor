# -*- coding: utf-8 -*-
"""probe_chirality_violation.py
TASK_20261110_V28_PHYSICS_PURIFIER —— 手性悖论探针（L7.6）。

职责：检测非甘氨酸残基 Cα 的立体构型是否为 D 型（右旋），即镜像废品。
  前置排除：忽略甘氨酸（GLY，无 Cβ，不具手性中心）。
  算法：取 N, Cα, C, Cβ 四原子坐标，计算标量三重积
        T = (N - Cα) · ((C - Cα) × (Cβ - Cα))（即三个向量的有向体积，单位 A^3）。
        符号标定（白盒实测，非查表）：天然 L 型蛋白 T > 0，
        实测 134l 119/119 全正、2KL8 75/75 全正；镜像结构 T < 0。
  击杀：T < 0 判为 D 型 -> [REJECTED_AI_STEREOCHEMICAL_PARADOX]

红线：纯 stdlib；只读坐标，不修改；确定性遍历，双跑一致。
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import op_full_atom_grid as fag

# 天然 L 型氨基酸的三重积符号（实测标定：L -> 正）
L_SIGN = 1
CHIRAL_ATOMS = ("N", "CA", "C", "CB")
SKIP_RESNAMES = ("GLY",)

STATUS_PARADOX = "[REJECTED_AI_STEREOCHEMICAL_PARADOX]"
STATUS_OK = "[OK_STEREOCHEMISTRY_VALID]"
STATUS_EMPTY = "[NOT_APPLICABLE_NO_CHIRAL_CENTER]"


def _triple_product(n, ca, c, cb):
    a = (n[0] - ca[0], n[1] - ca[1], n[2] - ca[2])
    b = (c[0] - ca[0], c[1] - ca[1], c[2] - ca[2])
    e = (cb[0] - ca[0], cb[1] - ca[1], cb[2] - ca[2])
    cx = b[1] * e[2] - b[2] * e[1]
    cy = b[2] * e[0] - b[0] * e[2]
    cz = b[0] * e[1] - b[1] * e[0]
    return a[0] * cx + a[1] * cy + a[2] * cz


def scan_chirality_violation(pdb_path, chain=None):
    """返回手性检测结果 dict。只读，不修改坐标。"""
    atoms = fag.parse_full_atoms(pdb_path, chain=chain)
    groups = {}
    order = []
    for a in atoms:
        key = (a["chain"], a["resseq"])
        if key not in groups:
            groups[key] = {}
            order.append(key)
        groups[key][a["name"]] = a

    checked = 0
    skipped_gly = 0
    d_residues = []
    for key in order:
        names = groups[key]
        resname = None
        for nm in CHIRAL_ATOMS:
            if nm in names:
                resname = names[nm]["resname"]
                break
        if resname is None:
            continue
        if resname in SKIP_RESNAMES:
            skipped_gly += 1
            continue
        if not all(nm in names for nm in CHIRAL_ATOMS):
            continue
        n = names["N"]["xyz"]
        ca = names["CA"]["xyz"]
        c = names["C"]["xyz"]
        cb = names["CB"]["xyz"]
        t = _triple_product(n, ca, c, cb)
        checked += 1
        if t * L_SIGN < 0.0:
            d_residues.append({
                "chain": key[0],
                "resseq": key[1],
                "resname": resname,
                "triple_product": round(t, 4),
                "assigned_chirality": "D",
            })

    if checked == 0:
        status = STATUS_EMPTY
    elif d_residues:
        status = STATUS_PARADOX
    else:
        status = STATUS_OK
    d_residues.sort(key=lambda r: (r["chain"] or "", r["resseq"] or ""))
    return {
        "probe": "probe_chirality_violation",
        "status": status,
        "triple_product_definition":
            "T = (N-CA) . ((C-CA) x (CB-CA)); L_SIGN=+1 (calibrated on 134l/2KL8)",
        "l_sign": L_SIGN,
        "n_chiral_centers_checked": checked,
        "n_glycine_skipped": skipped_gly,
        "n_d_residues": len(d_residues),
        "d_residues": d_residues,
    }


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "1ubq.pdb")
    r = scan_chirality_violation(p)
    print("input:", os.path.basename(p))
    print("status:", r["status"])
    print("checked:", r["n_chiral_centers_checked"],
          "gly_skipped:", r["n_glycine_skipped"],
          "D:", r["n_d_residues"])
    for d in r["d_residues"][:10]:
        print("  D %s%s  T=%.3f" % (d["resname"], d["resseq"], d["triple_product"]))
