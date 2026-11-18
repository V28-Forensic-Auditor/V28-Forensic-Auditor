# -*- coding: utf-8 -*-
"""op_full_atom_grid.py
TASK_20261102_V28_MILGRADE_PHYSICS_AUDIT —— 全原子加载与空间哈希网格。

职责：
  1. 解析 PDB，提取全部非氢原子（HEAVY ATOMS），带元素与 vdW 半径。
  2. 硬编码 vdW 半径字典（不查外部表，不黑盒）。
  3. 三维空间哈希网格，边长固定 4.0 A，把邻域查询从 O(N^2) 降到 O(N)。
  4. 单帧非氢原子数 <= 2500 的算力防爆闸；超出则确定性降采样并返回
     [WARNING_DOWNSAMPLED_FOR_PERFORMANCE]。

红线：纯 stdlib（math, os, sys）；零随机；所有输出经 round(...,6) 固定。
"""
import math
import os

# 硬编码 vdW 半径 (A)，来源：Bondi/Alvarez 常用取值，冻结在本文件内。
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98,
    "Se": 1.90, "Zn": 1.39, "Fe": 1.32, "Mg": 1.73, "Mn": 1.39,
    "Cu": 1.40, "Ni": 1.63, "Co": 1.26, "Na": 2.27, "K": 2.75,
}
DEFAULT_VDW = 1.70
MAX_VDW = max(VDW_RADII.values())

# 共价半径 (A)：仅用于建共价键图（不能用 vdW，否则会把非键接触误判成键）。
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05,
    "P": 1.07, "F": 0.57, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Se": 1.20, "Zn": 1.22, "Fe": 1.32, "Mg": 1.41, "Mn": 1.39,
    "Cu": 1.32, "Ni": 1.24, "Co": 1.26, "Na": 1.66, "K": 2.03,
}
DEFAULT_COVALENT = 0.76
MAX_COVALENT = max(COVALENT_RADII.values())

CELL_SIZE = 4.0              # 空间哈希网格边长 (A)
MAX_ATOMS_PER_FRAME = 2500   # 单帧非氢原子硬上限（纯 Python 防 OOM）
BOND_MARGIN = 0.45           # 共价键判据余量 (A)

WARNING_DOWNSAMPLED = "[WARNING_DOWNSAMPLED_FOR_PERFORMANCE]"


# ----------------------------- PDB 解析 -----------------------------

def normalize_element(raw, atom_name):
    """元素符号归一化：'ZN'->'Zn'，'SE'->'Se'。"""
    e = (raw or "").strip()
    if e:
        return e[0].upper() + e[1:].lower() if len(e) > 1 else e.upper()
    # 无元素列时从原子名推导（PDB 规范：原子名前两字符，去数字）
    name = "".join(ch for ch in atom_name if ch.isalpha())
    return name[:1].upper() if name else "C"


def parse_full_atoms(pdb_path, chain="A", keep_hetatm=True):
    """解析 PDB，返回重原子列表（按文件顺序，确定性）。

    每个元素：{"index", "name", "element", "vdw", "resname", "chain",
               "resseq", "xyz": (x,y,z)}
    跳过氢原子；跳过 altLoc 非空且非 'A' 的重复构象记录。
    """
    atoms = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            is_atom = line.startswith("ATOM")
            is_het = line.startswith("HETATM")
            if not (is_atom or (keep_hetatm and is_het)):
                continue
            altloc = line[16:17]
            if altloc not in (" ", "", "A"):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            ch = line[21].strip() if len(line) > 21 else ""
            if chain is not None and ch != chain:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = normalize_element(line[76:78] if len(line) >= 78 else "",
                                        atom_name)
            if element == "H":
                continue
            vdw = VDW_RADII.get(element, DEFAULT_VDW)
            cov = COVALENT_RADII.get(element, DEFAULT_COVALENT)
            atoms.append({
                "index": len(atoms),
                "name": atom_name,
                "element": element,
                "vdw": vdw,
                "cov": cov,
                "resname": resname,
                "chain": ch,
                "resseq": line[22:26].strip(),
                "xyz": (round(x, 6), round(y, 6), round(z, 6)),
            })
    return atoms


def atom_coords(atoms):
    return [a["xyz"] for a in atoms]


def atom_radii(atoms):
    return [a["vdw"] for a in atoms]


def atom_covalent(atoms):
    return [a["cov"] for a in atoms]


def enforce_atom_cap(atoms, cap=MAX_ATOMS_PER_FRAME):
    """算力防爆闸。返回 (atoms_used, warning or None, stride)。"""
    n = len(atoms)
    if n <= cap:
        return atoms, None, 1
    stride = (n + cap - 1) // cap
    return atoms[::stride], WARNING_DOWNSAMPLED, stride


# ----------------------------- 全原子拓扑点云选择 -----------------------------
#
# 红线：拓扑持久同调必须与位阻检测处在同一物理尺度，废弃「只取 Cα 骨架」
# 的 stride=2 点云。每残基选点规则：
#   - 主链 N, CA, C, O（+ OXT 若存在）
#   - 侧链 CB（若有）
#   - 侧链末端重原子（离 CA 最远的侧链重原子，避免长侧链被截断）
# 得到约 5 点/残基（GLY 约 4 点）。若超过 max_points，按体素网格自适应
# 降采样（每格保留一个原子，确定性），并返回 WARNING_TOPOLOGY_POINT_REDUCED。

TOPOLOGY_CANDIDATE_NAMES = ("N", "CA", "C", "O", "OXT")
TOPOLOGY_BASE_CELL = 1.5        # 降采样起始体素边长 (A)
WARNING_TOPOLOGY_POINT_REDUCED = "[WARNING_TOPOLOGY_POINT_REDUCED]"


def _residue_groups(atoms):
    """按 (chain, resseq) 分组，保持 PDB 顺序。"""
    by_res = {}
    order = []
    for a in atoms:
        key = (a["chain"], a["resseq"])
        if key not in by_res:
            by_res[key] = []
            order.append(key)
        by_res[key].append(a)
    return [(key, by_res[key]) for key in order]


def select_topology_candidates(atoms):
    """返回拓扑候选点原子索引（确定性，升序）。"""
    out = []
    for _key, res in _residue_groups(atoms):
        ca = None
        for a in res:
            if a["name"] == "CA":
                ca = a
                break
        for a in res:
            if a["name"] in TOPOLOGY_CANDIDATE_NAMES:
                out.append(a["index"])
        for a in res:
            if a["name"] == "CB":
                out.append(a["index"])
                break
        if ca is not None:
            ca_xyz = ca["xyz"]
            side = [a for a in res
                    if a["name"] not in TOPOLOGY_CANDIDATE_NAMES
                    and a["name"] != "CB"
                    and a["element"] != "H"]
            if side:
                def _d2(a):
                    p = a["xyz"]
                    return ((p[0] - ca_xyz[0]) ** 2 + (p[1] - ca_xyz[1]) ** 2
                            + (p[2] - ca_xyz[2]) ** 2)
                far = max(side, key=_d2)
                out.append(far["index"])
    return sorted(set(out))


def select_topology_points(atoms, max_points=128, base_cell=TOPOLOGY_BASE_CELL):
    """全原子降采样点云（拓扑持久同调输入）。

    返回 dict：indices / coords / point_count / candidate_count /
               voxel_cell / downsampled / warning / ratio。
    """
    coords_full = atom_coords(atoms)
    cand = select_topology_candidates(atoms)
    n_cand = len(cand)
    cell = float(base_cell)
    kept = list(cand)
    down_steps = 0
    while len(kept) > max_points and cell < 200.0:
        seen = {}
        order = []
        for idx in cand:
            p = coords_full[idx]
            key = (int(math.floor(p[0] / cell)),
                   int(math.floor(p[1] / cell)),
                   int(math.floor(p[2] / cell)))
            if key in seen:
                continue
            seen[key] = idx
            order.append(idx)
        kept = order
        down_steps += 1
        if len(kept) <= max_points:
            break
        cell *= 1.15
    warning = (WARNING_TOPOLOGY_POINT_REDUCED
               if len(kept) < n_cand else None)
    return {
        "indices": kept,
        "coords": [coords_full[i] for i in kept],
        "radii": [atoms[i]["vdw"] for i in kept],
        "point_count": len(kept),
        "candidate_count": n_cand,
        "atom_count": len(atoms),
        "voxel_cell": round(cell, 6),
        "downsample_steps": down_steps,
        "downsampled": warning is not None,
        "warning": warning,
        "ratio": round(len(kept) / float(n_cand), 6) if n_cand else 0.0,
    }


# ----------------------------- 空间哈希网格 -----------------------------

class SpatialHashGrid(object):
    """三维空间哈希网格：cell_size=4.0 A，邻域查询 O(1) 均摊。"""

    __slots__ = ("cell", "cells")

    def __init__(self, cell=CELL_SIZE):
        self.cell = float(cell)
        self.cells = {}

    def _key(self, p):
        c = self.cell
        return (int(math.floor(p[0] / c)),
                int(math.floor(p[1] / c)),
                int(math.floor(p[2] / c)))

    def insert(self, idx, p):
        k = self._key(p)
        bucket = self.cells.get(k)
        if bucket is None:
            self.cells[k] = [idx]
        else:
            bucket.append(idx)

    def build(self, coords):
        for i, p in enumerate(coords):
            self.insert(i, p)
        return self

    def neighbor_keys(self, p, radius):
        """返回覆盖 radius 球的所有网格键（确定性顺序：z->y->x）。"""
        c = self.cell
        span = int(math.ceil(radius / c))
        kx, ky, kz = self._key(p)
        keys = []
        for dz in range(-span, span + 1):
            for dy in range(-span, span + 1):
                for dx in range(-span, span + 1):
                    keys.append((kx + dx, ky + dy, kz + dz))
        return keys

    def query(self, p, radius, coords):
        """返回 (jid, dist, dist2) 列表，仅含 radius 内原子（确定性顺序）。"""
        r2 = radius * radius
        out = []
        for k in self.neighbor_keys(p, radius):
            bucket = self.cells.get(k)
            if not bucket:
                continue
            for j in bucket:
                q = coords[j]
                dx = p[0] - q[0]
                dy = p[1] - q[1]
                dz = p[2] - q[2]
                d2 = dx * dx + dy * dy + dz * dz
                if d2 <= r2:
                    out.append((j, math.sqrt(d2), d2))
        return out


def build_grid(coords, cell=CELL_SIZE):
    g = SpatialHashGrid(cell)
    g.build(coords)
    return g


# ----------------------------- 共价键图 -----------------------------

def build_bonds(coords, cov_radii, cell=CELL_SIZE, margin=BOND_MARGIN):
    """按共价半径判据建共价键图：d < covA + covB + margin。

    返回 (bond_set, adjacency)：
      bond_set: {(i,j)}，i<j
      adjacency: [[j,...], ...] 按索引升序
    """
    grid = build_grid(coords, cell)
    n = len(coords)
    adjacency = [[] for _ in range(n)]
    bond_set = set()
    for i in range(n):
        cutoff = cov_radii[i] + MAX_COVALENT + margin
        for j, d, _d2 in grid.query(coords[i], cutoff, coords):
            if j <= i:
                continue
            if d < cov_radii[i] + cov_radii[j] + margin:
                bond_set.add((i, j))
                adjacency[i].append(j)
                adjacency[j].append(i)
    for a in adjacency:
        a.sort()
    return bond_set, adjacency


def excluded_pair_set(bond_set, adjacency, max_graph_distance=3):
    """按共价键图距离 <= 3 排除配对：1-2 成键、1-3（键角）、1-4（二面角）。

    位阻/能量判定必须跳过这些配对，否则会把正常键长/键角/环内几何
    当作碰撞。1-4 在共价几何上被二面角约束，非自由非键接触，故一并排除
    （这是 clash 检测的通行口径；排除范围在报告中如实披露）。
    """
    n = len(adjacency)
    excl = set()
    for i in range(n):
        frontier = [i]
        visited = {i}
        depth = 0
        while depth < max_graph_distance and frontier:
            nxt = []
            for u in frontier:
                for v in adjacency[u]:
                    if v not in visited:
                        visited.add(v)
                        nxt.append(v)
            frontier = nxt
            depth += 1
        for v in visited:
            if v != i:
                excl.add((i, v) if i < v else (v, i))
    return excl


# ----------------------------- 自检 -----------------------------

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    pdb = os.path.join(here, "23WN.pdb")
    if not os.path.exists(pdb):
        print("23WN.pdb not found")
        raise SystemExit(1)
    atoms = parse_full_atoms(pdb, chain="A")
    used, warn, stride = enforce_atom_cap(atoms)
    elements = {}
    for a in used:
        elements[a["element"]] = elements.get(a["element"], 0) + 1
    print("heavy atoms (chain A):", len(atoms), "-> used:", len(used),
          "stride:", stride, "warn:", warn)
    print("element counts:", dict(sorted(elements.items())))
    coords = atom_coords(used)
    radii = atom_radii(used)
    cov = atom_covalent(used)
    g = build_grid(coords)
    print("grid cells:", len(g.cells), "cell size:", g.cell)
    bonds, adj = build_bonds(coords, cov)
    excl = excluded_pair_set(bonds, adj)
    print("bonds:", len(bonds), "excluded(1-2/1-3):", len(excl))
    # 邻域查询一致性抽查：暴力距离 vs 网格查询
    probe = coords[0]
    r = 6.0
    brute = [(j, math.dist(probe, coords[j])) for j in range(len(coords))]
    brute_in = {j for j, d in brute if d <= r}
    grid_in = {j for j, _d, _d2 in g.query(probe, r, coords)}
    print("grid query consistency:", brute_in == grid_in,
          "| neighbors:", len(grid_in))
