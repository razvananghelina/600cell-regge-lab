"""
verify_s08_edge_fibration_uniformity.py
=======================================

Focused verification of the weak S08 claims on the full six-fibration class.

This script verifies exactly the following statements:

  Q1. There are exactly 6 discrete Hopf fibrations of the required left-coset
      type.
  Q2. For every such fibration F, the edge-space operator
      Box_1(F) = L_cross(F) - a1 * L_fiber(F)
      has
          dim ker(Box_1(F)) = 13
      and the full 2I-character decomposition
          ker(Box_1(F)) = rho_0 + 2 rho_5.
  Q3. For every such fibration F, the action of 2I on the 12 fiber labels
      factors through A5 = 2I / {+/-1}.
  Q4. For every such fibration F, the 12-dimensional fiber permutation module
      decomposes as
          1 + 3 + 3' + 5.

No claim is made here about any canonical identification between the 12-d
nontrivial edge-kernel sector and the 12-d fiber permutation module.
"""

import numpy as np
from scipy.linalg import eigh
from itertools import permutations, product as cartesian_product
from collections import defaultdict


a1 = 5
phi = (1.0 + np.sqrt(a1)) / 2.0
EDGE_TOL = 1e-6
TOL = 1e-8

N_PASS = 0
N_FAIL = 0


def check(condition, label, detail=""):
    global N_PASS, N_FAIL
    if condition:
        N_PASS += 1
        print(f"[PASS] {label}")
    else:
        N_FAIL += 1
        print(f"[FAIL] {label}")
    if detail:
        print(f"       {detail}")


def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def build_2I():
    verts = set()

    def add_vert(v):
        arr = np.array(v, dtype=float)
        n = np.linalg.norm(arr)
        if n > 1e-12:
            arr = arr / n
        verts.add(tuple(np.round(arr, 10)))

    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add_vert(v)

    for signs in cartesian_product([0.5, -0.5], repeat=4):
        add_vert(list(signs))

    base = [0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)]
    even_perms = []
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            even_perms.append(p)

    for perm in even_perms:
        coords = [base[perm[i]] for i in range(4)]
        nz_indices = [i for i in range(4) if abs(coords[i]) > 1e-12]
        for signs in cartesian_product([1, -1], repeat=len(nz_indices)):
            v = list(coords)
            for idx, s in zip(nz_indices, signs):
                v[idx] *= s
            add_vert(v)

    return np.array(sorted(verts))


def find_vertex_index(verts, q, tol=1e-6):
    dists = np.linalg.norm(verts - q, axis=1)
    idx = np.argmin(dists)
    if dists[idx] < tol:
        return idx
    return -1


def build_adjacency(verts):
    n = len(verts)
    dots = verts @ verts.T
    np.clip(dots, -1.0, 1.0, out=dots)
    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(dots[i, j] - phi / 2.0) < EDGE_TOL:
                A[i, j] = 1.0
                A[j, i] = 1.0
    return A


def find_all_hopf_fibrations(verts):
    n = len(verts)
    fibrations = []
    seen = set()

    order10_indices = []
    for i in range(n):
        g = verts[i]
        power = g.copy()
        for k in range(2, 121):
            power = quat_mult(power, g)
            if np.allclose(power, [1, 0, 0, 0], atol=1e-6):
                if k == 10:
                    order10_indices.append(i)
                break

    for gen_idx in order10_indices:
        g = verts[gen_idx]
        subgroup_indices = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        valid = True
        for _ in range(10):
            idx = find_vertex_index(verts, power)
            if idx < 0:
                valid = False
                break
            subgroup_indices.append(idx)
            power = quat_mult(power, g)
        if not valid or len(set(subgroup_indices)) != 10:
            continue

        assigned = np.full(n, -1, dtype=int)
        fib_list = []
        fid = 0
        for i in range(n):
            if assigned[i] >= 0:
                continue
            coset = []
            for si in subgroup_indices:
                q_prod = quat_mult(verts[i], verts[si])
                idx = find_vertex_index(verts, q_prod)
                if idx >= 0 and assigned[idx] < 0:
                    coset.append(idx)
                    assigned[idx] = fid
            if len(coset) != 10:
                valid = False
                break
            fib_list.append(tuple(sorted(coset)))
            fid += 1

        if not valid or fid != 12:
            continue

        signature = tuple(sorted(fib_list))
        if signature not in seen:
            seen.add(signature)
            fibrations.append(fib_list)

    return fibrations


def build_edges(A):
    n = A.shape[0]
    edges = []
    edge_to_idx = {}
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] > 0.5:
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))
    return edges, edge_to_idx


def build_line_graph_adj(n_vertices, edges):
    edge_inc = defaultdict(list)
    for e_idx, (i, j) in enumerate(edges):
        edge_inc[i].append(e_idx)
        edge_inc[j].append(e_idx)

    n_edges = len(edges)
    A_line = np.zeros((n_edges, n_edges), dtype=float)
    for v in range(n_vertices):
        inc = edge_inc[v]
        for a in range(len(inc)):
            for b in range(a + 1, len(inc)):
                ea = inc[a]
                eb = inc[b]
                A_line[ea, eb] = 1.0
                A_line[eb, ea] = 1.0
    return A_line


def build_fiber_edge_data(fib_list, A, edges, edge_to_idx):
    n_edges = len(edges)
    is_fiber = np.zeros(n_edges, dtype=bool)
    A_fiber_line = np.zeros((n_edges, n_edges), dtype=float)

    for fiber in fib_list:
        fiber_edges = []
        fiber_set = set(fiber)
        for i in fiber:
            for j in fiber:
                if i < j and A[i, j] > 0.5:
                    e_idx = edge_to_idx[(i, j)]
                    is_fiber[e_idx] = True
                    fiber_edges.append(e_idx)

        assert len(fiber_edges) == 10

        for a in range(len(fiber_edges)):
            i1, j1 = edges[fiber_edges[a]]
            for b in range(a + 1, len(fiber_edges)):
                i2, j2 = edges[fiber_edges[b]]
                if len({i1, j1} & {i2, j2}) == 1 and {i1, j1, i2, j2} <= fiber_set:
                    A_fiber_line[fiber_edges[a], fiber_edges[b]] = 1.0
                    A_fiber_line[fiber_edges[b], fiber_edges[a]] = 1.0

    return is_fiber, A_fiber_line


def build_group_actions(verts, edges, edge_to_idx):
    n = len(verts)
    all_vperm = []
    all_eperm = []
    for g_idx in range(n):
        vp = np.array([
            find_vertex_index(verts, quat_mult(verts[g_idx], verts[i]))
            for i in range(n)
        ])
        all_vperm.append(vp)
        ep = np.array([
            edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))]
            for i, j in edges
        ])
        all_eperm.append(ep)
    return all_vperm, all_eperm


def build_irrep_characters(A, all_vperm):
    lap = np.diag(np.sum(A, axis=1)) - A
    evals_lap, evecs_lap = eigh(lap)

    ev_sp = {}
    for i, val in enumerate(evals_lap):
        key = round(val, 8)
        ev_sp.setdefault(key, []).append(i)
    eig_list = sorted(ev_sp.items())

    char_vecs = np.zeros((len(all_vperm), len(eig_list)))
    for g in range(len(all_vperm)):
        for k, (_, indices) in enumerate(eig_list):
            space = evecs_lap[:, indices]
            char_vecs[g, k] = sum(
                np.dot(space[:, c], space[:, c][all_vperm[g]])
                for c in range(space.shape[1])
            )

    class_dict = defaultdict(list)
    for g in range(len(all_vperm)):
        class_dict[tuple(np.round(char_vecs[g], 8))].append(g)

    class_list = sorted(class_dict.items(), key=lambda x: (-len(x[1]), x[0]))
    class_sizes = np.array([len(members) for _, members in class_list], dtype=float)
    irrep_dims = [int(round(np.sqrt(len(indices)))) for _, indices in eig_list]

    chi_true = np.zeros((len(eig_list), len(class_list)))
    for ri in range(len(eig_list)):
        for ci, (key, _) in enumerate(class_list):
            chi_true[ri, ci] = key[ri] / irrep_dims[ri]

    return class_list, class_sizes, irrep_dims, chi_true


def rep_character_on_class(basis, perm):
    return sum(np.dot(basis[:, k], basis[:, k][perm]) for k in range(basis.shape[1]))


def decompose_character(rep_class_chars, class_sizes, chi_true):
    mults = {}
    group_order = np.sum(class_sizes)
    for ri in range(chi_true.shape[0]):
        n = np.sum(class_sizes * chi_true[ri] * rep_class_chars) / group_order
        if abs(n) > 0.01:
            mults[ri] = int(round(n))
    return mults


def main():
    print("=" * 72)
    print("VERIFY S08 EDGE/FIBER UNIFORMITY")
    print("=" * 72)

    verts = build_2I()
    A = build_adjacency(verts)
    edges, edge_to_idx = build_edges(A)
    A_line = build_line_graph_adj(len(verts), edges)
    all_vperm, all_eperm = build_group_actions(verts, edges, edge_to_idx)
    class_list, class_sizes, irrep_dims, chi_true = build_irrep_characters(A, all_vperm)

    minus1 = None
    for g in range(len(verts)):
        if np.allclose(verts[g], [-1, 0, 0, 0], atol=1e-6):
            minus1 = g
            break
    assert minus1 is not None

    fibrations = find_all_hopf_fibrations(verts)
    check(len(fibrations) == 6, "Q1: exactly 6 discrete Hopf fibrations",
          f"Found {len(fibrations)}")

    all_kernel_dim13 = True
    all_kernel_decomp = True
    all_factor_through_a5 = True
    all_perm_decomp = True

    expected_kernel = {0: 1, 5: 2}
    expected_perm = {0: 1, 2: 1, 4: 1, 6: 1}

    print()
    for idx, fib_list in enumerate(fibrations, start=1):
        is_fiber, A_fiber_line = build_fiber_edge_data(fib_list, A, edges, edge_to_idx)
        A_cross_line = A_line - A_fiber_line
        L_fiber = np.diag(np.sum(A_fiber_line, axis=1)) - A_fiber_line
        L_cross = np.diag(np.sum(A_cross_line, axis=1)) - A_cross_line
        Box1 = L_cross - a1 * L_fiber

        evals_box, evecs_box = eigh(Box1)
        ker_mask = np.abs(evals_box) < TOL
        ker_basis = evecs_box[:, ker_mask]
        if ker_basis.shape[1] != 13:
            all_kernel_dim13 = False

        ker_class_chars = np.zeros(len(class_list))
        for ci, (_, members) in enumerate(class_list):
            rep = members[0]
            ker_class_chars[ci] = rep_character_on_class(ker_basis, all_eperm[rep])
        kernel_decomp = decompose_character(ker_class_chars, class_sizes, chi_true)
        if kernel_decomp != expected_kernel:
            all_kernel_decomp = False

        vtf = {}
        for fi, fiber in enumerate(fib_list):
            for v in fiber:
                vtf[v] = fi

        minus1_trivial = True
        for fi, fiber in enumerate(fib_list):
            v_rep = fiber[0]
            gv = all_vperm[minus1][v_rep]
            if vtf[gv] != fi:
                minus1_trivial = False
                break
        if not minus1_trivial:
            all_factor_through_a5 = False

        perm_class_chars = np.zeros(len(class_list))
        for ci, (_, members) in enumerate(class_list):
            rep = members[0]
            n_fixed = 0
            for fi, fiber in enumerate(fib_list):
                v_rep = fiber[0]
                if vtf[all_vperm[rep][v_rep]] == fi:
                    n_fixed += 1
            perm_class_chars[ci] = n_fixed
        perm_decomp = decompose_character(perm_class_chars, class_sizes, chi_true)
        if perm_decomp != expected_perm:
            all_perm_decomp = False

        print(
            f"Fibration {idx}: "
            f"ker_dim={ker_basis.shape[1]}, "
            f"ker_decomp={kernel_decomp}, "
            f"minus1_trivial={minus1_trivial}, "
            f"perm_decomp={perm_decomp}"
        )

    print()
    check(all_kernel_dim13, "Q2a: dim ker(Box_1(F)) = 13 for all 6 fibrations")
    check(
        all_kernel_decomp,
        "Q2b: ker(Box_1(F)) = rho_0 + 2 rho_5 for all 6 fibrations",
    )
    check(
        all_factor_through_a5,
        "Q3: the fiber action factors through A5 for all 6 fibrations",
    )
    check(
        all_perm_decomp,
        "Q4: fiber permutation module = 1 + 3 + 3' + 5 for all 6 fibrations",
    )

    print()
    print(f"TOTAL: {N_PASS} passed, {N_FAIL} failed")
    return 0 if N_FAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
