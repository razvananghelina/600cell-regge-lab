"""
exp608: Candidate Lie brackets on the 12 gauge modes.

Goal:
  Attack the remaining gauge-algebra gap directly.

Two natural candidates are tested:
  1. The simplicial wedge bracket on edge 1-forms, projected back to edges.
  2. An oriented triangle bracket on the 12-node fiber graph (icosahedron).

If both fail, the result supports the current interpretation:
the discrete theory derives the 12-dimensional gauge skeleton, while the
continuous Lie bracket belongs to the continuum completion.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh, norm

sys.path.insert(0, ".")
from commons import build_600cell

PHI = (1.0 + math.sqrt(5.0)) / 2.0
A1 = 5.0
TOL = 1e-8


def qmul(p, q):
    return np.array(
        [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
        ]
    )


def find_idx(v, verts, tol=1e-6):
    dots = verts @ v
    idx = np.argmax(dots)
    return idx if dots[idx] > 1.0 - tol else -1


def hopf_map(q):
    """Quaternion q = (w, x, y, z) -> S^2."""
    w, x, y, z = q
    return np.array(
        [
            2.0 * (x * z + w * y),
            2.0 * (y * z - w * x),
            w * w + z * z - x * x - y * y,
        ]
    )


def build_complex():
    verts, adj, lap = build_600cell()
    nv = len(verts)

    adj_list = defaultdict(set)
    edges = []
    edge_to_idx = {}
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                adj_list[i].add(j)
                adj_list[j].add(i)
                edge_to_idx[(i, j)] = len(edges)
                edges.append((i, j))

    triangles = []
    for i in range(nv):
        for j in adj_list[i]:
            if j > i:
                for k in adj_list[i] & adj_list[j]:
                    if k > j:
                        triangles.append((i, j, k))

    d0 = np.zeros((len(edges), nv))
    for e_idx, (i, j) in enumerate(edges):
        d0[e_idx, i] = -1.0
        d0[e_idx, j] = +1.0

    d1 = np.zeros((len(triangles), len(edges)))
    tri_edges = []
    for f_idx, (i, j, k) in enumerate(triangles):
        e_ij = edge_to_idx[(i, j)]
        e_jk = edge_to_idx[(j, k)]
        e_ik = edge_to_idx[(i, k)]
        d1[f_idx, e_ij] = +1.0
        d1[f_idx, e_jk] = +1.0
        d1[f_idx, e_ik] = -1.0
        tri_edges.append((e_ij, e_jk, e_ik))

    return verts, adj, lap, edges, edge_to_idx, triangles, tri_edges, d0, d1


def find_fibration(verts):
    nv = len(verts)
    target_w = PHI / 2.0
    for i in range(nv):
        if abs(verts[i, 0] - target_w) < 1e-6:
            g = verts[i]
            p = g.copy()
            ok = True
            for k in range(2, 11):
                p = qmul(p, g)
                if k == 5 and not np.allclose(p, [-1, 0, 0, 0], atol=1e-6):
                    ok = False
                    break
                if k == 10 and not np.allclose(p, [1, 0, 0, 0], atol=1e-6):
                    ok = False
            if not ok:
                continue
            used = set()
            fibers = []
            subgroup = []
            pp = np.array([1.0, 0, 0, 0])
            for _ in range(10):
                subgroup.append(find_idx(pp, verts))
                pp = qmul(pp, g)
            for s in range(nv):
                if s in used:
                    continue
                fib = []
                for si in subgroup:
                    idx = find_idx(qmul(verts[s], verts[si]), verts)
                    if idx >= 0 and idx not in used:
                        fib.append(idx)
                        used.add(idx)
                if len(fib) == 10:
                    fibers.append(fib)
            if len(fibers) == 12:
                return fibers
    raise RuntimeError("Could not find Hopf fibration")


def build_edge_box(edges, triangles, d0, d1, fibers, edge_to_idx):
    ne = len(edges)
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    is_fiber_edge = np.zeros(ne, dtype=bool)
    for e_idx, (i, j) in enumerate(edges):
        if vertex_to_fiber[i] == vertex_to_fiber[j]:
            is_fiber_edge[e_idx] = True

    vte = defaultdict(list)
    for e_idx, (i, j) in enumerate(edges):
        vte[i].append(e_idx)
        vte[j].append(e_idx)

    a_line = np.zeros((ne, ne))
    for v, inc in vte.items():
        _ = v
        for a in range(len(inc)):
            for b in range(a + 1, len(inc)):
                ea = inc[a]
                eb = inc[b]
                a_line[ea, eb] = 1.0
                a_line[eb, ea] = 1.0

    a_fib_line = np.zeros((ne, ne))
    fiber_edge_lists = []
    for fiber in fibers:
        fe = []
        for k in range(10):
            u = fiber[k]
            v = fiber[(k + 1) % 10]
            e_idx = edge_to_idx[(min(u, v), max(u, v))]
            fe.append(e_idx)
        fiber_edge_lists.append(fe)
        for k in range(10):
            ea = fe[k]
            eb = fe[(k + 1) % 10]
            a_fib_line[ea, eb] = 1.0
            a_fib_line[eb, ea] = 1.0

    lf = np.diag(np.sum(a_fib_line, axis=1)) - a_fib_line
    lc = np.diag(np.sum(a_line - a_fib_line, axis=1)) - (a_line - a_fib_line)
    box1 = lc - A1 * lf
    return box1, is_fiber_edge, vertex_to_fiber, fiber_edge_lists


def build_group_edge_perms(verts, edges, edge_to_idx):
    all_eperm = []
    for g_idx in range(len(verts)):
        vp = np.array([find_idx(qmul(verts[g_idx], verts[i]), verts) for i in range(len(verts))])
        ep = np.array([edge_to_idx[(min(vp[i], vp[j]), max(vp[i], vp[j]))] for i, j in edges])
        all_eperm.append(ep)
    return all_eperm


def build_gauge_basis(box1, all_eperm):
    evals, evecs = eigh(box1)
    ker = evecs[:, np.abs(evals) < 1e-7]
    if ker.shape[1] != 13:
        raise RuntimeError(f"Expected 13 kernel vectors, found {ker.shape[1]}")

    p0_ker = np.zeros((13, 13))
    for ep in all_eperm:
        for ki in range(13):
            gv = ker[:, ki][ep]
            for kj in range(13):
                p0_ker[ki, kj] += np.dot(ker[:, kj], gv)
    p0_ker /= float(len(all_eperm))

    p5_ker = np.eye(13) - p0_ker
    p5e, p5v = eigh(p5_ker)
    rho5_basis_ker = p5v[:, p5e > 0.5]
    gauge = ker @ rho5_basis_ker
    return gauge, p5e, ker


def wedge_to_faces(u, v, tri_edges):
    out = np.zeros(len(tri_edges))
    for t_idx, (e1, e2, e3) in enumerate(tri_edges):
        out[t_idx] = (
            u[e1] * v[e2]
            - u[e2] * v[e1]
            + u[e2] * v[e3]
            - u[e3] * v[e2]
            + u[e3] * v[e1]
            - u[e1] * v[e3]
        )
    return out


def build_fiber_graph(edges, is_fiber_edge, vertex_to_fiber):
    w = np.zeros((12, 12))
    for e_idx, (i, j) in enumerate(edges):
        if is_fiber_edge[e_idx]:
            continue
        fi = vertex_to_fiber[i]
        fj = vertex_to_fiber[j]
        if fi != fj:
            w[fi, fj] += 1.0
            w[fj, fi] += 1.0
    a = (w > 0).astype(float)
    return a, w


def oriented_fiber_triangles(a_fiber, fiber_xyz):
    triangles = []
    for i in range(12):
        for j in range(i + 1, 12):
            if a_fiber[i, j] < 0.5:
                continue
            for k in range(j + 1, 12):
                if a_fiber[i, k] > 0.5 and a_fiber[j, k] > 0.5:
                    det = np.linalg.det(
                        np.column_stack([fiber_xyz[i], fiber_xyz[j], fiber_xyz[k]])
                    )
                    if det >= 0:
                        triangles.append((i, j, k))
                    else:
                        triangles.append((i, k, j))
    return triangles


def build_triangle_bracket(oriented_triangles):
    c = np.zeros((12, 12, 12))
    for i, j, k in oriented_triangles:
        c[i, j, k] += 1.0
        c[j, k, i] += 1.0
        c[k, i, j] += 1.0
        c[j, i, k] -= 1.0
        c[k, j, i] -= 1.0
        c[i, k, j] -= 1.0
    return c


def bracket(c_tensor, x, y):
    return np.einsum("ijk,i,j->k", c_tensor, x, y)


def jacobi_max(c_tensor, basis):
    max_norm = 0.0
    worst = None
    for i in range(basis.shape[1]):
        for j in range(basis.shape[1]):
            for k in range(basis.shape[1]):
                xi = basis[:, i]
                xj = basis[:, j]
                xk = basis[:, k]
                term = (
                    bracket(c_tensor, xi, bracket(c_tensor, xj, xk))
                    + bracket(c_tensor, xj, bracket(c_tensor, xk, xi))
                    + bracket(c_tensor, xk, bracket(c_tensor, xi, xj))
                )
                nrm = norm(term)
                if nrm > max_norm:
                    max_norm = nrm
                    worst = (i, j, k)
    return max_norm, worst


def closure_report(c_tensor, basis_space, name):
    total = 0.0
    leak = 0.0
    for i in range(basis_space.shape[1]):
        for j in range(i + 1, basis_space.shape[1]):
            z = bracket(c_tensor, basis_space[:, i], basis_space[:, j])
            coeffs = basis_space.T @ z
            z_in = basis_space @ coeffs
            total += norm(z) ** 2
            leak += norm(z - z_in) ** 2
    frac = math.sqrt(leak / total) if total > 0 else 0.0
    print(f"  {name:<18} leakage = {frac:.6f}")


def commutator_report(c_tensor, left_basis, right_basis, target_basis, name):
    total = 0.0
    leak = 0.0
    for i in range(left_basis.shape[1]):
        for j in range(right_basis.shape[1]):
            z = bracket(c_tensor, left_basis[:, i], right_basis[:, j])
            coeffs = target_basis.T @ z
            z_in = target_basis @ coeffs
            total += norm(z) ** 2
            leak += norm(z - z_in) ** 2
    frac = math.sqrt(leak / total) if total > 0 else 0.0
    print(f"  {name:<18} leakage = {frac:.6f}")


def main():
    print("=" * 72)
    print("EXP608: GAUGE BRACKET CANDIDATES")
    print("=" * 72)

    print("\n[1] Building 600-cell complex and edge Box...")
    verts, adj, lap, edges, edge_to_idx, triangles, tri_edges, d0, d1 = build_complex()
    fibers = find_fibration(verts)
    box1, is_fiber_edge, vertex_to_fiber, fiber_edge_lists = build_edge_box(
        edges, triangles, d0, d1, fibers, edge_to_idx
    )
    all_eperm = build_group_edge_perms(verts, edges, edge_to_idx)
    gauge_basis, sv, ker = build_gauge_basis(box1, all_eperm)
    _ = adj, lap, ker
    print(f"  Kernel rho_5 projector eigenvalues: {np.round(sv, 8)}")
    support_outside = np.sum(gauge_basis[~is_fiber_edge] ** 2)
    print(f"  Gauge support outside fiber edges: {support_outside:.3e}")

    print("\n[2] Candidate A: simplicial wedge bracket on edge 1-forms")
    fiber_counts = []
    for i, j, k in triangles:
        count = int(vertex_to_fiber[i] == vertex_to_fiber[j])
        count += int(vertex_to_fiber[j] == vertex_to_fiber[k])
        count += int(vertex_to_fiber[i] == vertex_to_fiber[k])
        fiber_counts.append(count)
    unique_counts = sorted(set(fiber_counts))
    print(f"  Fiber-edge counts per triangle: {unique_counts}")
    print(f"  Max fiber edges in any triangle: {max(fiber_counts)}")

    wedge_max = 0.0
    bracket_max = 0.0
    for a in range(12):
        for b in range(a + 1, 12):
            w = wedge_to_faces(gauge_basis[:, a], gauge_basis[:, b], tri_edges)
            wedge_max = max(wedge_max, norm(w))
            back = d1.T @ w
            bracket_max = max(bracket_max, norm(back))
    print(f"  Max ||u ^ v|| on faces: {wedge_max:.3e}")
    print(f"  Max ||d1^T(u ^ v)|| on edges: {bracket_max:.3e}")
    if bracket_max < 1e-12:
        print("  Result: simplicial bracket vanishes exactly on the gauge sector.")

    print("\n[3] Candidate B: oriented triangle bracket on the 12 fibers")
    a_fiber, w_fiber = build_fiber_graph(edges, is_fiber_edge, vertex_to_fiber)
    degree = np.sum(a_fiber, axis=1)
    print(f"  Fiber graph degrees: {degree.astype(int)}")
    unique_weights = sorted({int(x) for x in w_fiber[w_fiber > 0].astype(int)})
    print(f"  Unique cross-edge weights: {unique_weights}")

    fiber_xyz = np.array([hopf_map(verts[fiber[0]]) for fiber in fibers])
    fiber_xyz /= np.linalg.norm(fiber_xyz, axis=1, keepdims=True)
    oriented_triangles = oriented_fiber_triangles(a_fiber, fiber_xyz)
    print(f"  Oriented fiber triangles: {len(oriented_triangles)} (expect 20 for icosahedron)")

    c_tensor = build_triangle_bracket(oriented_triangles)

    lf = np.diag(np.sum(a_fiber, axis=1)) - a_fiber
    evals_f, evecs_f = eigh(lf)
    rounded = np.round(evals_f, 10)
    print(f"  Fiber Laplacian eigenvalues: {rounded}")

    idx_u1 = np.where(np.abs(evals_f - 0.0) < 1e-8)[0]
    idx_su2 = np.where(np.abs(evals_f - (5.0 - math.sqrt(5.0))) < 1e-8)[0]
    idx_su3a = np.where(np.abs(evals_f - (5.0 + math.sqrt(5.0))) < 1e-8)[0]
    idx_su3b = np.where(np.abs(evals_f - 6.0) < 1e-8)[0]
    basis_u1 = evecs_f[:, idx_u1]
    basis_su2 = evecs_f[:, idx_su2]
    basis_su3 = evecs_f[:, np.concatenate([idx_su3a, idx_su3b])]
    full_basis = evecs_f

    central_norms = []
    for j in range(full_basis.shape[1]):
        central_norms.append(norm(bracket(c_tensor, basis_u1[:, 0], full_basis[:, j])))
    print(f"  U(1) centrality max ||[u1, x]||: {max(central_norms):.6f}")

    closure_report(c_tensor, basis_su2, "su(2) candidate")
    closure_report(c_tensor, basis_su3, "su(3) candidate")
    commutator_report(c_tensor, basis_su2, basis_su3, basis_su2, "[su2, su3] -> su2")
    commutator_report(c_tensor, basis_su2, basis_su3, basis_su3, "[su2, su3] -> su3")

    jacobi, worst = jacobi_max(c_tensor, full_basis)
    print(f"  Jacobi max norm on Laplacian basis: {jacobi:.6f} at {worst}")

    print("\n[4] Verdict")
    print("  Candidate A is exactly zero on fiber-localized gauge modes.")
    if max(central_norms) < 1e-8 and jacobi < 1e-8:
        print("  Candidate B passes centrality and Jacobi at numerical precision.")
    else:
        print("  Candidate B does not define the needed U(1)+SU(2)+SU(3) Lie algebra.")
    print("  If Candidate B fails, the 1+3+8 split remains a gauge-dimension skeleton,")
    print("  not a fully derived discrete Lie bracket.")

    print("\n" + "=" * 72)
    print("EXP608 COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
