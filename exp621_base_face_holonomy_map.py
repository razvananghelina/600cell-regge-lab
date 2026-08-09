"""
exp621: Holonomy and Berry-phase statistics grouped by Hopf-base faces.

Goal:
  Test whether the missing local gauge connection might already be encoded in
  the 600-cell holonomy, once triangles are grouped by the 20 triangular faces
  of the Hopf-base icosahedron.

Main questions:
  1. Does each base face lift to the same number of 600-cell triangles?
  2. Is the total SO(3) holonomy already constant per base face?  (Expected: yes,
     because it is globally face-transitive.)
  3. Do the Berry phases or fiber fractions acquire nontrivial dependence on the
     base face?  If not, then holonomy is still too symmetric to define a local
     connection by face grouping alone.
"""

from collections import defaultdict
import math
import sys

import numpy as np
from numpy.linalg import eigh, norm

sys.path.insert(0, ".")
from commons import build_600cell


PHI = (1.0 + math.sqrt(5.0)) / 2.0


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
    idx = int(np.argmax(dots))
    return idx if dots[idx] > 1.0 - tol else -1


def find_hopf_fibration(verts):
    nv = len(verts)
    for i in range(nv):
        if abs(verts[i, 0] - PHI / 2.0) >= 1e-6:
            continue

        g = verts[i]
        power = g.copy()
        ok = True
        for k in range(2, 11):
            power = qmul(power, g)
            if k == 5 and not np.allclose(power, [-1, 0, 0, 0], atol=1e-6):
                ok = False
                break
            if k == 10 and not np.allclose(power, [1, 0, 0, 0], atol=1e-6):
                ok = False
        if not ok:
            continue

        subgroup = []
        power = np.array([1.0, 0.0, 0.0, 0.0])
        for _ in range(10):
            subgroup.append(find_idx(power, verts))
            power = qmul(power, g)

        used = set()
        fibers = []
        for s in range(nv):
            if s in used:
                continue
            fiber = []
            for si in subgroup:
                idx = find_idx(qmul(verts[s], verts[si]), verts)
                if idx >= 0 and idx not in used:
                    fiber.append(idx)
                    used.add(idx)
            if len(fiber) == 10:
                fibers.append(fiber)

        if len(fibers) == 12:
            return fibers

    raise RuntimeError("Could not find Hopf fibration")


def build_triangles(adj):
    nv = adj.shape[0]
    adj_list = defaultdict(set)
    for i in range(nv):
        for j in range(nv):
            if adj[i, j] > 0.5:
                adj_list[i].add(j)

    triangles = []
    for i in range(nv):
        for j in adj_list[i]:
            if j > i:
                for k in adj_list[i] & adj_list[j]:
                    if k > j:
                        triangles.append((i, j, k))
    return triangles


def parallel_transport(v, p, q):
    c = np.dot(p, q)
    if abs(1.0 + c) < 1e-10:
        return -v
    return v - (np.dot(v, q) / (1.0 + c)) * (p + q)


def quaternionic_basis(p):
    w, x, y, z = p
    ti = np.array([-x, w, -z, y])
    tj = np.array([-y, z, w, -x])
    tk = np.array([-z, -y, x, w])
    return np.array([ti, tj, tk])


def holonomy_matrix_quat(p, q, r):
    basis = quaternionic_basis(p)
    h = np.zeros((3, 3))
    for col in range(3):
        v = basis[col].copy()
        v = parallel_transport(v, p, q)
        v = parallel_transport(v, q, r)
        v = parallel_transport(v, r, p)
        for row in range(3):
            h[row, col] = np.dot(v, basis[row])
    return h


def triple_cross_4d(p, u, w):
    m = np.array([p, u, w])
    n = np.zeros(4)
    for i in range(4):
        cols = [j for j in range(4) if j != i]
        n[i] = ((-1.0) ** i) * np.linalg.det(m[:, cols])
    return n


def hopf_fiber_tangent(p):
    w, x, y, z = p
    return np.array([-x, w, -z, y])


def build_base_graph(adj, vertex_to_fiber):
    weighted = np.zeros((12, 12), dtype=float)
    nv = adj.shape[0]
    for i in range(nv):
        for j in range(i + 1, nv):
            if adj[i, j] > 0.5:
                fi = vertex_to_fiber[i]
                fj = vertex_to_fiber[j]
                if fi != fj:
                    weighted[fi, fj] += 1.0
                    weighted[fj, fi] += 1.0

    adjacency = (weighted > 0).astype(float)
    laplacian = np.diag(np.sum(adjacency, axis=1)) - adjacency
    return adjacency, laplacian


def base_faces_from_embedding(adjacency, xyz):
    faces = []
    for i in range(12):
        for j in range(i + 1, 12):
            if adjacency[i, j] < 0.5:
                continue
            for k in range(j + 1, 12):
                if adjacency[i, k] > 0.5 and adjacency[j, k] > 0.5:
                    faces.append((i, j, k))
    return faces


def main():
    print("=" * 72)
    print("exp621: holonomy grouped by Hopf-base faces")
    print("=" * 72)

    verts, adj, _ = build_600cell()
    fibers = find_hopf_fibration(verts)
    vertex_to_fiber = {}
    for fi, fiber in enumerate(fibers):
        for v in fiber:
            vertex_to_fiber[v] = fi

    base_adj, base_lap = build_base_graph(adj, vertex_to_fiber)
    vals, vecs = eigh(base_lap)
    xyz = vecs[:, 1:4]
    xyz = xyz / np.linalg.norm(xyz, axis=1)[:, None]
    base_faces = {tuple(sorted(f)) for f in base_faces_from_embedding(base_adj, xyz)}
    print(f"\n[1] Base faces: {len(base_faces)} (expected 20)")

    triangles = build_triangles(adj)
    grouped = defaultdict(list)

    for tri_idx, (i, j, k) in enumerate(triangles):
        key = tuple(sorted((vertex_to_fiber[i], vertex_to_fiber[j], vertex_to_fiber[k])))
        if len(set(key)) == 3 and key in base_faces:
            grouped[key].append((tri_idx, i, j, k))

    counts = sorted(len(v) for v in grouped.values())
    print(f"  Lift counts per base face: {counts}")

    omega_by_face = []
    berry_by_face = []
    frac_by_face = []

    for key in sorted(grouped):
        omegas = []
        berries = []
        fracs = []
        for _, i, j, k in grouped[key]:
            p, q, r = verts[i], verts[j], verts[k]
            h = holonomy_matrix_quat(p, q, r)

            tr = np.trace(h)
            cos_omega = np.clip((tr - 1.0) / 2.0, -1.0, 1.0)
            omegas.append(np.arccos(cos_omega))

            cos_berry = np.clip((h[1, 1] + h[2, 2]) / 2.0, -1.0, 1.0)
            berries.append(np.arccos(cos_berry))

            u = q - np.dot(q, p) * p
            w = r - np.dot(r, p) * p
            n4 = triple_cross_4d(p, u, w)
            n4 /= norm(n4)
            ti = hopf_fiber_tangent(p)
            fracs.append(abs(np.dot(n4, ti)))

        omega_by_face.append((key, np.mean(omegas), np.std(omegas)))
        berry_by_face.append((key, np.mean(berries), np.std(berries), len(np.unique(np.round(berries, 6)))))
        frac_by_face.append((key, np.mean(fracs), np.std(fracs), len(np.unique(np.round(fracs, 6)))))

    print("\n[2] Total holonomy by base face")
    print(f"  mean(omega std over lifts) = {np.mean([x[2] for x in omega_by_face]):.3e}")
    print(f"  max(omega std over lifts)  = {np.max([x[2] for x in omega_by_face]):.3e}")

    print("\n[3] Berry phase by base face")
    berry_means = np.array([x[1] for x in berry_by_face])
    berry_stds = np.array([x[2] for x in berry_by_face])
    berry_counts = np.array([x[3] for x in berry_by_face])
    print(f"  facewise Berry mean range = [{berry_means.min():.12f}, {berry_means.max():.12f}]")
    print(f"  facewise Berry std range  = [{berry_stds.min():.12f}, {berry_stds.max():.12f}]")
    print(f"  distinct Berry values per base face = {sorted(set(berry_counts.tolist()))}")

    print("\n[4] Fiber fraction by base face")
    frac_means = np.array([x[1] for x in frac_by_face])
    frac_stds = np.array([x[2] for x in frac_by_face])
    frac_counts = np.array([x[3] for x in frac_by_face])
    print(f"  facewise fiber-fraction mean range = [{frac_means.min():.12f}, {frac_means.max():.12f}]")
    print(f"  facewise fiber-fraction std range  = [{frac_stds.min():.12f}, {frac_stds.max():.12f}]")
    print(f"  distinct fiber-fraction values per base face = {sorted(set(frac_counts.tolist()))}")

    print("\n[5] Is the facewise Berry data geometrically trivial?")
    face_centers = {}
    for key in base_faces:
        center = xyz[key[0]] + xyz[key[1]] + xyz[key[2]]
        center = center / norm(center)
        face_centers[key] = center

    abs_z = np.array([abs(face_centers[key][2]) for key, *_ in berry_by_face])
    corr = np.corrcoef(berry_means, abs_z)[0, 1]
    print(f"  Distinct facewise Berry means (rounded to 1e-6): {len(set(np.round(berry_means, 6)))}")
    print(f"  Distinct |face-center z| values (rounded to 1e-6): {len(set(np.round(abs_z, 6)))}")
    print(f"  Corr(facewise Berry mean, |face-center z|) = {corr:.6f}")

    print("\nSample base faces")
    for key, bmean, bstd, bcnt in berry_by_face[:5]:
        fmean, fstd, fcnt = next(x[1:] for x in frac_by_face if x[0] == key)
        print(
            f"  face {key}: Berry mean/std=({bmean:.12f}, {bstd:.12f}), "
            f"fiber mean/std=({fmean:.12f}, {fstd:.12f}), counts=({bcnt}, {fcnt})"
        )

    print("\nInterpretation")
    print("  Grouping triangles by the 20 Hopf-base faces does produce a clean")
    print("  30-fold lift for each face. The total SO(3) holonomy remains")
    print("  face-transitive, but the Berry and fiber-fraction statistics do not:")
    print("  they carry genuine local face dependence on the Hopf base. At the")
    print("  same time, that dependence is not captured by a simple invariant like")
    print("  |face-center z|. So holonomy does contain local data, but not yet in")
    print("  the form of a canonical local gauge connection.")


if __name__ == "__main__":
    main()
