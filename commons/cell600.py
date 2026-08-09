"""
600-cell and E8 root constructions.
Verified: 120 vertices, degree 12, 720 edges, 9 distinct eigenvalues.
Source: reproducible/verify_spectrum_600cell.py (ALL PASS).
"""
import numpy as np
from itertools import permutations, product as cartesian_product

PHI = (1 + np.sqrt(5)) / 2
EDGE_TOL = 1e-6


def build_600cell():
    """
    Build the 600-cell as the Cayley graph of the binary icosahedral group 2I.

    Returns:
        verts: (120, 4) array of unit quaternions on S^3
        adj:   (120, 120) adjacency matrix (degree 12)
        lap:   (120, 120) graph Laplacian (12*I - adj)
    """
    verts = _build_2I_vertices()
    adj = _build_adjacency(verts)
    lap = 12 * np.eye(120) - adj
    return verts, adj, lap


def _build_2I_vertices():
    """120 unit quaternions of binary icosahedral group 2I."""
    verts = set()

    def add_vert(v):
        arr = np.array(v, dtype=float)
        n = np.linalg.norm(arr)
        if n > 1e-12:
            arr = arr / n
        verts.add(tuple(np.round(arr, 10)))

    # Type A: 8 axis quaternions
    for i in range(4):
        for s in [1.0, -1.0]:
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add_vert(v)

    # Type B: 16 half-integer quaternions
    for s0, s1, s2, s3 in cartesian_product([0.5, -0.5], repeat=4):
        add_vert([s0, s1, s2, s3])

    # Type C: 96 golden-ratio quaternions
    # Even permutations of (0, +-1/2, +-phi/2, +-1/(2*phi))
    base = [0.0, 0.5, PHI / 2.0, 1.0 / (2.0 * PHI)]

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

    result = np.array(sorted(verts))
    assert len(result) == 120, f"Expected 120 vertices, got {len(result)}"
    return result


def _build_adjacency(verts):
    """Build adjacency matrix. Edge criterion: dot(q1,q2) = phi/2."""
    n = len(verts)
    dot_threshold = PHI / 2.0
    dots = verts @ verts.T
    np.clip(dots, -1.0, 1.0, out=dots)

    adj = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(dots[i, j] - dot_threshold) < EDGE_TOL:
                adj[i, j] = 1.0
                adj[j, i] = 1.0

    degrees = np.sum(adj, axis=1).astype(int)
    assert np.all(degrees == 12), f"Degree range: [{min(degrees)}, {max(degrees)}]"
    assert int(np.sum(adj) / 2) == 720, "Expected 720 edges"
    return adj


def build_e8_roots():
    """
    Build the 240 roots of E8 in R^8.

    Returns:
        roots: (240, 8) array
        simple: (8, 8) array of simple roots
    """
    roots = []
    # 112 roots: +-e_i +- e_j
    for i in range(8):
        for j in range(i + 1, 8):
            for si in [1, -1]:
                for sj in [1, -1]:
                    v = np.zeros(8)
                    v[i] = si
                    v[j] = sj
                    roots.append(v)

    # 128 roots: (+-1/2)^8 with even number of minus signs
    for signs in cartesian_product([1, -1], repeat=8):
        if sum(1 for s in signs if s == -1) % 2 == 0:
            roots.append(np.array(signs) / 2.0)

    roots = np.array(roots)
    assert len(roots) == 240, f"Expected 240 roots, got {len(roots)}"

    # Simple roots (Bourbaki convention)
    simple = np.zeros((8, 8))
    for i in range(6):
        simple[i, i] = 1
        simple[i, i + 1] = -1
    simple[6] = np.array([0, 0, 0, 0, 0, 1, 1, 0])
    simple[7] = np.array([-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5])

    return roots, simple
