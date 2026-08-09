"""
verify_discrete_scalar_response.py
==================================
Self-contained verification of the exact static scalar-response theorem on the
600-cell graph.

What this script verifies:
  1. The 600-cell graph and boundary maps d0, d1 are built directly.
  2. The Moore-Penrose identity
       B^+ d0 = d0 Delta_0^+
     holds numerically on the full complex, where
       Delta_0 = d0^T d0
       B = d0 d0^T.
  3. For every point source delta_v on vertices, the induced edge response
       h_v = B^+ d0 delta_v
     is exactly a gradient:
       h_v = d0 Phi_v,   Phi_v = Delta_0^+ delta_v.
  4. The response has no coexact content:
       d1 h_v = 0
     for every source vertex.
  5. The scalar potential is recovered from h_v without loss:
       d0^+ h_v = Phi_v.

This is the exact discrete theorem behind the paper's gamma_disc = 1 language.
It does NOT claim a full continuum PPN derivation.
"""

import time
from collections import defaultdict
from itertools import permutations, product

import numpy as np


a1 = 5
b1 = a1 + 1
phi = (1 + np.sqrt(a1)) / 2
N = 120

results = []


def record(name, passed, detail=""):
    results.append((name, passed, detail))
    tag = "PASS" if passed else "FAIL"
    print(f"  [{tag}] {name}")
    if detail:
        print(f"         {detail}")


print("=" * 72)
print("VERIFY DISCRETE SCALAR RESPONSE ON THE 600-CELL")
print("=" * 72)
print()
print(f"Constants: a1={a1}, b1={b1}, phi={phi:.10f}, N={N}")

# ============================================================================
# SECTION 1: BUILD THE 600-CELL GRAPH
# ============================================================================
print()
print("-" * 72)
print("SECTION 1: Build the 600-cell graph")
print("-" * 72)
t0 = time.time()

verts_set = set()

for i in range(4):
    for s in [1.0, -1.0]:
        v = [0.0, 0.0, 0.0, 0.0]
        v[i] = s
        verts_set.add(tuple(v))

for s0 in [0.5, -0.5]:
    for s1 in [0.5, -0.5]:
        for s2 in [0.5, -0.5]:
            for s3 in [0.5, -0.5]:
                verts_set.add((s0, s1, s2, s3))

base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
even_perms = []
for p in permutations(range(4)):
    inv_count = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
    if inv_count % 2 == 0:
        even_perms.append(p)

for perm in even_perms:
    coords = [base[perm[i]] for i in range(4)]
    nonzero_idx = [i for i in range(4) if abs(coords[i]) > 1e-12]
    for signs in product([1, -1], repeat=len(nonzero_idx)):
        v = list(coords)
        for idx, s in zip(nonzero_idx, signs):
            v[idx] *= s
        verts_set.add(tuple(round(x, 10) for x in v))

verts = np.array(sorted(verts_set))
Nv = len(verts)

dots = verts @ verts.T
edge_thresh = phi / 2
edges = []
for i in range(Nv):
    for j in range(i + 1, Nv):
        if abs(dots[i, j] - edge_thresh) < 0.001:
            edges.append((i, j))
Ne = len(edges)

edge_to_idx = {}
for idx, (i, j) in enumerate(edges):
    edge_to_idx[(i, j)] = idx
    edge_to_idx[(j, i)] = idx

adj_list = defaultdict(set)
for i, j in edges:
    adj_list[i].add(j)
    adj_list[j].add(i)

triangles = []
for i in range(Nv):
    for j in adj_list[i]:
        if j > i:
            common = adj_list[i] & adj_list[j]
            for k in common:
                if k > j:
                    triangles.append((i, j, k))
Nf = len(triangles)

build_time = time.time() - t0
print(f"  Built graph in {build_time:.1f}s")
print(f"  Vertices={Nv}, Edges={Ne}, Faces={Nf}")

record("Vertices = 120", Nv == 120, f"got {Nv}")
record("Edges = 720", Ne == 720, f"got {Ne}")
record("Faces = 1200", Nf == 1200, f"got {Nf}")

# ============================================================================
# SECTION 2: d0, d1, Laplacians, pseudoinverses
# ============================================================================
print()
print("-" * 72)
print("SECTION 2: Boundary maps and pseudoinverse identity")
print("-" * 72)
t0 = time.time()

d0 = np.zeros((Ne, Nv))
for e_idx, (i, j) in enumerate(edges):
    d0[e_idx, i] = -1.0
    d0[e_idx, j] = 1.0

d1 = np.zeros((Nf, Ne))
for f_idx, (i, j, k) in enumerate(triangles):
    d1[f_idx, edge_to_idx[(i, j)]] = 1.0
    d1[f_idx, edge_to_idx[(j, k)]] = 1.0
    d1[f_idx, edge_to_idx[(i, k)]] = -1.0

Delta0 = d0.T @ d0
B = d0 @ d0.T

Delta0_pinv = np.linalg.pinv(Delta0, rcond=1e-12)
B_pinv = np.linalg.pinv(B, rcond=1e-12)
d0_pinv = np.linalg.pinv(d0, rcond=1e-12)

lhs = B_pinv @ d0
rhs = d0 @ Delta0_pinv
identity_residual = np.max(np.abs(lhs - rhs))

compute_time = time.time() - t0
print(f"  Built operators and pseudoinverses in {compute_time:.1f}s")
print(f"  max |B^+ d0 - d0 Delta0^+| = {identity_residual:.2e}")

record("Pseudoinverse identity B^+ d0 = d0 Delta0^+",
       identity_residual < 1e-10,
       f"max residual = {identity_residual:.2e}")

# ============================================================================
# SECTION 3: Static point sources
# ============================================================================
print()
print("-" * 72)
print("SECTION 3: Static response to point masses")
print("-" * 72)
t0 = time.time()

P0 = np.eye(Nv) - np.ones((Nv, Nv)) / Nv

max_poisson_residual = 0.0
max_gradient_residual = 0.0
max_coexact_residual = 0.0
max_recovery_residual = 0.0
max_mean_phi = 0.0

for v in range(Nv):
    delta_v = np.zeros(Nv)
    delta_v[v] = 1.0

    Phi = Delta0_pinv @ delta_v
    h = B_pinv @ (d0 @ delta_v)
    h_grad = d0 @ Phi
    Phi_rec = d0_pinv @ h

    poisson_residual = np.max(np.abs(Delta0 @ Phi - P0 @ delta_v))
    gradient_residual = np.max(np.abs(h - h_grad))
    coexact_residual = np.max(np.abs(d1 @ h))
    recovery_residual = np.max(np.abs(Phi_rec - Phi))
    mean_phi = abs(np.mean(Phi))

    max_poisson_residual = max(max_poisson_residual, poisson_residual)
    max_gradient_residual = max(max_gradient_residual, gradient_residual)
    max_coexact_residual = max(max_coexact_residual, coexact_residual)
    max_recovery_residual = max(max_recovery_residual, recovery_residual)
    max_mean_phi = max(max_mean_phi, mean_phi)

elapsed = time.time() - t0
print(f"  Tested all {Nv} point sources in {elapsed:.1f}s")
print(f"  max |Delta0 Phi - P0 delta_v| = {max_poisson_residual:.2e}")
print(f"  max |h - d0 Phi|             = {max_gradient_residual:.2e}")
print(f"  max |d1 h|                   = {max_coexact_residual:.2e}")
print(f"  max |d0^+ h - Phi|           = {max_recovery_residual:.2e}")
print(f"  max |mean(Phi)|              = {max_mean_phi:.2e}")

record("Poisson equation holds on the zero-mean subspace",
       max_poisson_residual < 1e-10,
       f"max residual = {max_poisson_residual:.2e}")
record("Static edge response is exactly a gradient",
       max_gradient_residual < 1e-10,
       f"max residual = {max_gradient_residual:.2e}")
record("Static response has no coexact component",
       max_coexact_residual < 1e-10,
       f"max |d1 h| = {max_coexact_residual:.2e}")
record("Potential is recovered uniquely from the edge field",
       max_recovery_residual < 1e-10,
       f"max residual = {max_recovery_residual:.2e}")
record("Pseudoinverse potentials are zero-mean",
       max_mean_phi < 1e-12,
       f"max |mean(Phi)| = {max_mean_phi:.2e}")

# ============================================================================
# SUMMARY
# ============================================================================
print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)

n_pass = sum(1 for _, ok, _ in results if ok)
n_fail = len(results) - n_pass

for name, ok, detail in results:
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}")

print()
print(f"  Total: {n_pass}/{len(results)} passed, {n_fail} failed")
print()

if n_fail == 0:
    print("  ALL TESTS PASSED.")
    print("  Exact discrete scalar-response theorem verified:")
    print("    - B^+ d0 = d0 Delta0^+")
    print("    - each point source induces h = d0 Phi")
    print("    - the response is purely exact, with d1 h = 0")
    print("    - d0^+ h recovers the same scalar potential")
    print()
    print("  Interpretation:")
    print("    Static vertex sources do not excite an independent coexact")
    print("    scalar mode on the tested 600-cell sector.")
else:
    print("  WARNING: at least one scalar-response check failed.")

print("=" * 72)
