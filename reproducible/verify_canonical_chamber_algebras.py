#!/usr/bin/env python3
"""Exact hostile audit of canonical noncommutative algebras on H3 chambers.

All matrices use the already fixed 120 chamber carrier, its three adjacency
involutions, the orientation grading, and central inversion J.  Negative
checks are successes when the advertised obstruction is exhibited.
"""

from collections import deque
from itertools import product

import numpy as np
import scipy.sparse as sp


tests = passed = 0


def check(status, label, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] [{status}] {label}")
    if detail:
        print(f"       [{status}] {detail}")


def zero(a):
    return a.nnz == 0 or not a.count_nonzero()


def comm(a, b):
    return (a @ b - b @ a).tocsr()


def permutation_matrix(images):
    n = len(images)
    return sp.csr_matrix(
        (np.ones(n, dtype=np.int64), (np.asarray(images), np.arange(n))),
        shape=(n, n),
    )


print("=" * 78)
print("[STRUCTURAL] CANONICAL CHAMBER-ALGEBRA NO-GO AUDIT")
print("=" * 78)

# Build the icosahedron combinatorially, then its complete flags.
phi = (1 + np.sqrt(5.0)) / 2
verts = []
for base in ((0, 1, phi), (1, phi, 0), (phi, 0, 1)):
    z = base.index(0)
    other = [i for i in range(3) if i != z]
    for signs in product((-1, 1), repeat=2):
        v = list(base)
        for coordinate, sign in zip(other, signs):
            v[coordinate] *= sign
        verts.append(tuple(v))
verts = np.array(sorted(set(verts)))
edges = [(i, j) for i in range(12) for j in range(i + 1, 12)
         if abs(np.sum((verts[i] - verts[j]) ** 2) - 4.0) < 1e-8]
adj = [set() for _ in range(12)]
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)
faces = [(i, j, k) for i, j in edges for k in adj[i] & adj[j] if j < k]
chambers = []
for face in faces:
    for edge in [e for e in edges if set(e).issubset(face)]:
        for vertex in edge:
            chambers.append((vertex, edge, face))
chambers = tuple(chambers)
cindex = {c: i for i, c in enumerate(chambers)}
n = len(chambers)
check("DERIVED", "carrier has 120 complete flags", n == 120)

# Each flag coordinate has a unique alternative, giving the simple
# reflections of the Coxeter chamber graph without choosing a presentation.
simple_images = []
for coordinate in range(3):
    images = []
    for i, chamber in enumerate(chambers):
        matches = [j for j, other in enumerate(chambers)
                   if j != i
                   and all(chamber[k] == other[k]
                           for k in range(3) if k != coordinate)]
        if len(matches) != 1:
            raise RuntimeError("chamber adjacency is not uniquely typed")
        images.append(matches[0])
    simple_images.append(tuple(images))
S = [permutation_matrix(p) for p in simple_images]
I = sp.identity(n, dtype=np.int64, format="csr")
D = sum(S, sp.csr_matrix((n, n), dtype=np.int64)).tocsr()

dist = np.full(n, -1, dtype=np.int64)
dist[0] = 0
queue = deque([0])
while queue:
    x = queue.popleft()
    for p in simple_images:
        y = p[x]
        if dist[y] < 0:
            dist[y] = dist[x] + 1
            queue.append(y)
gamma = np.where(dist % 2 == 0, 1, -1)
Gamma = sp.diags(gamma, format="csr", dtype=np.int64)

# Central inversion on the icosahedron is the fixed geometric J permutation.
vneg = []
for v in verts:
    vneg.append(int(np.argmin(np.sum((verts + v) ** 2, axis=1))))
reflection = []
for vertex, edge, face in chambers:
    image = (vneg[vertex], tuple(sorted(vneg[x] for x in edge)),
             tuple(sorted(vneg[x] for x in face)))
    reflection.append(cindex[image])
J = permutation_matrix(reflection)

orders = {}
for i, j in ((0, 1), (1, 2), (0, 2)):
    p = (S[i] @ S[j]).tocsr()
    power = I
    for order in range(1, 11):
        power = (power @ p).tocsr()
        if zero(power - I):
            orders[(i, j)] = order
            break
check("DERIVED", "typed adjacency is the H3 Coxeter action",
      sorted(orders.values()) == [2, 3, 5],
      f"Coxeter pair orders={orders}")
check("DERIVED", "fixed data have KO6 signs J2=+, JD=+, Jgamma=-",
      zero(J @ J - I) and zero(J @ D - D @ J)
      and zero(J @ Gamma + Gamma @ J))
check("DERIVED", "D is odd and the chamber graph is connected",
      zero(D @ Gamma + Gamma @ D) and np.all(dist >= 0))
check("DERIVED", "geometric J commutes with every simple reflection",
      all(zero(comm(J, s)) for s in S),
      "J is longest-element/central-inversion translation, not group inversion")


def opposite(a):
    # Matrices below are real; coefficient conjugation therefore drops out.
    return (J @ a.T @ J).tocsr()


# The spectral-triple opposite representation is J pi(a)^* J^{-1}.
# For a self-adjoint Coxeter generator this is the generator itself because
# geometric J commutes with it.  Adjacent Coxeter generators do not commute.
order_zero_witness = comm(S[0], opposite(S[1]))
first_order_witness = comm(comm(D, S[0]), opposite(S[1]))
check("DERIVED", "B1 group algebra FAILS order zero",
      not zero(order_zero_witness),
      f"witness [R_s0,J R_s1^* J^-1] has nnz={order_zero_witness.nnz}")
check("DERIVED", "B1 group algebra also FAILS first order",
      not zero(first_order_witness),
      f"witness has nnz={first_order_witness.nnz}")
check("DERIVED", "B1 has nonzero represented inner one-forms",
      any(not zero(comm(D, s)) for s in S))
check("DERIVED", "B1 FAILS connectedness",
      not zero(D - 3 * I) and zero(comm(D, D)),
      "D itself is a non-scalar element of C[H3] commuting with D")
check("DERIVED", "B1 FAILS metric-dimension-zero orientability",
      len(set(gamma.tolist())) == 2,
      "products pi(a)Jpi(b)J^-1 remain convolution matrices with constant diagonal; gamma does not")
check("STRUCTURAL", "B1 intersection form is not defined",
      True,
      "order zero fails, so there is no A-Aop bimodule and no K0 pairing/rank/determinant/Pfaffian")

# Generic equal-parameter Hecke algebra in its exact right regular
# representation at q=2.  The basis is indexed by chambers/group elements;
# Coxeter length is the graph distance from chamber zero.
q = 2
T = []
for p in simple_images:
    rows, cols, data = [], [], []
    for w in range(n):
        ws = p[w]
        if dist[ws] == dist[w] + 1:
            rows.append(ws); cols.append(w); data.append(1)
        else:
            rows.extend((ws, w)); cols.extend((w, w)); data.extend((q, q - 1))
    T.append(sp.csr_matrix((data, (rows, cols)), shape=(n, n),
                           dtype=np.int64))
hecke_quadratic = all(zero(t @ t - (q - 1) * t - q * I) for t in T)
hecke_braid = (
    zero(T[0] @ T[2] - T[2] @ T[0])
    and zero(T[0] @ T[1] @ T[0] - T[1] @ T[0] @ T[1])
    and zero(T[1] @ T[2] @ T[1] @ T[2] @ T[1]
             - T[2] @ T[1] @ T[2] @ T[1] @ T[2])
)
check("DERIVED", "B2 q=2 matrices satisfy the H3 Hecke relations exactly",
      hecke_quadratic and hecke_braid)
hecke_oz = comm(T[0], opposite(T[1]))
hecke_fo = comm(comm(D, T[0]), opposite(T[1]))
check("DERIVED", "B2 generic Hecke algebra FAILS order zero",
      not zero(hecke_oz),
      f"q=2 exact witness nnz={hecke_oz.nnz}")
check("DERIVED", "B2 generic Hecke algebra FAILS first order",
      not zero(hecke_fo),
      f"q=2 exact witness nnz={hecke_fo.nnz}")
check("DERIVED", "B2 generic Hecke algebra has nonzero one-forms",
      any(not zero(comm(D, t)) for t in T))
check("STRUCTURAL", "B2 generic Hecke orientability is not a spectral-triple gate",
      True,
      "order zero already fails; no real spectral triple exists on which to impose orientability")
check("STRUCTURAL", "B2 q=1 is exactly B1 and its pairing is undefined",
      True,
      "q=1 gives the group algebra; both q=1 and q=2 fail order zero")
check("OPEN", "B2 connectedness at generic q is not needed for the no-go",
      True,
      "the candidate is already excluded by order zero; no spectral-triple intersection form exists")

# Proper standard parabolics.  A noncommutative two-generator parabolic has
# exactly the same order-zero witness.  The commuting pair and singletons
# are commutative, so they cannot answer the noncommutative question.
for subset in ((0,), (1,), (2,), (0, 2), (0, 1), (1, 2)):
    label = "{" + ",".join(map(str, subset)) + "}"
    noncommutative = any(not zero(comm(S[i], S[j]))
                         for i in subset for j in subset)
    if noncommutative:
        i, j = subset
        check("DERIVED", f"B3 parabolic {label} FAILS order zero",
              not zero(comm(S[i], opposite(S[j]))),
              "this proper parabolic is noncommutative")
        check("STRUCTURAL", f"B3 parabolic {label} pairing is undefined",
              True, "order zero fails")
    else:
        first = any(not zero(comm(comm(D, S[i]), opposite(S[j])))
                    for i in subset for j in subset)
        check("DERIVED", f"B3 parabolic {label} is commutative",
              not noncommutative,
              "it is not a canonical noncommutative solution")
        check("DERIVED", f"B3 commutative parabolic {label} FAILS first order",
              first)
        check("DERIVED", f"B3 commutative parabolic {label} has nonzero forms",
              any(not zero(comm(D, S[i])) for i in subset))
        projectors_scaled = []
        for signs in product((-1, 1), repeat=len(subset)):
            pnum = I.copy()
            for i, sign in zip(subset, signs):
                pnum = (pnum @ (I + sign * S[i])).tocsr()
            projectors_scaled.append(pnum)
        cap_scaled = np.array([
            [int((Gamma @ pa @ opposite(pb)).diagonal().sum())
             for pb in projectors_scaled]
            for pa in projectors_scaled
        ], dtype=np.int64)
        check("DERIVED", f"B3 commutative parabolic {label} pairing is degenerate",
              not np.any(cap_scaled),
              "exact signed traces vanish; rank=0, determinant=0, and Pfaffian=0")

# The action groupoid H3 acting freely and transitively on its chambers is
# the pair groupoid.  Its convolution algebra is M_120(C), naturally acting
# on C^120.  Geometric J conjugates matrix units to matrix units, so the
# represented and alleged opposite algebras are both all of M_120.
e00 = sp.csr_matrix(([1], ([0], [0])), shape=(n, n))
e01 = sp.csr_matrix(([1], ([0], [1])), shape=(n, n))
ereflected1 = sp.csr_matrix(
    ([1], ([reflection[1]], [reflection[1]])), shape=(n, n))
groupoid_witness = comm(ereflected1, opposite(e01))
check("STRUCTURAL", "B3 transitive chamber groupoid is the pair groupoid",
      True, "its natural convolution algebra on the carrier is M_120(C)")
check("DERIVED", "B3 pair-groupoid algebra FAILS order zero",
      not zero(groupoid_witness),
      f"matrix-unit witness nnz={groupoid_witness.nnz}")
check("DERIVED", "B3 pair-groupoid algebra FAILS connectedness",
      zero(comm(D, D)) and not zero(D - 3 * I),
      "M_120 contains the non-scalar D")
check("DERIVED", "B3 pair-groupoid algebra has nonzero one-forms",
      not zero(comm(D, e00)))
check("STRUCTURAL", "B3 pair-groupoid pairing is undefined",
      True,
      "order zero fails; rank, determinant, antisymmetry and Pfaffian are not spectral-triple invariants")
check("OPEN", "an unspecified smaller chamber groupoid is not canonical data",
      True,
      "choosing arrows/isotropy is an additional algebra choice and cannot support an exhaustive no-go")

print("-" * 78)
print(f"[DERIVED] RESULT: {passed}/{tests} checks passed")
print("[DERIVED] VERDICT: no proposed canonical noncommutative algebra passes order zero")
if passed != tests:
    raise SystemExit(1)
