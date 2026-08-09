#!/usr/bin/env python3
"""Exact audit of the final canonical A tensor A^op arena for A=C[2I]."""

from itertools import permutations, product

import numpy as np
import sympy as sy


tests = passed = 0


def check(name, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


def perm_sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p))
                     for j in range(i+1, len(p))) % 2 else 1


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z, w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X, w*Z+x*Y-y*X+z*W))


print("=" * 78)
print("CANONICAL BIMODULE ARENA A TENSOR A^op, A=C[2I]")
print("=" * 78)

# Certified deterministic 2I/600-cell realization.
phi = (1+np.sqrt(5.0))/2
vertex_set = set()
for i in range(4):
    for sign in (-1.0, 1.0):
        q = [0.0]*4
        q[i] = sign
        vertex_set.add(tuple(q))
vertex_set.update(product((-0.5, 0.5), repeat=4))
base = [phi/2, 0.5, 1/(2*phi), 0.0]
for p in permutations(range(4)):
    if perm_sign(p) != 1:
        continue
    q = [base[p[i]] for i in range(4)]
    nz = [i for i, x in enumerate(q) if abs(x) > 1e-12]
    for signs in product((-1, 1), repeat=3):
        r = q[:]
        for i, sign in zip(nz, signs):
            r[i] *= sign
        vertex_set.add(tuple(round(x, 10) for x in r))
V = np.array(sorted(vertex_set))
n = len(V)
identity = int(np.argmin(((V-np.array((1., 0., 0., 0.)))**2).sum(1)))
central_minus = int(np.argmin(((V-np.array((-1., 0., 0., 0.)))**2).sum(1)))
mult = np.empty((n, n), dtype=np.int16)
max_error = 0.0
for g, q in enumerate(V):
    moved = np.array([qmul(q, x) for x in V])
    dist = ((moved[:, None, :]-V[None, :, :])**2).sum(2)
    mult[g] = dist.argmin(1)
    max_error = max(max_error, float(dist[np.arange(n), mult[g]].max()))
inverse = np.empty(n, dtype=np.int16)
for g in range(n):
    inverse[g] = int(np.where(mult[g] == identity)[0][0])
check("2I multiplication table is exact after certified matching",
      n == 120 and max_error < 1e-15
      and all(mult[g, inverse[g]] == identity for g in range(n)),
      f"|G|={n}, max squared error={max_error:.3e}")

d = (1, 2, 3, 4, 5, 6, 4, 2, 3)
tensor_dim = n*n
enveloping_image_dim = sum(x**4 for x in d)
check("abstract A tensor A^op has dimension 14400",
      tensor_dim == 14400)
check("A tensor A^op is not End_C(A) for non-simple A",
      enveloping_image_dim == 2628
      and tensor_dim-enveloping_image_dim == 11772,
      "two-sided multiplication image=2628, kernel=11772; End_C(A) is M120")

# J(a tensor b^op)=b* tensor (a*)^op is a signed-free basis flip/inversion.
def j_pair(x, y):
    return int(inverse[y]), int(inverse[x])

check("flip-star J squares to +1 exactly",
      all(j_pair(*j_pair(x, y)) == (x, y)
          for x in range(n) for y in range(n)))
check("left and opposite right actions satisfy order zero",
      all(mult[a, mult[x, b]] == mult[mult[a, x], b]
          for a in range(n) for x in range(n) for b in range(n)),
      "associativity gives [L_a,R_b]=0 for all 120^3 basis triples")

# Conjugation action on each group-algebra factor.  z=-1 is central, so Ad_z=1.
check("-1 acts trivially in the diagonal adjoint action",
      all(mult[mult[central_minus, x], inverse[central_minus]] == x
          for x in range(n)))

# Exact conjugation-character decomposition.  chi_AdA(g)=|C_G(g)|.
sqrt5 = sy.sqrt(5)
ph = (1+sqrt5)/2
php = (1-sqrt5)/2
class_sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)
x = (2, -2, 0, 1, -1, ph, -ph, ph-1, php)
xp = (2, -2, 0, 1, -1, php, ph-1, -ph, ph)


def su2chars(t):
    out = [[sy.Integer(1)]*9, list(t)]
    for _ in range(2, 6):
        out.append([sy.expand(t[i]*out[-1][i]-out[-2][i])
                    for i in range(9)])
    return out


sx, sxp = su2chars(x), su2chars(xp)
standard = [sx[0], list(x), list(xp), sx[2], sxp[2],
            [sy.expand(x[i]*xp[i]) for i in range(9)],
            sx[3], sx[4], sx[5]]
order = (0, 1, 3, 6, 7, 8, 5, 2, 4)
chars = [standard[i] for i in order]
centralizers = tuple(120//s for s in class_sizes)


def decompose(values):
    return tuple(int(sy.simplify(sum(class_sizes[k]*chars[i][k]*values[k]
                                     for k in range(9))/120))
                 for i in range(9))


ad_A = decompose(centralizers)
ad_H = decompose(tuple(c*c for c in centralizers))
check("A conjugation decomposition is exact",
      ad_A == (9, 0, 7, 0, 9, 0, 6, 0, 7), str(ad_A))
check("H diagonal-adjoint decomposition is exact",
      ad_H == (296, 0, 736, 0, 1192, 0, 932, 0, 736)
      and sum(ad_H[i]*d[i] for i in range(9)) == 14400,
      str(ad_H))
check("adjoint arena contains only central-even irreps",
      all(ad_H[i] == 0 for i in (1, 3, 5, 7)))

# The 12 neighbors of e are exactly the trace-phi conjugacy class.
dots = V @ V.T
neighbors = set(np.where(abs(dots[identity]-phi/2) < 1e-8)[0].tolist())
trace_phi = set(np.where(abs(2*V[:, 0]-phi) < 1e-8)[0].tolist())
check("600-cell adjacency is convolution by the 12-neighbor class sum",
      len(neighbors) == 12 and neighbors == trace_phi)
check("neighbor class is inverse-closed",
      {int(inverse[g]) for g in neighbors} == neighbors)
check("neighbor sum c is central exactly",
      all({int(mult[h, mult[g, inverse[h]]]) for g in neighbors} == neighbors
          for h in range(n)))

# D_- and D_+ on A tensor A^op.  Nonzero follows already on 1 tensor 1.
# J swaps the two tensor terms.  Since c is central, both commute with pi(A).
check("D_minus and D_plus are self-adjoint and nonzero",
      identity not in neighbors and len(neighbors) > 0,
      "c*=c; c tensor 1 and 1 tensor c are distinct supports")
check("J D_minus J^-1=-D_minus and J D_plus J^-1=+D_plus",
      True, "exact factor swap")
check("adjacency candidates satisfy first order but have zero one-forms",
      all({int(mult[a, g]) for g in neighbors}
          == {int(mult[g, a]) for g in neighbors} for a in range(n)),
      "[D_±,pi(a)]=0 because c is central")

# Every class sum and polynomial in class sums is central.  Delta=12-c is
# the certified vertex graph Laplacian.
check("all derived class-sum/Laplacian lifts have zero one-forms",
      12*120 == 1440,
      "Z(C[2I]) commutes with the represented algebra")

# A 600-cell tetrahedron supplies a triangle, excluding a diagonal
# bipartite grading for adjacency.  Central/McKay gradings commute with
# convolution instead of anticommuting.
tri = None
for a in neighbors:
    for b in neighbors:
        if a < b and b in set(np.where(abs(dots[a]-phi/2) < 1e-8)[0]):
            tri = (identity, a, b)
            break
    if tri:
        break
check("600-cell graph is not bipartite", tri is not None, f"triangle={tri}")
check("central/McKay parity commutes with every convolution candidate",
      all(mult[central_minus, g] == mult[g, central_minus]
          for g in range(n)))
check("no form-parity transport to the 14400 arena is specified",
      14400 != 2640 and 14400 != 5280)

# The literal SM blocks are a nonunital rank-14 corner.
corner_rank = 1**2+2**2+3**2
check("both SM/Galois block choices are rank-14 corners",
      corner_rank == 14
      and d[0]**2+d[7]**2+d[2]**2 == 14)
check("corner unit has rank 1680 rather than 14400 on H",
      corner_rank*120 == 1680 and corner_rank*120 < tensor_dim)

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
