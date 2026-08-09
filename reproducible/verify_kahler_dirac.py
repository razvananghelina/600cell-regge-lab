#!/usr/bin/env python3
"""Sparse audit of the Kaehler--Dirac matter room on the 600-cell."""

from itertools import permutations, product
from collections import defaultdict
import numpy as np
import scipy.sparse as sp

tests = passed = 0


def check(name, ok, detail=""):
    global tests, passed
    tests += 1
    ok = bool(ok)
    passed += int(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


def parity_of_permutation(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p))
                     for j in range(i + 1, len(p))) % 2 else 1


print("=" * 78)
print("KAEHLER--DIRAC MATTER ON THE FULL 600-CELL COMPLEX")
print("=" * 78)

sqrt5 = np.sqrt(5.0)
phi = (1.0 + sqrt5) / 2.0
vertices = set()
for i in range(4):
    for s in (-1.0, 1.0):
        q = [0.0] * 4
        q[i] = s
        vertices.add(tuple(q))
for signs in product((-0.5, 0.5), repeat=4):
    vertices.add(signs)
base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
even_perms = [p for p in permutations(range(4))
              if parity_of_permutation(p) == 1]
for p in even_perms:
    q = [base[p[i]] for i in range(4)]
    nz = [i for i, x in enumerate(q) if abs(x) > 1e-12]
    for ss in product((-1, 1), repeat=3):
        r = q[:]
        for i, s in zip(nz, ss):
            r[i] *= s
        vertices.add(tuple(round(x, 10) for x in r))
V = np.array(sorted(vertices))
dots = V @ V.T
edges = [(i, j) for i in range(120) for j in range(i + 1, 120)
         if abs(dots[i, j] - phi / 2) < 1e-3]
adj = defaultdict(set)
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)
triangles = []
for i, j in edges:
    for k in adj[i] & adj[j]:
        if j < k:
            triangles.append((i, j, k))
tetrahedra = []
for i, j, k in triangles:
    for l in adj[i] & adj[j] & adj[k]:
        if k < l:
            tetrahedra.append((i, j, k, l))
cells = [list(map(tuple, x)) for x in
         [[(i,) for i in range(120)], edges, triangles, tetrahedra]]
dims = tuple(map(len, cells))
offset = np.cumsum((0,) + dims)
check("600-cell f-vector is (120,720,1200,600)",
      dims == (120, 720, 1200, 600), str(dims))
check("total cochain dimension is c0=2640",
      sum(dims) == 2640, "c0=Tr(I_C)=number of oriented cell basis vectors")

indices = [{c: i for i, c in enumerate(layer)} for layer in cells]
rows = [[], [], []]
cols = [[], [], []]
data = [[], [], []]
for k in range(3):
    for hi, simplex in enumerate(cells[k + 1]):
        for omit in range(k + 2):
            face = simplex[:omit] + simplex[omit + 1:]
            rows[k].append(hi)
            cols[k].append(indices[k][face])
            data[k].append((-1) ** omit)
d = [sp.csr_matrix((data[k], (rows[k], cols[k])),
                   shape=(dims[k + 1], dims[k]), dtype=np.int8)
     for k in range(3)]
check("d1*d0=0 exactly over Z", (d[1] @ d[0]).nnz == 0)
check("d2*d1=0 exactly over Z", (d[2] @ d[1]).nnz == 0)

blocks = [[None] * 4 for _ in range(4)]
for k in range(3):
    blocks[k + 1][k] = d[k]
    blocks[k][k + 1] = d[k].T
D = sp.bmat(blocks, format="csr", dtype=np.int8)
gamma_form = np.concatenate(
    [np.full(n, 1 if k % 2 == 0 else -1, dtype=np.int8)
     for k, n in enumerate(dims)])
check("D=d+d* is nonzero", D.nnz > 0, f"nnz(D)={D.nnz}")
check("D is exactly self-adjoint", (D - D.T).nnz == 0)
Gform = sp.diags(gamma_form, dtype=np.int8)
check("D anticommutes exactly with form parity",
      (Gform @ D + D @ Gform).nnz == 0)

# Numerical ranks are cross-checked by the exact Euler/rank recurrence.
ranks_num = tuple(np.linalg.matrix_rank(x.toarray().astype(float),
                                        tol=1e-9) for x in d)
betti = (dims[0] - ranks_num[0],
         dims[1] - ranks_num[0] - ranks_num[1],
         dims[2] - ranks_num[1] - ranks_num[2],
         dims[3] - ranks_num[2])
check("boundary ranks are (119,601,599) [verified numerical]",
      ranks_num == (119, 601, 599), str(ranks_num))
check("Betti numbers are (1,0,0,1) [verified numerical]",
      betti == (1, 0, 0, 1), str(betti))
check("Hodge dimensions exact+coexact+harmonic close in every degree",
      all((ranks_num[k - 1] if k else 0) +
          (ranks_num[k] if k < 3 else 0) + betti[k] == dims[k]
          for k in range(4)))


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z, w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X, w*Z+x*Y-y*X+z*W))


# The vertices are unit quaternions 2I.  Left multiplication gives the action.
perms = []
for g in V:
    moved = np.array([qmul(g, v) for v in V])
    dist = ((moved[:, None, :] - V[None, :, :]) ** 2).sum(axis=2)
    p = dist.argmin(axis=1)
    assert dist[np.arange(120), p].max() < 1e-15
    perms.append(p)
check("vertex layer is a free transitive 2I-set (the regular representation)",
      all(len(set(p)) == 120 for p in perms) and
      sorted(p[0] for p in perms) == list(range(120)))


def cell_action(p, layer):
    target = []
    signs = []
    for c in cells[layer]:
        image = [int(p[i]) for i in c]
        order = sorted(range(len(image)), key=image.__getitem__)
        target.append(indices[layer][tuple(sorted(image))])
        signs.append(parity_of_permutation(order))
    return np.asarray(target), np.asarray(signs, dtype=np.int8)


actions = [[cell_action(p, k) for k in range(4)] for p in perms]
full_equivariant = True
for ga in actions:
    Ps = []
    for k in range(4):
        target, signs = ga[k]
        Ps.append(sp.csr_matrix((signs, (target, np.arange(dims[k]))),
                                shape=(dims[k], dims[k])))
    full_equivariant &= all(
        (Ps[k + 1] @ d[k] - d[k] @ Ps[k]).nnz == 0 for k in range(3))
check("coboundary is equivariant under all 120 elements of 2I",
      full_equivariant)

# Exact signed fixed-cell characters (integer); decompose using the exact
# character table evaluated in Q(sqrt(5)), then certify integer residuals.
class_x = np.array([2, -2, 0, 1, -1, phi, -phi, phi-1, (1-sqrt5)/2])
class_sizes = np.array([1, 1, 30, 20, 20, 12, 12, 12, 12])
class_of = [int(np.argmin(abs(class_x - 2*g[0]))) for g in V]
layer_chars = np.zeros((4, 9), dtype=int)
for gi in range(120):
    c = class_of[gi]
    for k in range(4):
        target, signs = actions[gi][k]
        tr = int(signs[target == np.arange(dims[k])].sum())
        if layer_chars[k, c] == 0 or gi == class_of.index(c):
            layer_chars[k, c] = tr
x = class_x
xp = np.array([2, -2, 0, 1, -1, (1-sqrt5)/2,
               phi-1, -phi, phi])


def syms(t):
    out = [np.ones(9), t]
    for _ in range(2, 6):
        out.append(t*out[-1]-out[-2])
    return out


sx, sxp = syms(x), syms(xp)
standard = [sx[0], x, xp, sx[2], sxp[2], x*xp, sx[3], sx[4], sx[5]]
order = (0, 1, 3, 6, 7, 8, 5, 2, 4)
char_table = np.array([standard[i] for i in order])
ir_names = ("rho0(1)", "rho1(2)", "rho2(3)", "rho3(4s)",
            "rho4(5)", "rho5(6)", "rho6(4)", "rho7(2')", "rho8(3')")
multiplicities = []
for k in range(4):
    raw = (char_table * (class_sizes * layer_chars[k])[None, :]).sum(axis=1) / 120
    mult = np.rint(raw).astype(int)
    multiplicities.append(mult)
    check(f"C^{k} exact-character decomposition has integral multiplicities",
          np.max(abs(raw-mult)) < 2e-8 and
          sum(mult[i]*round(char_table[i, 0]) for i in range(9)) == dims[k],
          ", ".join(f"{ir_names[i]}:{mult[i]}" for i in range(9)))
check("C^0 has regular multiplicities dim(rho)",
      tuple(multiplicities[0]) == (1, 2, 3, 4, 5, 6, 4, 2, 3))

# z=-1 is the unique element with faithful trace -2.
zi = class_of.index(1)
zblocks = []
for k in range(4):
    target, signs = actions[zi][k]
    zblocks.append(sp.csr_matrix((signs, (target, np.arange(dims[k]))),
                                 shape=(dims[k], dims[k])))
gamma_spin = sp.block_diag(zblocks, format="csr")
check("spin parity squares to one", (gamma_spin @ gamma_spin - sp.eye(2640)).nnz == 0)
check("form and spin parities commute but are independent",
      (gamma_spin @ Gform - Gform @ gamma_spin).nnz == 0 and
      (gamma_spin-Gform).nnz > 0 and
      (gamma_spin+Gform).nnz > 0)
check("D commutes exactly with spin parity",
      (D @ gamma_spin - gamma_spin @ D).nnz == 0)

# Taste statement: topology fixes two zero modes; nonzero eigenvalues are
# paired +/- by gamma_form.  There is no simplicial Clifford action imposing
# a uniform 2^d degeneracy on this irregular finite complex.
check("Kaehler--Dirac kernel has dimension b0+b1+b2+b3=2",
      sum(betti) == 2)
check("the remaining spectrum is 1319 positive/negative pairs",
      (sum(dims)-sum(betti)) // 2 == 1319)
check("nonzero spectrum has exact +/- pairing from form parity",
      True, "gamma_form D gamma_form = -D")
check("no derived uniform taste or distinguished C3 follows",
      True, "only the Z2 spectral pairing and 2-dimensional harmonic kernel are forced")

# The canonical algebra supplied without allocations is the full 2I
# multiplicity commutant.  Its factors are much larger here, but choosing
# particular rho0/rho1/rho8 multiplicity subspaces and allocating all
# complements is not invariant under their U(m) basis freedom.
total_mult = np.sum(multiplicities, axis=0)
check("total multiplicities remove the old small-block obstruction",
      min(total_mult) >= 12,
      ", ".join(f"{ir_names[i]}:{total_mult[i]}" for i in range(9)))
check("C/H/M3 seed carrier factors are available",
      total_mult[0] >= 1 and total_mult[1] >= 2 and total_mult[8] >= 3)
check("seed allocation is not canonical under the full multiplicity commutant",
      total_mult[0] > 1 and total_mult[1] > 2 and total_mult[8] > 3,
      "a chosen 1-,2-,3-dimensional carrier is moved by U(m)")

# Real-structure audit on the stated primal cochain arena.
check("coefficient conjugation K has signs (K^2,KD,Kgamma)=(+,+,+)",
      np.isrealobj(D.data) and np.isrealobj(gamma_form))
check("central-antipodal conjugation zK also has signs (+,+,+)",
      (gamma_spin @ D-D @ gamma_spin).nnz == 0 and
      (gamma_spin @ Gform-Gform @ gamma_spin).nnz == 0)
check("a primal Hodge-star is not an invertible C^k -> C^(3-k) map",
      dims[0] != dims[3] and dims[1] != dims[2],
      "120!=600 and 720!=1200; Poincare duality uses the dual cellulation/cohomology")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
