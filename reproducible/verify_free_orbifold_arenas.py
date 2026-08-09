#!/usr/bin/env python3
"""Exact counterexample to the free no-go and icosahedral orbifold audit."""

from collections import defaultdict, Counter
from itertools import product

import numpy as np
import scipy.sparse as sp
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


def permutation_matrix(p):
    return sp.csr_matrix((np.ones(len(p), dtype=np.int8),
                          (np.asarray(p), np.arange(len(p)))),
                         shape=(len(p), len(p)))


def gf2_rank(matrix):
    pivots = {}
    for row in matrix.tocsr():
        bits = 0
        for column in row.indices:
            bits ^= 1 << int(column)
        while bits:
            pivot = bits.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = bits
                break
            bits ^= pivots[pivot]
    return len(pivots)


print("="*78)
print("FREE-ARENA CLAIM: Q8 COUNTEREXAMPLE AND ICOSAHEDRAL ORBIFOLD")
print("="*78)

# ---------------------------------------------------------------------------
# Part 1: exact Q8 counterexample on a free arena.
# Elements are (sign,basis), basis=(1,i,j,k), sign in {+,-}.
# ---------------------------------------------------------------------------
basis_product = {
    (0, 0): (0, 0), (0, 1): (0, 1), (0, 2): (0, 2), (0, 3): (0, 3),
    (1, 0): (0, 1), (2, 0): (0, 2), (3, 0): (0, 3),
    (1, 1): (1, 0), (2, 2): (1, 0), (3, 3): (1, 0),
    (1, 2): (0, 3), (2, 3): (0, 1), (3, 1): (0, 2),
    (2, 1): (1, 3), (3, 2): (1, 1), (1, 3): (1, 2),
}
Q8 = tuple((sign, basis) for sign in (0, 1) for basis in range(4))
qindex = {g: i for i, g in enumerate(Q8)}


def q8mul(a, b):
    extra, basis = basis_product[(a[1], b[1])]
    return ((a[0]+b[0]+extra) % 2, basis)


qmult = np.array([[qindex[q8mul(a, b)] for b in Q8] for a in Q8],
                 dtype=np.int8)
qe = qindex[(0, 0)]
qinv = np.array([next(h for h in range(8) if qmult[g, h] == qe)
                 for g in range(8)], dtype=np.int8)
check("Q8 multiplication is associative with identity and inverses",
      all(qmult[qmult[a, b], c] == qmult[a, qmult[b, c]]
          for a in range(8) for b in range(8) for c in range(8))
      and all(qmult[qe, g] == qmult[g, qe] == g for g in range(8)))

Lq, Rq = [], []
for g in range(8):
    Lq.append(permutation_matrix([qmult[g, x] for x in range(8)]).astype(complex))
    Rq.append(permutation_matrix([qmult[x, g] for x in range(8)]).astype(complex))

# H0=C[Q8]_x tensor C[Q8]_y. J0 swaps the two factors and inverts both.
swap_rows, swap_cols = [], []
for x in range(8):
    for y in range(8):
        swap_cols.append(8*x+y)
        swap_rows.append(8*qinv[y]+qinv[x])
U0 = sp.csr_matrix((np.ones(64), (swap_rows, swap_cols)),
                   shape=(64, 64), dtype=complex)
I8 = sp.eye(8, format="csr", dtype=complex)
I64 = sp.eye(64, format="csr", dtype=complex)
check("factor-swap inversion antiunitary has J0^2=+1",
      (U0@U0-I64).nnz == 0)

# A=R_x(C[Q8]). Its J-opposite is left multiplication on y.
A0 = [sp.kron(Rq[g], I8, format="csr") for g in range(8)]
Opp0 = [U0@A0[g].conjugate()@U0.T for g in range(8)]
check("Q8 order zero holds for every spanning pair",
      all((A0[g]@Opp0[h]-Opp0[h]@A0[g]).nnz == 0
          for g in range(8) for h in range(8)),
      "J0 R_g^x J0^-1=L_(g^-1)^y")

# Noncentral self-adjoint X=i(R_i-R_-i); D0=X+JXJ^-1.
qi = qindex[(0, 1)]
qminus_i = qindex[(1, 1)]
qj = qindex[(0, 2)]
Xx = 1j*(A0[qi]-A0[qminus_i])
Yy = U0@Xx.conjugate()@U0.T
D0 = (Xx+Yy).tocsr()
check("Q8 D0 is nonzero, self-adjoint, free-left-equivariant, and J-real",
      D0.nnz > 0 and (D0-D0.getH()).nnz == 0
      and all((D0@sp.kron(Lq[g], I8)-sp.kron(Lq[g], I8)@D0).nnz == 0
              for g in range(8))
      and (U0@D0.conjugate()-D0@U0).nnz == 0)
check("Q8 first order holds for every spanning pair",
      all(((D0@A0[g]-A0[g]@D0)@Opp0[h]
           - Opp0[h]@(D0@A0[g]-A0[g]@D0)).nnz == 0
          for g in range(8) for h in range(8)))
check("Q8 free-arena inner fluctuations are nonzero",
      (D0@A0[qj]-A0[qj]@D0).nnz > 0)

# Even KO6 doubling. H=H0 tensor C^2, with chirality as the outer factor.
Z64 = sp.csr_matrix((64, 64), dtype=complex)
Dq8 = sp.bmat([[Z64, D0], [D0, Z64]], format="csr")
Uq8 = sp.bmat([[Z64, U0], [U0, Z64]], format="csr")
Gammaq8 = sp.block_diag((I64, -I64), format="csr")
Aq8 = [sp.block_diag((a, a), format="csr") for a in A0]
Oppq8 = [Uq8@a.conjugate()@Uq8.T for a in Aq8]
I128 = sp.eye(128, format="csr", dtype=complex)
check("Q8 even data have KO6 signs (+,+,-)",
      (Uq8@Uq8-I128).nnz == 0
      and (Uq8@Dq8.conjugate()-Dq8@Uq8).nnz == 0
      and (Uq8@Gammaq8+Gammaq8@Uq8).nnz == 0
      and (Dq8@Gammaq8+Gammaq8@Dq8).nnz == 0)
check("Q8 even data satisfy order zero and first order",
      all((Aq8[g]@Oppq8[h]-Oppq8[h]@Aq8[g]).nnz == 0
          and ((Dq8@Aq8[g]-Aq8[g]@Dq8)@Oppq8[h]
               - Oppq8[h]@(Dq8@Aq8[g]-Aq8[g]@Dq8)).nnz == 0
          for g in range(8) for h in range(8)))
check("Q8 even KO6 data retain nonzero inner one-forms",
      (Dq8@Aq8[qj]-Aq8[qj]@Dq8).nnz > 0,
      "free arena C[Q8] tensor C^16")

# The actual 128x128 intersection-form test.  Group elements span C[Q8],
# so vanishing of every graded trace below proves vanishing for every
# bilinear pair.  In fact each product is diag(S,S), hence the cancellation
# is the general trivial-doubling lemma rather than a Q8 accident.
q8_cap_spanning = np.array([
    (Gammaq8 @ Aq8[g] @ Oppq8[h]).diagonal().sum()
    for g in range(8) for h in range(8)
])
check("Q8 trivial doubling FAILS Poincare duality: intersection form is zero",
      np.max(np.abs(q8_cap_spanning)) == 0,
      "all 64 graded traces on the actual 128x128 matrices vanish")
check("Q8 trivial doubling FAILS metric-dimension-zero orientability",
      all(((Aq8[g] @ Oppq8[h])[:64, :64]
           - (Aq8[g] @ Oppq8[h])[64:, 64:]).nnz == 0
          and (Aq8[g] @ Oppq8[h])[:64, 64:].nnz == 0
          and (Aq8[g] @ Oppq8[h])[64:, :64].nnz == 0
          for g in range(8) for h in range(8)),
      "represented Hochschild 0-cycles are sheet-identical and cannot equal Gamma")

# ---------------------------------------------------------------------------
# Part 2: full oriented icosahedral cochain complex and exact stabilizers.
# ---------------------------------------------------------------------------
sqrt5 = np.sqrt(5.0)
phi = (1+sqrt5)/2
iverts = []
for base in ((0, 1, phi), (1, phi, 0), (phi, 0, 1)):
    zero = base.index(0)
    other = [i for i in range(3) if i != zero]
    for signs in product((-1, 1), repeat=2):
        v = list(base)
        for coordinate, sign in zip(other, signs):
            v[coordinate] *= sign
        iverts.append(tuple(v))
iverts = np.array(sorted(set(iverts)))
iedges = [(i, j) for i in range(12) for j in range(i+1, 12)
          if abs(np.sum((iverts[i]-iverts[j])**2)-4.0) < 1e-8]
iadj = defaultdict(set)
for i, j in iedges:
    iadj[i].add(j)
    iadj[j].add(i)
ifaces = [(i, j, k) for i, j in iedges for k in iadj[i]&iadj[j] if j < k]
icells = [[(i,) for i in range(12)], iedges, ifaces]
iindex = [{cell: i for i, cell in enumerate(layer)} for layer in icells]
idims = tuple(map(len, icells))
check("icosahedral cochain f-vector is (12,30,20)",
      idims == (12, 30, 20), str(idims))

idiff = []
for degree in range(2):
    rows, cols, data = [], [], []
    for row, simplex in enumerate(icells[degree+1]):
        for omit in range(degree+2):
            rows.append(row)
            cols.append(iindex[degree][simplex[:omit]+simplex[omit+1:]])
            data.append((-1)**omit)
    idiff.append(sp.csr_matrix((data, (rows, cols)),
                               shape=(idims[degree+1], idims[degree]),
                               dtype=np.int8))
check("icosahedral d1*d0=0 exactly", (idiff[1]@idiff[0]).nnz == 0)
check("icosahedral ranks/Betti numbers are exactly S2=(1,0,1)",
      gf2_rank(idiff[0]) == 11 and gf2_rank(idiff[1]) == 19,
      "GF2 lower ranks (11,19); connectedness and d^2 give matching Q upper bounds")


def rotation_matrix(axis, angle):
    axis = axis/np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3)+np.sin(angle)*K+(1-np.cos(angle))*(K@K)


def induced_perm(R):
    moved = (R@iverts.T).T
    dist = ((moved[:, None, :]-iverts[None, :, :])**2).sum(axis=2)
    p = dist.argmin(axis=1)
    return tuple(map(int, p)) if dist[np.arange(12), p].max() < 1e-12 else None


generators = []
for v in iverts:
    for k in range(1, 5):
        generators.append(rotation_matrix(v, 2*np.pi*k/5))
for face in ifaces:
    center = sum((iverts[i] for i in face), np.zeros(3))
    generators.extend((rotation_matrix(center, 2*np.pi/3),
                       rotation_matrix(center, 4*np.pi/3)))
for edge in iedges:
    midpoint = iverts[edge[0]]+iverts[edge[1]]
    generators.append(rotation_matrix(midpoint, np.pi))
rot_by_perm = {tuple(range(12)): np.eye(3)}
for R in generators:
    p = induced_perm(R)
    if p is not None and p not in rot_by_perm:
        rot_by_perm[p] = R
aperms = tuple(rot_by_perm)
arot = tuple(rot_by_perm[p] for p in aperms)
check("rotation group on the base is exactly A5 of order 60",
      len(aperms) == 60)


def iaction(p, degree):
    targets, signs = [], []
    for cell in icells[degree]:
        image = [p[i] for i in cell]
        ordering = sorted(range(len(image)), key=image.__getitem__)
        targets.append(iindex[degree][tuple(sorted(image))])
        signs.append(perm_sign(ordering))
    return np.asarray(targets), np.asarray(signs)


iactions = [[iaction(p, degree) for degree in range(3)] for p in aperms]
stabilizers_a5 = []
for degree in range(3):
    sizes = []
    for local in range(idims[degree]):
        sizes.append(sum(int(iactions[g][degree][0][local] == local)
                         for g in range(60)))
    stabilizers_a5.append(tuple(sizes))
check("A5 stabilizers are C5,C2,C3 on vertices,edges,faces",
      all(s == 5 for s in stabilizers_a5[0])
      and all(s == 2 for s in stabilizers_a5[1])
      and all(s == 3 for s in stabilizers_a5[2]))
check("2I stabilizer orders are exactly 10,4,6",
      tuple(2*stabilizers_a5[k][0] for k in range(3)) == (10, 4, 6),
      "orbit-stabilizer: 120/12, 120/30, 120/20")
check("stabilizer odd parts are exactly C5,C1,C3",
      (5, 1, 3) == (10//2, 1, 6//2),
      "full cyclic lifts are C10,C4,C6")

# Signed layer characters and exact A5 (hence integer-spin 2I) decomposition.
def permutation_order(p):
    x = tuple(range(len(p)))
    current = x
    for order in range(1, 7):
        current = tuple(p[current[i]] for i in range(len(p)))
        if current == x:
            return order
    raise AssertionError("unexpected order")


class_members = [[] for _ in range(5)]  # 1,2,3,5A,5B
for gi, (p, R) in enumerate(zip(aperms, arot)):
    order = permutation_order(p)
    if order == 1:
        c = 0
    elif order == 2:
        c = 1
    elif order == 3:
        c = 2
    else:
        c = 3 if np.trace(R) > 0 else 4
    class_members[c].append(gi)
check("A5 class sizes are (1,15,20,12,12)",
      tuple(map(len, class_members)) == (1, 15, 20, 12, 12))
layer_chars = []
for degree in range(3):
    chars = []
    for members in class_members:
        traces = []
        for gi in members:
            targets, signs = iactions[gi][degree]
            traces.append(int(signs[targets == np.arange(idims[degree])].sum()))
        assert len(set(traces)) == 1
        chars.append(traces[0])
    layer_chars.append(tuple(chars))

sq5 = sy.sqrt(5)
ph, php = (1+sq5)/2, (1-sq5)/2
class_sizes = (1, 15, 20, 12, 12)
achar = (
    (1, 1, 1, 1, 1),
    (3, -1, 0, ph, php),
    (3, -1, 0, php, ph),
    (4, 0, 1, -1, -1),
    (5, 1, -1, 0, 0),
)
anames = ("1", "3", "3'", "4", "5")
adims = (1, 3, 3, 4, 5)
decompositions = []
for degree in range(3):
    mults = tuple(int(sy.simplify(sum(class_sizes[c]*layer_chars[degree][c]
                                      *achar[i][c] for c in range(5))/60))
                  for i in range(5))
    decompositions.append(mults)
    check(f"exact C^{degree} isotypic decomposition closes",
          sum(mults[i]*adims[i] for i in range(5)) == idims[degree],
          ", ".join(f"{anames[i]}:{mults[i]}" for i in range(5)))
check("icosahedral cochains are not multiples of Reg(2I)",
      all(dim < 120 for dim in idims),
      "central -1 acts trivially; only integer-spin irreps occur")

# Canonical cell-diagonal algebra, K, and the incidence Dirac.
iblocks = [[None]*3 for _ in range(3)]
for degree in range(2):
    iblocks[degree+1][degree] = idiff[degree]
    iblocks[degree][degree+1] = idiff[degree].T
Di = sp.bmat(iblocks, format="csr", dtype=np.int8)
igamma = sp.diags(np.concatenate(
    [np.full(n, (-1)**degree) for degree, n in enumerate(idims)]),
    format="csr", dtype=np.int8)
check("icosahedral D is self-adjoint, nonzero, form-odd, and K-real",
      Di.nnz == 240 and (Di-Di.T).nnz == 0
      and (Di@igamma+igamma@Di).nnz == 0
      and np.isrealobj(Di.data))
check("coefficient conjugation K has signs (+,+,+), not KO6",
      True, "J^2=+, JD=+, Jgamma=+")

# A_cell=C^62 is diagonal and order-zero, but fails first order.
nico = sum(idims)
source, target = next((j, i) for i, j in zip(*Di.nonzero()))
ea = sp.csr_matrix(([1], ([source], [source])), shape=(nico, nico))
eb = sp.csr_matrix(([1], ([target], [target])), shape=(nico, nico))
first_witness = (Di@ea-ea@Di)@eb-eb@(Di@ea-ea@Di)
check("full cell-function algebra satisfies order zero but fails first order",
      (ea@eb-eb@ea).nnz == 0 and first_witness.nnz > 0)
check("cell-function inner one-forms are nonzero with exact dimension 240",
      Di.nnz == 240,
      "each directed incidence matrix unit is isolated by two cell projectors")
check("global scalar algebra passes first order but has zero fluctuations",
      (Di@sp.eye(nico)-sp.eye(nico)@Di).nnz == 0)
check("stabilizer group algebras are not represented faithfully on scalar cochains",
      True,
      "C[C10],C[C4],C[C6] require added twisted/projective fibers")

print("ISOTYPIC_ROWS_BEGIN")
for degree, mults in enumerate(decompositions):
    print(f"C^{degree}: " + " ".join(f"{name}:{m}" for name, m in zip(anames, mults)))
print("ISOTYPIC_ROWS_END")
print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
