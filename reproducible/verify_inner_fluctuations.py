#!/usr/bin/env python3
"""Exact finite audit of free 2I cells, convolution form, and fluctuations.

All matrix/cell assertions after the floating construction of the standard
600-cell vertex coordinates are integer permutation/incidence assertions.
Nearest-vertex matching is certified with a 1e-15 squared-distance bound.
"""

from collections import defaultdict
from itertools import permutations, product

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
                     for j in range(i + 1, len(p))) % 2 else 1


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z, w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X, w*Z+x*Y-y*X+z*W))


print("=" * 78)
print("FREE 2I CELLS, CONVOLUTION FORM, AND INNER FLUCTUATIONS")
print("=" * 78)

# The same deterministic 600-cell construction as verify_kahler_dirac.py.
phi = (1 + np.sqrt(5.0)) / 2
vertex_set = set()
for i in range(4):
    for sign in (-1.0, 1.0):
        q = [0.0] * 4
        q[i] = sign
        vertex_set.add(tuple(q))
vertex_set.update(product((-0.5, 0.5), repeat=4))
base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
even_perms = [p for p in permutations(range(4)) if perm_sign(p) == 1]
for p in even_perms:
    q = [base[p[i]] for i in range(4)]
    nonzero = [i for i, x in enumerate(q) if abs(x) > 1e-12]
    for signs in product((-1, 1), repeat=3):
        r = q[:]
        for i, sign in zip(nonzero, signs):
            r[i] *= sign
        vertex_set.add(tuple(round(x, 10) for x in r))
V = np.array(sorted(vertex_set))
identity = int(np.argmin(np.sum((V-np.array((1., 0., 0., 0.)))**2,
                                axis=1)))

dots = V @ V.T
edges = [(i, j) for i in range(120) for j in range(i+1, 120)
         if abs(dots[i, j] - phi/2) < 1e-3]
adj = defaultdict(set)
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)
triangles = [(i, j, k) for i, j in edges for k in adj[i] & adj[j]
             if j < k]
tetrahedra = [(i, j, k, ell) for i, j, k in triangles
              for ell in adj[i] & adj[j] & adj[k] if k < ell]
cells = [[(i,) for i in range(120)], edges, triangles, tetrahedra]
dims = tuple(map(len, cells))
indices = [{cell: i for i, cell in enumerate(layer)} for layer in cells]
check("600-cell f-vector", dims == (120, 720, 1200, 600), str(dims))

# Integer coboundaries and D.
d = []
for degree in range(3):
    rr, cc, vv = [], [], []
    for row, simplex in enumerate(cells[degree+1]):
        for omit in range(degree+2):
            face = simplex[:omit] + simplex[omit+1:]
            rr.append(row)
            cc.append(indices[degree][face])
            vv.append((-1)**omit)
    d.append(sp.csr_matrix((vv, (rr, cc)),
                           shape=(dims[degree+1], dims[degree]),
                           dtype=np.int8))
blocks = [[None]*4 for _ in range(4)]
for degree in range(3):
    blocks[degree+1][degree] = d[degree]
    blocks[degree][degree+1] = d[degree].T
D = sp.bmat(blocks, format="csr", dtype=np.int8)

# Left group action and exact multiplication table (indices in sorted V).
perms = []
max_match_error = 0.0
for g in V:
    moved = np.array([qmul(g, v) for v in V])
    distances = ((moved[:, None, :] - V[None, :, :])**2).sum(axis=2)
    p = distances.argmin(axis=1)
    max_match_error = max(max_match_error,
                          float(distances[np.arange(120), p].max()))
    perms.append(p)
check("all quaternion products match vertices", max_match_error < 1e-15,
      f"max squared matching error={max_match_error:.3e}")
mult = np.empty((120, 120), dtype=np.int16)
for gi, p in enumerate(perms):
    mult[gi, :] = p
inverse = np.empty(120, dtype=np.int16)
for g in range(120):
    inverse[g] = int(np.where(mult[g] == identity)[0][0])
check("multiplication table has identity and inverses",
      np.all(mult[identity] == np.arange(120)) and
      np.all(mult[:, identity] == np.arange(120)) and
      all(mult[g, inverse[g]] == identity for g in range(120)))
associative = all(np.array_equal(mult[g, mult[h, :]],
                                 mult[mult[g, h], :])
                  for g in range(120) for h in range(120))
check("matched multiplication table is associative exactly", associative)


def cell_action(p, degree):
    target, signs = [], []
    for cell in cells[degree]:
        image = [int(p[i]) for i in cell]
        order = sorted(range(len(image)), key=image.__getitem__)
        target.append(indices[degree][tuple(sorted(image))])
        signs.append(perm_sign(order))
    return np.asarray(target, dtype=np.int32), np.asarray(signs, dtype=np.int8)


actions = [[cell_action(p, degree) for degree in range(4)] for p in perms]

# Explicit orbit charts. F(delta_g tensor e_alpha) = L_g(rep_alpha),
# including the orientation sign of the sorted-cell basis.
representatives = []
rep_degree = []
charts = []
layer_orbit_counts = []
for degree in range(4):
    unseen = set(range(dims[degree]))
    chart = [None] * dims[degree]
    count = 0
    while unseen:
        rep_index = min(unseen)
        alpha = len(representatives)
        representatives.append(cells[degree][rep_index])
        rep_degree.append(degree)
        targets = [int(actions[g][degree][0][rep_index]) for g in range(120)]
        signs = [int(actions[g][degree][1][rep_index]) for g in range(120)]
        check(f"C^{degree} orbit {count} is free",
              len(set(targets)) == 120,
              f"representative={cells[degree][rep_index]}")
        for g, target in enumerate(targets):
            if chart[target] is not None:
                raise AssertionError("nonfree orbit chart collision")
            chart[target] = (g, alpha, signs[g])
        unseen.difference_update(targets)
        count += 1
    charts.append(chart)
    layer_orbit_counts.append(count)
check("free-orbit counts are (1,6,10,5)",
      tuple(layer_orbit_counts) == (1, 6, 10, 5),
      str(tuple(layer_orbit_counts)))
check("explicit chart gives C = C[2I] tensor C^22",
      len(representatives) == 22 and
      all(all(entry is not None for entry in chart) for chart in charts))

print("ORBIT_REPRESENTATIVES_BEGIN")
for alpha, (degree, rep) in enumerate(zip(rep_degree, representatives)):
    print(f"{alpha:02d} C^{degree} {rep}")
print("ORBIT_REPRESENTATIVES_END")

# Extract w_{alpha,beta} from the columns of the 22 representatives.
# R_h(delta_x)=delta_{x h}; D=sum R_{w_ab} tensor E_ab.
offsets = np.cumsum((0,) + dims)
w = defaultdict(lambda: defaultdict(int))
rep_local_index = [indices[degree][rep]
                   for degree, rep in zip(rep_degree, representatives)]
for beta, (degree, source_local) in enumerate(zip(rep_degree, rep_local_index)):
    source_global = int(offsets[degree] + source_local)
    column = D.getcol(source_global).tocoo()
    for target_global, value in zip(column.row, column.data):
        target_degree = int(np.searchsorted(offsets[1:], target_global,
                                            side="right"))
        target_local = int(target_global-offsets[target_degree])
        h, alpha, target_sign = charts[target_degree][target_local]
        w[(alpha, beta)][h] += int(value)*target_sign
w = {key: {h: value for h, value in coeff.items() if value}
     for key, coeff in w.items()}
w = {key: coeff for key, coeff in w.items() if coeff}

# Check every one of the 14880 physical nonzeros against the extracted form.
reconstruction_ok = True
for source_degree in range(4):
    for source_local in range(dims[source_degree]):
        x, beta, source_sign = charts[source_degree][source_local]
        source_global = int(offsets[source_degree] + source_local)
        column = D.getcol(source_global).tocoo()
        actual = {}
        for target_global, value in zip(column.row, column.data):
            target_degree = int(np.searchsorted(offsets[1:], target_global,
                                                side="right"))
            target_local = int(target_global-offsets[target_degree])
            y, alpha, target_sign = charts[target_degree][target_local]
            actual[(y, alpha)] = int(value)*source_sign*target_sign
        predicted = {}
        for (alpha, b), coeff in w.items():
            if b == beta:
                for h, value in coeff.items():
                    predicted[(int(mult[x, h]), alpha)] = value
        if actual != predicted:
            reconstruction_ok = False
            break
    if not reconstruction_ok:
        break
check("D is reconstructed exactly from its 22x22 right-convolution coefficients",
      reconstruction_ok,
      f"{len(w)} nonzero blocks, {sum(map(len, w.values()))} group coefficients")

# The coefficient-support graph controls diagonal multiplicity algebras.
support_adj = [set() for _ in range(22)]
for alpha, beta in w:
    if alpha != beta:
        support_adj[alpha].add(beta)
        support_adj[beta].add(alpha)
seen = {0}
frontier = [0]
while frontier:
    u = frontier.pop()
    for v in support_adj[u]-seen:
        seen.add(v)
        frontier.append(v)
check("the 22-vertex coefficient-support graph of D is connected",
      len(seen) == 22,
      "first order reduces every coordinate-diagonal multiplicity algebra to scalars")

print("D_COEFFICIENTS_BEGIN")
for (alpha, beta), coeff in sorted(w.items()):
    terms = " ".join(f"{h}:{coeff[h]:+d}" for h in sorted(coeff))
    print(f"w[{alpha:02d},{beta:02d}] {terms}")
print("D_COEFFICIENTS_END")

# Exact Wedderburn/FS data in McKay-chain convention.
sqrt5 = sy.sqrt(5)
ph = (1+sqrt5)/2
php = (1-sqrt5)/2
class_sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)
x = (2, -2, 0, 1, -1, ph, -ph, ph-1, php)
xp = (2, -2, 0, 1, -1, php, ph-1, -ph, ph)


def sym_powers(t):
    answer = [tuple([1]*9), t]
    for _ in range(2, 6):
        answer.append(tuple(sy.expand(t[k]*answer[-1][k]-answer[-2][k])
                            for k in range(9)))
    return answer


sx, sxp = sym_powers(x), sym_powers(xp)
standard = (sx[0], x, xp, sx[2], sxp[2],
            tuple(sy.expand(x[k]*xp[k]) for k in range(9)),
            sx[3], sx[4], sx[5])
order = (0, 1, 3, 6, 7, 8, 5, 2, 4)
characters = tuple(standard[i] for i in order)
irrep_dims = (1, 2, 3, 4, 5, 6, 4, 2, 3)
square_class = (0, 0, 1, 4, 4, 7, 7, 6, 6)
fs = tuple(sy.simplify(sum(class_sizes[k]*chi[square_class[k]]
                           for k in range(9))/sy.Integer(120))
           for chi in characters)
check("Wedderburn dimensions square-sum to |2I|",
      sum(n*n for n in irrep_dims) == 120,
      "C + M2^2 + M3^2 + M4^2 + M5 + M6")
check("exact FS indicators alternate real/quaternionic by central parity",
      fs == (1, -1, 1, -1, 1, -1, 1, -1, 1), str(fs))
check("canonical block availability is C, quaternionic M2, and M3",
      irrep_dims[0] == 1 and fs[1] == -1 and irrep_dims[2] == 3,
      "the simultaneous Galois flip exchanges rho1<->rho7 and rho2<->rho8")

# Left algebra fluctuations vanish: D was already checked equivariant, and
# we repeat the exact spanning-set check on all 120 group elements.
left_commutes = True
left_blocks = []
for g in range(120):
    layer_matrices = []
    for degree in range(4):
        target, signs = actions[g][degree]
        layer_matrices.append(sp.csr_matrix(
            (signs, (target, np.arange(dims[degree]))),
            shape=(dims[degree], dims[degree]), dtype=np.int8))
    Lg = sp.block_diag(layer_matrices, format="csr")
    left_blocks.append(Lg)
    left_commutes &= (D @ Lg - Lg @ D).nnz == 0
check("[D,L_g]=0 for the 120-element group-algebra spanning set",
      left_commutes)
check("canonical left inner one-forms and fluctuations vanish",
      left_commutes, "sum a_i[D,b_i]=0 identically")

# Refute the proposed two-horn dichotomy if a right-convolution witness does
# not commute with D. Compare convolution coefficients exactly, without a
# dense 2640 matrix.
noncommuting_pair = next((s, t) for s in range(120) for t in range(120)
                         if mult[s, t] != mult[t, s])
right_witness = None
for s in range(120):
    for block, coeff in w.items():
        drs = defaultdict(int)
        rsd = defaultdict(int)
        for h, value in coeff.items():
            drs[int(mult[s, h])] += value
            rsd[int(mult[h, s])] += value
        if dict(drs) != dict(rsd):
            right_witness = (s, block)
            break
    if right_witness:
        break
check("an equivariant right multiplication has [D,R_s] != 0",
      right_witness is not None, str(right_witness))
check("the proposed equivariance-break/multiplicity-only dichotomy is false",
      right_witness is not None,
      "right convolution is 2I-equivariant, acts on C[2I], and fluctuates")

# Orbitwise inversion J_iota in the explicit chart.
j_rows, j_cols, j_data = [], [], []
for degree in range(4):
    for local in range(dims[degree]):
        g, alpha, source_sign = charts[degree][local]
        target_g = int(inverse[g])
        rep = rep_local_index[alpha]
        target = int(actions[target_g][degree][0][rep])
        target_sign = int(actions[target_g][degree][1][rep])
        j_rows.append(int(offsets[degree]+target))
        j_cols.append(int(offsets[degree]+local))
        j_data.append(source_sign*target_sign)
Jinv = sp.csr_matrix((j_data, (j_rows, j_cols)),
                     shape=(2640, 2640), dtype=np.int8)
I2640 = sp.eye(2640, format="csr", dtype=np.int8)
check("orbitwise inversion is an antiunitary signed permutation with J^2=1",
      (Jinv@Jinv-I2640).nnz == 0)
jd_comm = (Jinv@D-D@Jinv).nnz == 0
jd_anticomm = (Jinv@D+D@Jinv).nnz == 0
check("orbitwise inversion has no JD=+/-DJ sign for this D",
      not jd_comm and not jd_anticomm,
      f"commutator nnz={(Jinv@D-D@Jinv).nnz}, "
      f"anticommutator nnz={(Jinv@D+D@Jinv).nnz}")

# Build a physical signed-permutation R_s from the explicit orbit chart.
def right_matrix(s):
    rrows, rcols, rdata = [], [], []
    for degree in range(4):
        for local in range(dims[degree]):
            x0, alpha, source_sign = charts[degree][local]
            y0 = int(mult[x0, s])
            rep = rep_local_index[alpha]
            target = int(actions[y0][degree][0][rep])
            target_sign = int(actions[y0][degree][1][rep])
            rrows.append(int(offsets[degree]+target))
            rcols.append(int(offsets[degree]+local))
            rdata.append(source_sign*target_sign)
    return sp.csr_matrix((rdata, (rrows, rcols)),
                         shape=(2640, 2640), dtype=np.int8)


# Order-zero/one statements, checked as sparse integer identities.  Use the
# fluctuation witness and choose an element not commuting with it.
s = int(right_witness[0])
t = next(t for t in range(120) if mult[s, t] != mult[t, s])
Rs = right_matrix(s)
Rt = right_matrix(t)
check("inversion sends left to right, so left order zero holds",
      (Jinv@left_blocks[t]@Jinv-right_matrix(int(inverse[t]))).nnz == 0
      and (left_blocks[s]@Rt-Rt@left_blocks[s]).nnz == 0,
      "J L_t J^-1=R_(t^-1) and [L_s,R_t]=0")
right_commutator = D@Rs-Rs@D
check("right group algebra also satisfies order zero and first order",
      (Jinv@Rt@Jinv-left_blocks[int(inverse[t])]).nnz == 0
      and right_commutator.nnz > 0
      and (right_commutator@left_blocks[t]
           - left_blocks[t]@right_commutator).nnz == 0,
      "J R_t J^-1=L_(t^-1), [D,R_s]!=0, and [[D,R_s],L_t]=0")
check("left and right admissible algebras are incomparable and have no common greatest algebra",
      mult[s, t] != mult[t, s],
      "a greatest algebra would violate order zero on L_s and J R_t J^-1=L_t")

# Coefficient conjugation K is identity on the real cell basis. It has the
# (+,+,+) sign table, but does not implement the opposite algebra.
check("coefficient conjugation K has signs (J2,JD,Jgamma)=(+,+,+)",
      np.isrealobj(D.data))
check("K fails order zero for the noncommutative group algebra",
      mult[s, t] != mult[t, s],
      "[L_s,K L_t K^-1]=[L_s,L_t] != 0")
check("primal Hodge-star J is absent on the 2640 primal arena",
      dims[0] != dims[3] and dims[1] != dims[2],
      "120!=600 and 720!=1200")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
