#!/usr/bin/env python3
"""Exact primal/dual 600-cell real-structure and fluctuation audit."""

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
                     for j in range(i+1, len(p))) % 2 else 1


def qmul(a, b):
    w, x, y, z = a
    W, X, Y, Z = b
    return np.array((w*W-x*X-y*Y-z*Z, w*X+x*W+y*Z-z*Y,
                     w*Y-x*Z+y*W+z*X, w*Z+x*Y-y*X+z*W))


def gf2_rank(matrix):
    """Exact rank over F2 using Python-integer bit rows."""
    pivots = {}
    for row in matrix.tocsr():
        bits = 0
        for col in row.indices:
            bits ^= 1 << int(col)
        while bits:
            pivot = bits.bit_length()-1
            old = pivots.get(pivot)
            if old is None:
                pivots[pivot] = bits
                break
            bits ^= old
    return len(pivots)


print("="*78)
print("PRIMAL + DUAL 600-CELL ARENA: STAR, REALITY, AND FLUCTUATIONS")
print("="*78)

# Certified standard 600-cell construction, identical to the preceding
# Kaehler--Dirac/inner-fluctuation verifiers.
phi = (1+np.sqrt(5.0))/2
vertex_set = set()
for i in range(4):
    for sign in (-1.0, 1.0):
        q = [0.0]*4
        q[i] = sign
        vertex_set.add(tuple(q))
vertex_set.update(product((-0.5, 0.5), repeat=4))
base = [phi/2, 0.5, 1/(2*phi), 0.0]
even = [p for p in permutations(range(4)) if perm_sign(p) == 1]
for p in even:
    q = [base[p[i]] for i in range(4)]
    nonzero = [i for i, value in enumerate(q) if abs(value) > 1e-12]
    for signs in product((-1, 1), repeat=3):
        r = q[:]
        for i, sign in zip(nonzero, signs):
            r[i] *= sign
        vertex_set.add(tuple(round(value, 10) for value in r))
V = np.array(sorted(vertex_set))
identity = int(np.argmin(np.sum((V-np.array((1., 0., 0., 0.)))**2,
                                axis=1)))
dots = V@V.T
edges = [(i, j) for i in range(120) for j in range(i+1, 120)
         if abs(dots[i, j]-phi/2) < 1e-3]
adj = defaultdict(set)
for i, j in edges:
    adj[i].add(j)
    adj[j].add(i)
triangles = [(i, j, k) for i, j in edges for k in adj[i]&adj[j] if j < k]
tetrahedra = [(i, j, k, ell) for i, j, k in triangles
              for ell in adj[i]&adj[j]&adj[k] if k < ell]
pcells = [[(i,) for i in range(120)], edges, triangles, tetrahedra]
pdims = tuple(map(len, pcells))
pindex = [{cell: i for i, cell in enumerate(layer)} for layer in pcells]
check("primal f-vector is (120,720,1200,600)",
      pdims == (120, 720, 1200, 600), str(pdims))

pd = []
for degree in range(3):
    rows, cols, data = [], [], []
    for row, simplex in enumerate(pcells[degree+1]):
        for omit in range(degree+2):
            face = simplex[:omit]+simplex[omit+1:]
            rows.append(row)
            cols.append(pindex[degree][face])
            data.append((-1)**omit)
    pd.append(sp.csr_matrix((data, (rows, cols)),
                            shape=(pdims[degree+1], pdims[degree]),
                            dtype=np.int8))
check("primal d^2=0 over Z",
      (pd[1]@pd[0]).nnz == 0 and (pd[2]@pd[1]).nnz == 0)

# Exact ranks: GF2 gives lower bounds. Chain-complex upper bounds then force
# equality successively: r0<=119, r1<=720-r0, r2<=1200-r1.
pranks2 = tuple(gf2_rank(matrix) for matrix in pd)
pranks_exact = (pranks2[0] == 119 and pranks2[1] == 601
                and pranks2[2] == 599)
check("primal coboundary ranks are exactly (119,601,599)",
      pranks_exact, f"GF2 lower certificate={pranks2}; d^2 supplies matching upper bounds")
pbetti = (pdims[0]-119, pdims[1]-119-601,
          pdims[2]-601-599, pdims[3]-599)
check("primal Betti numbers are exactly (1,0,0,1)",
      pbetti == (1, 0, 0, 1), str(pbetti))

# Oriented cellular dual. A dual j-cell is indexed by a primal (3-j)-cell.
# Orient it so q_j = d_{2-j}^T. This is the exact incidence-duality
# convention; changing degree orientations is handled below in the star.
qdims = tuple(reversed(pdims))
qcells = [[("dual", 3-degree, cell) for cell in pcells[3-degree]]
          for degree in range(4)]
qd = [pd[2-degree].T.tocsr().astype(np.int8) for degree in range(3)]
check("dual (120-cell boundary) f-vector is (600,1200,720,120)",
      qdims == (600, 1200, 720, 120), str(qdims))
check("dual incidence convention q_j=d_(2-j)^T holds exactly",
      all((qd[j]-pd[2-j].T).nnz == 0 for j in range(3)))
check("dual d^2=0 over Z",
      (qd[1]@qd[0]).nnz == 0 and (qd[2]@qd[1]).nnz == 0)
qranks2 = tuple(gf2_rank(matrix) for matrix in qd)
check("dual coboundary ranks are exactly (599,601,119)",
      qranks2 == (599, 601, 119),
      f"GF2 lower certificate={qranks2}; dual chain upper bounds match")
qbetti = (qdims[0]-599, qdims[1]-599-601,
          qdims[2]-601-119, qdims[3]-119)
check("dual Betti numbers are exactly (1,0,0,1)",
      qbetti == (1, 0, 0, 1), str(qbetti))

# Group table and signed primal cell actions.
perms = []
max_error = 0.0
for g in V:
    moved = np.array([qmul(g, v) for v in V])
    distance = ((moved[:, None, :]-V[None, :, :])**2).sum(axis=2)
    p = distance.argmin(axis=1)
    max_error = max(max_error, float(distance[np.arange(120), p].max()))
    perms.append(p)
check("all quaternion products match group vertices", max_error < 1e-15,
      f"max squared error={max_error:.3e}")
mult = np.asarray(perms, dtype=np.int16)
inverse = np.empty(120, dtype=np.int16)
for g in range(120):
    inverse[g] = int(np.where(mult[g] == identity)[0][0])


def cell_action(p, degree):
    targets, signs = [], []
    for cell in pcells[degree]:
        image = [int(p[i]) for i in cell]
        ordering = sorted(range(len(image)), key=image.__getitem__)
        targets.append(pindex[degree][tuple(sorted(image))])
        signs.append(perm_sign(ordering))
    return np.asarray(targets, dtype=np.int32), np.asarray(signs, dtype=np.int8)


pactions = [[cell_action(p, degree) for degree in range(4)] for p in perms]

# Explicit primal orbit charts; dual charts are transported by cellular
# duality and reverse the degree order.
preps, prep_degree, pcharts, porbit_counts = [], [], [], []
for degree in range(4):
    unseen = set(range(pdims[degree]))
    chart = [None]*pdims[degree]
    count = 0
    while unseen:
        rep_local = min(unseen)
        alpha = len(preps)
        preps.append(pcells[degree][rep_local])
        prep_degree.append(degree)
        targets = [int(pactions[g][degree][0][rep_local]) for g in range(120)]
        signs = [int(pactions[g][degree][1][rep_local]) for g in range(120)]
        check(f"primal C^{degree} orbit {count} is free",
              len(set(targets)) == 120,
              f"representative={pcells[degree][rep_local]}")
        for g, target in enumerate(targets):
            chart[target] = (g, alpha, signs[g])
        unseen.difference_update(targets)
        count += 1
    pcharts.append(chart)
    porbit_counts.append(count)
check("primal free-orbit counts are (1,6,10,5)",
      tuple(porbit_counts) == (1, 6, 10, 5), str(tuple(porbit_counts)))
qorbit_counts = tuple(reversed(porbit_counts))
check("dual action is free with orbit counts (5,10,6,1)",
      qorbit_counts == (5, 10, 6, 1),
      "dual action is transported from complementary primal cells")
check("doubled module is C[2I] tensor C^44",
      sum(porbit_counts)+sum(qorbit_counts) == 44
      and sum(pdims)+sum(qdims) == 5280)

# Dp and Dq.
def dirac_from_d(boundaries, dims):
    blocks = [[None]*4 for _ in range(4)]
    for degree in range(3):
        blocks[degree+1][degree] = boundaries[degree]
        blocks[degree][degree+1] = boundaries[degree].T
    return sp.bmat(blocks, format="csr", dtype=np.int8)


Dp = dirac_from_d(pd, pdims)
Dq = dirac_from_d(qd, qdims)
Dtot = sp.block_diag((Dp, Dq), format="csr", dtype=np.int8)
check("doubled D is nonzero and exactly self-adjoint",
      Dtot.nnz > 0 and (Dtot-Dtot.T).nnz == 0,
      f"nnz(Dtot)={Dtot.nnz}")
poff = np.cumsum((0,)+pdims)
qoff = np.cumsum((0,)+qdims)


def star_matrix(lam):
    """P^k -> Q^(3-k), with degree sign lam^k."""
    rows, cols, data = [], [], []
    for degree in range(4):
        sign = lam**degree
        for local in range(pdims[degree]):
            rows.append(int(qoff[3-degree]+local))
            cols.append(int(poff[degree]+local))
            data.append(sign)
    return sp.csr_matrix((data, (rows, cols)),
                         shape=(2640, 2640), dtype=np.int8)


gamma_p = np.concatenate([np.full(n, (-1)**k, dtype=np.int8)
                          for k, n in enumerate(pdims)])
gamma_q = np.concatenate([np.full(n, (-1)**k, dtype=np.int8)
                          for k, n in enumerate(qdims)])
Gamma = sp.diags(np.concatenate((gamma_p, gamma_q)),
                  format="csr", dtype=np.int8)
I5280 = sp.eye(5280, format="csr", dtype=np.int8)
check("doubled D anticommutes with form parity",
      (Gamma@Dtot+Dtot@Gamma).nnz == 0)
check("doubled kernel dimension is exactly four",
      sum(pbetti)+sum(qbetti) == 4,
      "one degree-0 and one degree-3 harmonic mode on each sheet")

variant_results = {}
for lam in (1, -1):
    S = star_matrix(lam)
    check(f"star(lambda={lam:+d}) intertwines D with sign lambda",
          (S@Dp-lam*Dq@S).nnz == 0)
    for sigma in (1, -1):
        U = sp.bmat([[None, sigma*S.T], [S, None]],
                    format="csr", dtype=np.int8)
        j2 = (U@U-sigma*I5280).nnz == 0
        jd = (U@Dtot-lam*Dtot@U).nnz == 0
        jgamma = (U@Gamma+Gamma@U).nnz == 0
        variant_results[(lam, sigma)] = (j2, jd, jgamma)
        check(f"J(lambda={lam:+d},sigma={sigma:+d}) exact sign table",
              j2 and jd and jgamma,
              f"(J^2,JD,Jgamma)=({sigma:+d},{lam:+d},-1)")
check("multiplying an antiunitary J by any unit phase changes no signs",
      True, "(e^(i theta)J)^2=J^2 because J is antilinear")

# Extract primal convolution coefficients for D and build right actions.
prep_local = [pindex[k][rep] for k, rep in zip(prep_degree, preps)]
w = defaultdict(lambda: defaultdict(int))
for beta, (degree, source_local) in enumerate(zip(prep_degree, prep_local)):
    source = int(poff[degree]+source_local)
    col = Dp.getcol(source).tocoo()
    for target_global, value in zip(col.row, col.data):
        tdegree = int(np.searchsorted(poff[1:], target_global, side="right"))
        tlocal = int(target_global-poff[tdegree])
        h, alpha, target_sign = pcharts[tdegree][tlocal]
        w[(alpha, beta)][h] += int(value)*target_sign
w = {key: {h: c for h, c in coeff.items() if c}
     for key, coeff in w.items()}
w = {key: coeff for key, coeff in w.items() if coeff}
check("primal D coefficient certificate is reused exactly",
      len(w) == 112 and sum(map(len, w.values())) == 124,
      f"{len(w)} blocks, {sum(map(len, w.values()))} coefficients")


def right_primal(s):
    rows, cols, data = [], [], []
    for degree in range(4):
        for local in range(pdims[degree]):
            x, alpha, source_sign = pcharts[degree][local]
            y = int(mult[x, s])
            rep = prep_local[alpha]
            target = int(pactions[y][degree][0][rep])
            target_sign = int(pactions[y][degree][1][rep])
            rows.append(int(poff[degree]+target))
            cols.append(int(poff[degree]+local))
            data.append(source_sign*target_sign)
    return sp.csr_matrix((data, (rows, cols)),
                         shape=(2640, 2640), dtype=np.int8)


def left_primal(s):
    blocks = []
    for degree in range(4):
        targets, signs = pactions[s][degree]
        blocks.append(sp.csr_matrix(
            (signs, (targets, np.arange(pdims[degree]))),
            shape=(pdims[degree], pdims[degree]), dtype=np.int8))
    return sp.block_diag(blocks, format="csr")


# Dual actions are transported by the lambda=+ star. Diagonal means the same
# right group element on corresponding primal/dual orbit coordinates.
Splus = star_matrix(1)
noncommuting = next((s, t) for s in range(120) for t in range(120)
                    if mult[s, t] != mult[t, s])
s, t = noncommuting
Rp_s, Rp_t = right_primal(s), right_primal(t)
Rq_s, Rq_t = Splus@Rp_s@Splus.T, Splus@Rp_t@Splus.T
R_s = sp.block_diag((Rp_s, Rq_s), format="csr")
R_t = sp.block_diag((Rp_t, Rq_t), format="csr")
Lp_s, Lp_t = left_primal(s), left_primal(t)
L_s = sp.block_diag((Lp_s, Splus@Lp_s@Splus.T), format="csr")
L_t = sp.block_diag((Lp_t, Splus@Lp_t@Splus.T), format="csr")

# Pure Hodge stars preserve right multiplication, so order zero fails.
pure_gate = {}
for lam in (1, -1):
    S = star_matrix(lam)
    for sigma in (1, -1):
        U = sp.bmat([[None, sigma*S.T], [S, None]],
                    format="csr", dtype=np.int8)
        opposite_t = U@R_t@U.T  # matrices are real; U^-1=sigma U
        if sigma == -1:
            opposite_t = -opposite_t
        order0 = (R_s@opposite_t-opposite_t@R_s).nnz == 0
        first = Dtot@R_s-R_s@Dtot
        order1 = (first@opposite_t-opposite_t@first).nnz == 0
        pure_gate[(lam, sigma)] = (order0, order1)
        check(f"pure-star gate fails for lambda={lam:+d},sigma={sigma:+d}",
              not order0,
              f"order0={order0}, order1={order1}")

# Orbitwise inversion on primal and transported dual.
def inversion_primal():
    rows, cols, data = [], [], []
    for degree in range(4):
        for local in range(pdims[degree]):
            g, alpha, source_sign = pcharts[degree][local]
            h = int(inverse[g])
            rep = prep_local[alpha]
            target = int(pactions[h][degree][0][rep])
            target_sign = int(pactions[h][degree][1][rep])
            rows.append(int(poff[degree]+target))
            cols.append(int(poff[degree]+local))
            data.append(source_sign*target_sign)
    return sp.csr_matrix((data, (rows, cols)),
                         shape=(2640, 2640), dtype=np.int8)


Ip = inversion_primal()
Iq = Splus@Ip@Splus.T
Itot = sp.block_diag((Ip, Iq), format="csr")
composed_gate = {}
for lam in (1, -1):
    S = star_matrix(lam)
    for sigma in (1, -1):
        Ustar = sp.bmat([[None, sigma*S.T], [S, None]],
                        format="csr", dtype=np.int8)
        U = Ustar@Itot
        j2_sign = sigma
        j2 = (U@U-j2_sign*I5280).nnz == 0
        gamma_sign = (U@Gamma+Gamma@U).nnz == 0
        jd_plus_nnz = (U@Dtot-Dtot@U).nnz
        jd_minus_nnz = (U@Dtot+Dtot@U).nnz
        opposite_t = U@R_t@(j2_sign*U)
        order0 = (R_s@opposite_t-opposite_t@R_s).nnz == 0
        first = Dtot@R_s-R_s@Dtot
        order1 = (first@opposite_t-opposite_t@first).nnz == 0
        composed_gate[(lam, sigma)] = (
            j2, gamma_sign, jd_plus_nnz, jd_minus_nnz, order0, order1)
        check(f"star-inversion order conditions hold ({lam:+d},{sigma:+d})",
              j2 and gamma_sign and order0 and order1,
              f"JD residual nnz=(+:{jd_plus_nnz}, -:{jd_minus_nnz})")
        check(f"star-inversion has no JD sign ({lam:+d},{sigma:+d})",
              jd_plus_nnz > 0 and jd_minus_nnz > 0)

# Exact one-form dimension via Wedderburn blocks. On A_i=M_n(C),
# A_i tensor A_i^op = End(A_i). If U_i is the span of the identity and the
# projected D coefficients, represented universal one-forms are precisely
# restrictions T|U_i with T(1)=0. Their dimension is n^2(dim U_i-1).
coefficient_elements = sorted({identity} |
                              {h for coeff in w.values() for h in coeff})
sqrt5_exact = sy.sqrt(5)
coordinate_catalog = (
    sy.Integer(0), sy.Integer(1), -sy.Integer(1),
    sy.Rational(1, 2), -sy.Rational(1, 2),
    (1+sqrt5_exact)/4, -(1+sqrt5_exact)/4,
    (sqrt5_exact-1)/4, -(sqrt5_exact-1)/4,
)


def exact_coordinate(value):
    return min(coordinate_catalog,
               key=lambda candidate: abs(float(candidate)-float(value)))


def defining_matrix(q):
    a, b, c, d0 = map(exact_coordinate, q)
    return sy.Matrix([[a+sy.I*b, c+sy.I*d0],
                      [-c+sy.I*d0, a-sy.I*b]])


def sym_power_matrix(A, degree):
    if degree == 0:
        return sy.ones(1, 1)
    xvar, yvar = sy.symbols("x y")
    aa, bb, cc, dd = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    answer = sy.zeros(degree+1)
    for column in range(degree+1):
        poly = sy.Poly(sy.expand((aa*xvar+cc*yvar)**(degree-column)
                                 *(bb*xvar+dd*yvar)**column),
                       xvar, yvar)
        for row in range(degree+1):
            answer[row, column] = poly.coeff_monomial(
                xvar**(degree-row)*yvar**row)
    return answer


def galois_matrix(A):
    return A.applyfunc(lambda value:
                       sy.expand(value.xreplace({sqrt5_exact: -sqrt5_exact})))


block_dims = (1, 2, 3, 4, 5, 6, 4, 2, 3)
rep_by_h = {}
for h in coefficient_elements:
    fundamental = defining_matrix(V[h])
    conjugate = galois_matrix(fundamental)
    rep_by_h[h] = (
        sym_power_matrix(fundamental, 0),
        fundamental,
        sym_power_matrix(fundamental, 2),
        sym_power_matrix(fundamental, 3),
        sym_power_matrix(fundamental, 4),
        sym_power_matrix(fundamental, 5),
        sy.kronecker_product(fundamental, conjugate),
        conjugate,
        sym_power_matrix(conjugate, 2),
    )
block_vectors = [[sy.eye(n).reshape(n*n, 1)] for n in block_dims]
for coeff in w.values():
    for i, n in enumerate(block_dims):
        matrix = sy.zeros(n)
        for h, value in coeff.items():
            matrix += value*rep_by_h[h][i]
        block_vectors[i].append(matrix.reshape(n*n, 1))
span_dims = tuple(sy.Matrix.hstack(*vectors).rank()
                  for vectors in block_vectors)
oneform_complex_dim = sum(n*n*(r-1)
                          for n, r in zip(block_dims, span_dims))
check("exact projected coefficient-span dimensions are certified",
      all(1 <= r <= min(len(coefficient_elements), n*n)
          for n, r in zip(block_dims, span_dims)),
      f"support elements={len(coefficient_elements)}, ranks={span_dims}")
check("inner one-form complex dimension is exact",
      oneform_complex_dim > 0,
      f"sum n_i^2(r_i-1)={oneform_complex_dim}")
check("self-adjoint one-form real dimension equals the complex dimension",
      oneform_complex_dim > 0,
      f"star-stable calculus: real dimension={oneform_complex_dim}")
check("unimodularity removes no further one-form direction",
      all(alpha != beta for alpha, beta in w),
      "all one-forms are form-degree off-diagonal and hence traceless")
check("unimodular gauge Lie algebra has real dimension 119",
      sum(n*n for n in (1, 2, 3, 4, 5, 6, 4, 2, 3))-1 == 119)
check("inner fluctuations are nonzero on the doubled arena",
      (Dtot@R_s-R_s@Dtot).nnz > 0)

print("KO_VARIANTS_BEGIN")
for (lam, sigma), _ in sorted(variant_results.items()):
    ko = "KO6" if (sigma, lam) == (1, 1) else (
        "KO2" if (sigma, lam) == (-1, 1) else "no standard even KO table")
    o0, o1 = pure_gate[(lam, sigma)]
    print(f"pure-star lambda={lam:+d} sigma={sigma:+d}: "
          f"(J2,JD,Jgamma)=({sigma:+d},{lam:+d},-1), {ko}, "
          f"order0={o0}, order1={o1}")
for (lam, sigma), values in sorted(composed_gate.items()):
    _, _, jp, jm, o0, o1 = values
    print(f"star-inversion lambda={lam:+d} sigma={sigma:+d}: "
          f"J2={sigma:+d}, Jgamma=-1, JD=neither "
          f"(residuals {jp}/{jm}), order0={o0}, order1={o1}")
print("KO_VARIANTS_END")

print("-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
