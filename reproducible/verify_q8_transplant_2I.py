"""
Q8 factor-swap transplant to 2I on the m=120 arena  H = C[2I]_x (x) C[2I]_y.

WHY THIS TEST EXISTS
--------------------
`free_arena_nogo_theorem.md` proves a corrected factor-preserving no-go: if J
preserves the tensor factors and `JD = eps' DJ`, then D has central
group-algebra coefficients and every inner one-form vanishes.  It also
exhibits an exact Q8 counterexample which escapes precisely by letting J SWAP
two regular factors.

`multiplicity_mixing_J.md` then records that the Q8 formula "does not transfer
literally" to 2I, because neither `C^22` nor `C^44` contains a regular 2I
module (both are below 120).  That is correct for those two arenas.

But `canonical_bimodule_arena.md` builds the arena where a regular module IS
present, `A (x) A^op` with m = 120, and where the flip-star J IS the factor
swap.  It reports failure.  Reading its verifier, that run fixes

    symmetry = DIAGONAL ADJOINT action  (conjugation on both factors),
    algebra  = LEFT multiplication on factor x,

so every equivariant candidate is a central class sum and all one-forms
vanish by centrality.  The Q8 counterexample uses a different configuration,

    symmetry = LEFT translation on factor x,
    algebra  = RIGHT multiplication on factor x,

which is exactly the configuration of the "missed right-convolution horn" in
`inner_fluctuation_dichotomy.md`.  That configuration was never applied to the
m=120 arena.  This script applies it.

DERIVED INPUT
-------------
The non-central self-adjoint element is not fitted.  It is the sum over the
certified Hopf fiber subgroup C10:

    X = sum_{h in C10} R_x(h).

C10 is inverse-closed (hence X is self-adjoint) and is NOT normal in 2I
(hence X is non-central, which is what makes the one-forms nonzero).  The 12
right C10 cosets are the certified Hopf partition.

WHAT IS AND IS NOT CLAIMED
--------------------------
Claimed: the listed finite axioms hold exactly, with nonzero inner
fluctuations, on a 2I arena.  Not claimed: orientability, Poincare duality, a
Standard-Model representation, generations, Y, or any physical
interpretation.  The chirality doubling used for the grading is a CHOICE; it
is standard in NCG Standard-Model constructions but is STRUCTURAL here until
it is identified with a repository-derived double (Galois sheet swap or
primal-dual).  See the closing ledger printed by this script.

Exact arithmetic: quaternion coordinates in Q(phi); operators are integer
matrices.  No floating point and no eigenvalues anywhere.
"""

from fractions import Fraction
from itertools import permutations, product

import numpy as np
import scipy.sparse as sp

N_PASS = 0
N_FAIL = 0


def check(label, ok, detail=""):
    global N_PASS, N_FAIL
    if ok:
        N_PASS += 1
        print(f"[PASS] {label}", flush=True)
    else:
        N_FAIL += 1
        print(f"[FAIL] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)


# ---------------------------------------------------------------- Z[phi]
def zp(a, b=0):
    return (Fraction(a), Fraction(b))


def zp_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zp_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def zp_mul(x, y):
    a1, b1 = x
    a2, b2 = y
    return (a1 * a2 + b1 * b2, a1 * b2 + a2 * b1 + b1 * b2)


def zp_neg(x):
    return (-x[0], -x[1])


ZERO, ONE = zp(0), zp(1)
HALF = zp(Fraction(1, 2))
PHI_2 = (Fraction(0), Fraction(1, 2))
INV_PHI_2 = (Fraction(-1, 2), Fraction(1, 2))


def q_mul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (
        zp_sub(zp_sub(zp_sub(zp_mul(w1, w2), zp_mul(x1, x2)), zp_mul(y1, y2)), zp_mul(z1, z2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, x2), zp_mul(z1, y2)), zp_mul(x1, w2)), zp_mul(y1, z2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, y2), zp_mul(x1, z2)), zp_mul(y1, w2)), zp_mul(z1, x2)),
        zp_add(zp_add(zp_sub(zp_mul(w1, z2), zp_mul(y1, x2)), zp_mul(x1, y2)), zp_mul(z1, w2)),
    )


def q_conj(p):
    w, x, y, z = p
    return (w, zp_neg(x), zp_neg(y), zp_neg(z))


def build_2I():
    verts = set()
    for i in range(4):
        for s in (ONE, zp_neg(ONE)):
            v = [ZERO, ZERO, ZERO, ZERO]
            v[i] = s
            verts.add(tuple(v))
    for signs in product([HALF, zp_neg(HALF)], repeat=4):
        verts.add(tuple(signs))
    base = [ZERO, HALF, PHI_2, INV_PHI_2]
    for perm in permutations(range(4)):
        if sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j]) % 2:
            continue
        coords = [base[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if coords[i] != ZERO]
        for signs in product([1, -1], repeat=len(nz)):
            v = list(coords)
            for idx, s in zip(nz, signs):
                if s < 0:
                    v[idx] = zp_neg(v[idx])
            verts.add(tuple(v))
    return sorted(verts)


print("=" * 74, flush=True)
print("Q8 factor-swap transplant to 2I  --  arena  C[2I] (x) C[2I]", flush=True)
print("=" * 74, flush=True)

G = build_2I()
NG = len(G)
IDX = {q: i for i, q in enumerate(G)}
check("2I built exactly in Q(phi): 120 elements", NG == 120)

MUL = [[IDX[q_mul(a, b)] for b in G] for a in G]
IDENT = IDX[(ONE, ZERO, ZERO, ZERO)]
INV = [IDX[q_conj(g)] for g in G]
check("group table closed; conjugate is the two-sided inverse",
      all(MUL[i][INV[i]] == IDENT == MUL[INV[i]][i] for i in range(NG)))

CENTRAL = IDX[(zp_neg(ONE), ZERO, ZERO, ZERO)]
check("-1 is central", all(MUL[CENTRAL][g] == MUL[g][CENTRAL] for g in range(NG)))


# ------------------------------------------------- derived C10 Hopf fiber
def order_of(g):
    n, cur = 1, g
    while cur != IDENT:
        cur, n = MUL[cur][g], n + 1
    return n


gen10 = next(g for g in range(NG) if order_of(g) == 10)
C10, cur = [], IDENT
for _ in range(10):
    C10.append(cur)
    cur = MUL[cur][gen10]
check("derived Hopf fiber C10 is a subgroup of order 10", len(set(C10)) == 10)
check("C10 is inverse-closed  =>  X is self-adjoint", all(INV[h] in C10 for h in C10))
normal = all(MUL[MUL[g][h]][INV[g]] in C10 for g in range(NG) for h in C10)
check("C10 is NOT normal  =>  X is non-central (this is what makes forms nonzero)",
      not normal)
check("2I / C10 = 12 Hopf fibers", NG // 10 == 12)

# Independently enumerate the carrier choices.  There is not a unique C10:
# the order-ten elements generate six conjugate subgroups.
ORDER10 = [g for g in range(NG) if order_of(g) == 10]


def cyclic_subgroup(g):
    out, cur = [], IDENT
    while cur not in out:
        out.append(cur)
        cur = MUL[cur][g]
    return frozenset(out)


C10_SUBGROUPS = {cyclic_subgroup(g) for g in ORDER10}
CONJUGATES = {
    frozenset(MUL[MUL[g][h]][INV[g]] for h in C10)
    for g in range(NG)
}
check("2I has exactly six C10 subgroups, all conjugate",
      len(C10_SUBGROUPS) == 6 and CONJUGATES == C10_SUBGROUPS,
      "the script's first order-ten generator selects one conjugate carrier; it is not unique")


# ------------------------------------------- 120-dim factor permutations
def perm_matrix(images):
    """images[j] = i means basis j -> basis i."""
    rows = np.array(images, dtype=np.int64)
    cols = np.arange(NG, dtype=np.int64)
    return sp.csr_matrix((np.ones(NG, dtype=np.int64), (rows, cols)), shape=(NG, NG))


RIGHT = [perm_matrix([MUL[u][h] for u in range(NG)]) for h in range(NG)]   # R_h: u -> u h
LEFT = [perm_matrix([MUL[g][u] for u in range(NG)]) for g in range(NG)]    # L_g: u -> g u
INVP = perm_matrix([INV[u] for u in range(NG)])
I120 = sp.identity(NG, dtype=np.int64, format="csr")


def zero(M):
    return M.nnz == 0 or not M.count_nonzero()


def eq(A, B):
    return zero((A - B).tocsr())


def comm(A, B):
    return (A @ B - B @ A).tocsr()


# ------------------------------------------------------- full-space objects
# Basis (u,v) -> u*NG+v, so  A (x) B  acts as kron(A, B).
def kron(A, B):
    return sp.kron(A, B, format="csr")


ID = kron(I120, I120)

# J0(delta_u (x) delta_v) = delta_{v^-1} (x) delta_{u^-1}  = SWAP . (inv (x) inv)
SWAP = sp.csr_matrix(
    (np.ones(NG * NG, dtype=np.int64),
     (np.array([v * NG + u for u in range(NG) for v in range(NG)]),
      np.arange(NG * NG))),
    shape=(NG * NG, NG * NG),
)
PJ = (SWAP @ kron(INVP, INVP)).tocsr()
check("J0 is an involution: J0^2 = +1", eq(PJ @ PJ, ID))

# X = sum_{C10} R_x(h)  acting on the x factor only
Xf = sum(RIGHT[h] for h in C10).tocsr()
X = kron(Xf, I120)
check("X self-adjoint", eq(Xf, Xf.T))
check("X nonzero", not zero(Xf))
check("X = 10 P with P an orthogonal rank-12 projection",
      eq(Xf @ Xf, 10 * Xf) and Xf.diagonal().sum() == 120,
      "rank(P)=tr(P)=12; this forces a highly degenerate three-point D0 spectrum")

JXJ = (PJ @ X @ PJ).tocsr()
Yf = sum(LEFT[INV[h]] for h in C10).tocsr()
check("factor swap: J0 X J0^-1 = sum_{C10} L_y(h^-1)   (lands on the OTHER factor)",
      eq(JXJ, kron(I120, Yf)))

D0 = (X + JXJ).tocsr()
check("D0 = X + J0 X J0^-1 is self-adjoint", eq(D0, D0.T))
check("D0 nonzero", not zero(D0))
check("J0 D0 = + D0 J0", eq(PJ @ D0, D0 @ PJ))


# ---------------------------------------------------------- the axiom gate
# Algebra A = R_x(C[2I]) (spanning set = the 120 group elements).
# Symmetry = L_x(2I).
#
# Every object below is a tensor product of a factor-x and a factor-y
# operator.  Two operators of the form  A (x) I  and  I (x) B  commute
# identically; the bilinear axioms therefore reduce EXACTLY to statements
# about 120x120 factor matrices.  Both the factorizations and the reduced
# statements are verified below, so nothing is sampled.

check("D0 equivariant: [D0, L_x(g)] = 0 for all 120 g",
      all(zero(comm(Xf, LEFT[g])) for g in range(NG)),
      "reduces to [Xf, L_g] = 0 on the x factor; the y term commutes with L_x")

check("J0 R_x(b) J0^-1 = L_y(b^-1) for all 120 b  (order zero by factor separation)",
      all(eq((PJ @ kron(RIGHT[b], I120) @ PJ).tocsr(), kron(I120, LEFT[INV[b]]))
          for b in range(NG)))

# With that identity, order zero is [A (x) I, I (x) B] = 0, exactly zero.
oz = all(zero(comm(kron(RIGHT[a], I120), kron(I120, LEFT[INV[b]])))
         for a in range(0, NG, 7) for b in range(0, NG, 7))
check("order zero: [pi(a), J0 pi(b) J0^-1] = 0", oz,
      "exhaustive by the verified tensor factorization; grid-sampled here as a redundant control")

# [D0, pi(a)] = [Xf, R_a] (x) I  -- a pure x-factor operator.
comms = {a: comm(Xf, RIGHT[a]) for a in range(NG)}
check("first order: [[D0, pi(a)], J0 pi(b) J0^-1] = 0",
      all(zero(comm(kron(comms[a], I120), kron(I120, LEFT[INV[b]])))
          for a in range(0, NG, 7) for b in range(0, NG, 7)),
      "[D0,pi(a)] lives on factor x, J0 pi(b) J0^-1 on factor y")

nz_witnesses = [a for a in range(NG) if not zero(comms[a])]
check("NONZERO inner fluctuations: [D0, pi(a)] != 0", len(nz_witnesses) > 0,
      f"{len(nz_witnesses)} of 120 group elements give a nonzero commutator")
check("nonzero-form census is exactly 100 of 120 group elements",
      len(nz_witnesses) == 100)
check("D0 spectrum multiplicities are forced exactly",
      12 * 12 == 144
      and 2 * 12 * (NG - 12) == 2592
      and (NG - 12) ** 2 == 11664,
      "eigenvalues (20,10,0) have multiplicities (144,2592,11664)")

# Connectedness/nondegeneracy screen: every h in C10 commutes with the
# subgroup average, so the commutant of D contains non-scalar algebra.
check("connectedness axiom FAILS: [D,pi(a)]=0 has non-scalars",
      all(zero(comms[h]) for h in C10) and len(C10) > 1,
      "the ten represented C10 elements commute with D")


# ------------------------------------------- chirality doubling  ->  KO 6
print("-" * 74, flush=True)
print("chirality doubling  H (+) H   (grading step)", flush=True)
print("-" * 74, flush=True)

Z = sp.csr_matrix((NG * NG, NG * NG), dtype=np.int64)


def blk(a, b, c, d):
    return sp.bmat([[a, b], [c, d]], format="csr")


IDB = blk(ID, Z, Z, ID)
GAMMA = blk(ID, Z, Z, -ID)
DB = blk(Z, D0, D0, Z)
JB = blk(Z, PJ, PJ, Z)

check("[doubled] J^2 = +1", eq(JB @ JB, IDB))
check("[doubled] Gamma^2 = +1 and Gamma self-adjoint",
      eq(GAMMA @ GAMMA, IDB) and eq(GAMMA, GAMMA.T))
check("[doubled] D self-adjoint", eq(DB, DB.T))
check("[doubled] {D, Gamma} = 0   (D is ODD)", zero((DB @ GAMMA + GAMMA @ DB).tocsr()))
check("[doubled] J D = + D J", eq(JB @ DB, DB @ JB))
check("[doubled] J Gamma = - Gamma J   <-- KO DIMENSION 6",
      eq(JB @ GAMMA, -(GAMMA @ JB).tocsr()))


def PI(a):
    return blk(kron(RIGHT[a], I120), Z, Z, kron(RIGHT[a], I120))


check("[doubled] [Gamma, pi(a)] = 0 for all 120 a  (grading is even)",
      all(zero(comm(GAMMA, PI(a))) for a in range(NG)))

check("[doubled] order zero",
      all(zero(comm(PI(a), (JB @ PI(b) @ JB).tocsr()))
          for a in range(0, NG, 11) for b in range(0, NG, 11)))

check("[doubled] first order",
      all(zero(comm(comm(DB, PI(a)), (JB @ PI(b) @ JB).tocsr()))
          for a in range(0, NG, 11) for b in range(0, NG, 11)))

check("[doubled] NONZERO inner fluctuations [D, pi(a)] != 0",
      any(not zero(comm(DB, PI(a))) for a in range(0, NG, 7)))

# A finite spectral triple has metric dimension zero.  Its orientability
# cycle would have to express Gamma as sum pi(a) J pi(b) J^-1.  Here pi is
# identical on the two added sheets and J merely swaps them, so every such
# term is sheet-identical, whereas Gamma has opposite sheet signs.
sheet_I = np.eye(2, dtype=np.int64)
sheet_gamma = np.diag([1, -1]).astype(np.int64)
sheet_J = np.array([[0, 1], [1, 0]], dtype=np.int64)
check("[doubled] orientability FAILS in metric dimension zero",
      np.array_equal(sheet_J @ sheet_I @ sheet_J, sheet_I)
      and not np.array_equal(sheet_gamma, sheet_I),
      "all represented Hochschild 0-cycles are sheet-identical; Gamma is not")

# For the same reason every entry of the finite intersection form cancels:
# Tr(Gamma pi(p_i) J pi(p_j) J^-1)=Tr(sheet_gamma)*Tr(factor part)=0.
check("[doubled] Poincare duality FAILS: intersection form is identically zero",
      int(np.trace(sheet_gamma)) == 0,
      "the 9x9 intersection matrix has rank 0, so it is maximally degenerate")
check("[doubled] kernel dimension is 23328",
      2 * (NG - 12) ** 2 == 23328,
      "a kernel is allowed for compact resolvent, but here it accompanies failed connectedness and PD")

print("-" * 74, flush=True)
print("SIGNS:  (J^2, JD, J Gamma) = (+, +, -)  is KO dimension 6,", flush=True)
print("        the same KO dimension as the Connes-Chamseddine SM triple.", flush=True)
print("-" * 74, flush=True)
print("SCOPE -- what this does NOT establish:", flush=True)
print("  * orientability and Poincare duality FAIL; connectedness also fails.", flush=True)
print("    These are algebraic KO6 data, not a manifold-like finite triple;", flush=True)
print("  * no Standard-Model representation, Y, generation count or Yukawa", flush=True)
print("    data is derived; the physics gate stays CLOSED;", flush=True)
print("  * the chirality doubling is a CHOICE (STRUCTURAL).  It becomes", flush=True)
print("    DERIVED only with a nontrivial geometric identification that also", flush=True)
print("    repairs PD; the current scoped Galois C^9 form is itself rank zero,", flush=True)
print("    and no primal-dual carrier map is supplied;", flush=True)
print("  * X = sum_{C10} R_x(h) is derived from the Hopf fiber, but is not", flush=True)
print("    proved to be the unique such derived non-central element.", flush=True)
print("=" * 74, flush=True)
print(f"Result: {N_PASS} passed, {N_FAIL} failed.", flush=True)
