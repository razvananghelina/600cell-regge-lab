"""STEP 1 of the (2,3,5) orbifold index route: enumerate the canonical
incidence operators and record kernel/cokernel characters.

NO TARGET COMPARISON IS PERFORMED IN THIS FILE.  M16, M15 and "three
generations" are never mentioned outside this docstring, and the enumeration
below cannot branch on them.  The comparison is a separate, later step, so that
the git history proves the operator list preceded the target.

WHY THE INDEX IS NOT THE OBSERVABLE
-----------------------------------
For any linear map T : A -> B, [A] = [ker] + [im] and [B] = [im] + [coker], so
index(T) = [ker] - [coker] = [A] - [B] identically -- independent of T.  So an
index character carries no information about the operator; it only restates a
virtual identity among induced characters, and those form an affine sublattice
of rank 11 here.  The operator-dependent content is ker and coker SEPARATELY,
which is what this file records.

WHAT MAKES AN OPERATOR CANONICAL
--------------------------------
Vertices, edges and faces are 2I/C10, 2I/C4, 2I/C6 (12, 30, 20).  The
equivariant Hom space between two induced line modules has one dimension per
double coset, and every relevant intersection H ∩ xKx^-1 is exactly the centre
{+-1}, so the Mackey condition is only that the two characters agree on the
centre.  Geometric incidence is a single double coset -- verified here, not
assumed -- so the incidence operator is the unique (up to one scalar) element
supported on that coset.  Its five (or three, or nine) non-incidence components
are exactly zero.

Exact arithmetic for the group; the linear algebra is integer/rational.
"""
import io
import contextlib
import json
import runpy
from fractions import Fraction

import numpy as np

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


# ---------------------------------------------------------------- the group
Q = "verify_q8_transplant_2I.py"
g = {"__name__": "__inc__", "__file__": Q}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(Q).read(), Q, "exec"), g)

G, MUL, INV, IDENT, NG, CENTRAL = g["G"], g["MUL"], g["INV"], g["IDENT"], g["NG"], g["CENTRAL"]
PHI = (1 + 5 ** 0.5) / 2
check("2I built exactly: 120 elements with a central involution", NG == 120)


def order_of(x):
    o, c = 1, x
    while c != IDENT:
        c, o = MUL[c][x], o + 1
    return o


def cyclic(x):
    out, c = [], IDENT
    while c not in out:
        out.append(c)
        c = MUL[c][x]
    return out                      # ordered as powers of x


gen10 = g["gen10"]
gen4 = next(x for x in range(NG) if order_of(x) == 4)
gen6 = next(x for x in range(NG) if order_of(x) == 6)
SUBS = {"C10": (gen10, cyclic(gen10)), "C4": (gen4, cyclic(gen4)), "C6": (gen6, cyclic(gen6))}
for name, (_, H) in SUBS.items():
    check(f"{name} is a subgroup of order {len(H)} containing the centre",
          len(set(H)) == int(name[1:]) and CENTRAL in H)


# ------------------------------------------------------ geometric incidence
def rot(qf):
    w, x, y, z = qf
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


QF = np.array([[float(a) + float(b) * PHI for a, b in q] for q in G])
ROT = [rot(QF[i]) for i in range(NG)]


def axis(i):
    w, v = np.linalg.eig(ROT[i])
    a = np.real(v[:, int(np.argmin(np.abs(w - 1.0)))])
    return a / np.linalg.norm(a)


AXIS = {"C10": axis(gen10), "C4": axis(gen4), "C6": axis(gen6)}


def cosets(H):
    """Right cosets Hx, with one representative each."""
    reps, seen = [], set()
    for x in range(NG):
        if x in seen:
            continue
        seen |= {MUL[h][x] for h in H}
        reps.append(x)
    return reps


REPS = {n: cosets(H) for n, (_, H) in SUBS.items()}
for n in REPS:
    check(f"{n} has {NG // len(SUBS[n][1])} cosets", len(REPS[n]) == NG // len(SUBS[n][1]))

POINTS = {n: [ROT[INV[r]] @ AXIS[n] for r in REPS[n]] for n in REPS}
check("cell counts are the icosahedron f-vector (12, 30, 20)",
      (len(POINTS["C10"]), len(POINTS["C4"]), len(POINTS["C6"])) == (12, 30, 20))


def double_cosets(H, K):
    table, seen, idx = {}, set(), 0
    for x in range(NG):
        if x in seen:
            continue
        d = {MUL[MUL[h][x]][k] for h in H for k in K}
        for y in d:
            table[y] = idx
        seen |= d
        idx += 1
    return table, idx


def incidence_coset(nh, nk, n_incident):
    """The single double coset carrying maximal-cosine (= incident) pairs.

    Convention trap: for right cosets Hx and Ky under the right G-action the
    invariant is the double coset of x*y^-1, NOT x^-1*y.  The wrong one
    scatters incidence across every coset.
    """
    H, K = SUBS[nh][1], SUBS[nk][1]
    tab, ndc = double_cosets(H, K)
    A, B = POINTS[nh], POINTS[nk]
    cosv = np.array([[float(np.dot(a, b)) for b in B] for a in A])
    best = np.round(cosv, 6).max()
    inc = np.abs(cosv - best) < 1e-6
    tally = {}
    for i, ra in enumerate(REPS[nh]):
        for j, rb in enumerate(REPS[nk]):
            d = tab[MUL[ra][INV[rb]]]
            t = tally.setdefault(d, [0, 0])
            t[0] += 1
            t[1] += int(inc[i, j])
    pure = [d for d, (tot, ic) in tally.items() if ic == tot and tot]
    mixed = [d for d, (tot, ic) in tally.items() if 0 < ic < tot]
    check(f"{nh}/{nk} incidence is exactly one double coset of {ndc}, none mixed",
          len(pure) == 1 and not mixed and int(inc.sum()) == n_incident,
          f"incident pairs={int(inc.sum())} (expected {n_incident}), coset={pure}")
    return tab, pure[0]


INC = {
    ("C10", "C4"): incidence_coset("C10", "C4", 60),    # each edge has 2 vertices
    ("C10", "C6"): incidence_coset("C10", "C6", 60),    # each face has 3 vertices
    ("C4", "C6"): incidence_coset("C4", "C6", 60),      # each face has 3 edges
}


# -------------------------------------------------- induced modules and maps
def characters_of(name):
    """The |H| characters of the cyclic subgroup, as complex power tables."""
    gen, H = SUBS[name]
    m = len(H)
    pos = {h: k for k, h in enumerate(H)}
    return m, pos, [np.exp(2j * np.pi * q * np.arange(m) / m) for q in range(m)]


def induced_rep(name, q):
    """Matrices of Ind_H^G(chi_q) on basis f_i supported on H r_i."""
    m, pos, chars = characters_of(name)
    chi = chars[q]
    reps = REPS[name]
    n = len(reps)
    idx = {}
    for i, r in enumerate(reps):
        for h in SUBS[name][1]:
            idx[MUL[h][r]] = (i, pos[h])
    mats = []
    for gg in range(NG):
        M = np.zeros((n, n), dtype=complex)
        for i, r in enumerate(reps):
            j, k = idx[MUL[r][INV[gg]]]
            M[j, i] = np.conj(chi[k])
        mats.append(M)
    return mats


def central_parity(name, q):
    m, pos, chars = characters_of(name)
    return int(np.round(np.real(chars[q][pos[CENTRAL]])))


def incidence_operator(nh, qh, nk, qk):
    """The canonical map Ind_K(chi_qk) -> Ind_H(chi_qh), supported on the
    incidence double coset alone.  Returns None if the central parities differ,
    in which case the whole Hom space is zero."""
    if central_parity(nh, qh) != central_parity(nk, qk):
        return None
    tab, d0 = INC[(nh, nk)]
    mh, posh, ch = characters_of(nh)
    mk, posk, ck = characters_of(nk)
    A, B = REPS[nh], REPS[nk]
    T = np.zeros((len(A), len(B)), dtype=complex)
    Hs, Ks = SUBS[nh][1], SUBS[nk][1]
    for i, ra in enumerate(A):
        for j, rb in enumerate(B):
            s = 0j
            for h in Hs:
                for k in Ks:
                    y = MUL[MUL[h][MUL[ra][INV[rb]]]][k]
                    if tab[y] == d0:
                        s += np.conj(ch[qh][posh[h]]) * ck[qk][posk[k]]
            T[i, j] = s
    return T


def char_of_projector(P, mats):
    return np.array([np.trace(P @ M) for M in mats])


# Character table built here, from the SU(2) restriction, and CHECKED
# orthonormal before use.  Importing another script's table and re-deriving its
# class order was the source of a wrong decomposition on the first run.
def su2_char(n, w):
    th = np.arccos(np.clip(w, -1.0, 1.0))
    s = np.sin(th)
    safe = np.where(np.abs(s) < 1e-12, 1.0, s)
    return np.where(np.abs(s) < 1e-12,
                    (n + 1) * np.cos(th) ** n,
                    np.sin((n + 1) * th) / safe)


W = np.array([float(a) + float(b) * PHI for a, b in [q[0:1][0:1][0] and q[0] or q[0] for q in G]]) \
    if False else np.array([float(q[0][0]) + float(q[0][1]) * PHI for q in G])
WBAR = np.array([float(q[0][0]) + float(q[0][1]) * (1 - 5 ** 0.5) / 2 for q in G])
base = {n: su2_char(n, W) for n in range(6)}
CH = {"1": base[0], "2": base[1], "3": base[2], "4": base[3], "5": base[4], "6": base[5],
      "2p": su2_char(1, WBAR), "3p": su2_char(2, WBAR)}
CH["4p"] = base[1] * base[5] - base[4] - CH["3p"]
NAMES = ["1", "2", "2p", "3", "3p", "4", "4p", "5", "6"]
CLASS_CHAR = [CH[n] for n in NAMES]
DIMS = [int(round(CH[n][IDENT].real if np.iscomplexobj(CH[n]) else CH[n][IDENT])) for n in NAMES]

check("character table is orthonormal (built and verified here)",
      all(abs(float(np.sum(CLASS_CHAR[a] * np.conj(CLASS_CHAR[b]))) / NG
              - (1.0 if a == b else 0.0)) < 1e-8
          for a in range(9) for b in range(9)),
      f"dims in order = {DIMS}")


def decompose(chi_vals):
    return [int(np.round(np.real(np.sum(chi_vals * np.conj(CLASS_CHAR[r])) / NG)))
            for r in range(9)]


# ------------------------------------------------------------- enumeration
records = []
pairs = [("C10", "C4"), ("C10", "C6"), ("C4", "C6")]
for nh, nk in pairs:
    mh = len(SUBS[nh][1])
    mk = len(SUBS[nk][1])
    for qh in range(mh):
        for qk in range(mk):
            T = incidence_operator(nh, qh, nk, qk)
            if T is None:
                continue
            if np.allclose(T, 0):
                records.append(dict(target=nh, tq=qh, source=nk, sq=qk,
                                    status="zero operator"))
                continue
            matsH = induced_rep(nh, qh)
            matsK = induced_rep(nk, qk)
            equi = max(np.abs(matsH[y] @ T - T @ matsK[y]).max() for y in range(NG))
            u, s, vh = np.linalg.svd(T)
            tol = max(T.shape) * s.max() * 1e-12
            rank = int((s > tol).sum())
            Pker = vh[rank:].conj().T @ vh[rank:]
            Pcok = u[:, rank:] @ u[:, rank:].conj().T
            kc = decompose(char_of_projector(Pker, matsK))
            cc = decompose(char_of_projector(Pcok, matsH))
            records.append(dict(
                target=nh, tq=qh, source=nk, sq=qk,
                equivariance_residual=float(equi),
                dim_source=T.shape[1], dim_target=T.shape[0], rank=rank,
                dim_ker=T.shape[1] - rank, dim_coker=T.shape[0] - rank,
                ker_char=kc, coker_char=cc, status="ok"))

ok = [r for r in records if r.get("status") == "ok"]
check("every canonical operator is exactly equivariant",
      all(r["equivariance_residual"] < 1e-8 for r in ok),
      f"max residual = {max((r['equivariance_residual'] for r in ok), default=0):.2e}")
check("kernel and cokernel decompose into integer multiplicities",
      all(sum(m * d for m, d in zip(r["ker_char"], DIMS)) == r["dim_ker"]
          and sum(m * d for m, d in zip(r["coker_char"], DIMS)) == r["dim_coker"]
          for r in ok))

N = len(ok)
distinct_ker = {tuple(r["ker_char"]) for r in ok}
print()
print(f"N (canonical nonzero incidence operators) = {N}")
print(f"distinct kernel characters  = {len(distinct_ker)}")
print(f"distinct cokernel characters= {len({tuple(r['coker_char']) for r in ok})}")
print(f"irrep dimensions in order   = {DIMS}")
print()
print(f"{'tgt':>5} {'q':>2} {'src':>5} {'q':>2} {'rank':>5} {'ker':>4} {'cok':>4}  ker_char")
for r in sorted(ok, key=lambda r: (r["target"], r["source"], r["tq"], r["sq"])):
    print(f"{r['target']:>5} {r['tq']:>2} {r['source']:>5} {r['sq']:>2} "
          f"{r['rank']:>5} {r['dim_ker']:>4} {r['dim_coker']:>4}  {r['ker_char']}")

with open("incidence_operator_enumeration.json", "w") as fh:
    json.dump(dict(irrep_dims=DIMS, records=records), fh, indent=1, sort_keys=True)

print()
print("-" * 74)
print(f"RESULT: {N_PASS} passed, {N_FAIL} failed.")
print("NO TARGET COMPARISON PERFORMED IN THIS FILE.")
print("Table written to incidence_operator_enumeration.json")
