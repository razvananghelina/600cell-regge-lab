#!/usr/bin/env python3
"""Finite screens for the proposed mathematics-to-physics bridge.

Exact parts:
  * four-summand Krajewski design census on one 60-dimensional chamber sheet;
  * the C5-coset 12+48 chamber decomposition and its character obstruction;
  * induction from C10 versus all (C10,C4,C6) orbifold stabilizers;
  * the cubic-Cayley obstruction to a free 2I chamber symmetry.

The final inverse-spectral calculation is explicitly a numerical detector.
It is not used as an exact existence certificate.
"""

import contextlib
from fractions import Fraction
import io
import itertools
import runpy

import networkx as nx
import numpy as np
import sympy as sy


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("MATHEMATICS-TO-PHYSICS BRIDGE AUDIT")
print("=" * 78)

# ---------------------------------------------------------------------------
# 1. Necessary Krajewski-design screen for the four physical algebra types.
# ---------------------------------------------------------------------------

# Complex carrier sizes for C, C, H, M3(C).  H acts on C^2 in its defining
# representation.  A signed multiplicity x_ij chooses one orientation of an
# off-diagonal Krajewski cell; absence of diagonal and reverse-paired cells is
# the metric-dimension-zero orientability screen.
sizes = (1, 1, 2, 3)
pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
weights = tuple(sizes[i] * sizes[j] for i, j in pairs)


def cells_from_signed(signed):
    cells = []
    for (i, j), value, weight in zip(pairs, signed, weights):
        if value > 0:
            cells.append((i, j, value * weight))
        elif value < 0:
            cells.append((j, i, -value * weight))
    return cells


def structural_rank(signed):
    """Capacitated generic rank of first-order-allowed Dirac blocks."""
    cells = cells_from_signed(signed)
    flow = nx.DiGraph()
    source, sink = "source", "sink"
    for q, (_, _, dimension) in enumerate(cells):
        flow.add_edge(source, ("row", q), capacity=dimension)
        flow.add_edge(("column", q), sink, capacity=dimension)
    for q, (i, j, dim_q) in enumerate(cells):
        for r, (k, l, dim_r) in enumerate(cells):
            if q != r and (k == j or l == i):
                flow.add_edge(("row", q), ("column", r),
                              capacity=min(dim_q, dim_r))
    return nx.maximum_flow_value(flow, source, sink)


dimension_solutions = 0
unimodular_signed = 0
full_rank_signed = 0
first_full = None

# Exhaust every nonnegative magnitude solution of
#   a + 2b + 3c + 2d + 3e + 6f = 60,
# then every genuine sign choice on its nonzero entries.  For a 4x4
# antisymmetric cap matrix the Pfaffian is af-be+cd.
for a in range(61):
    remainder1 = 60 - a
    for b in range(remainder1 // 2 + 1):
        remainder2 = remainder1 - 2*b
        for c in range(remainder2 // 3 + 1):
            remainder3 = remainder2 - 3*c
            for d in range(remainder3 // 2 + 1):
                remainder4 = remainder3 - 2*d
                for e in range(remainder4 // 3 + 1):
                    remainder5 = remainder4 - 3*e
                    if remainder5 % 6:
                        continue
                    f = remainder5 // 6
                    magnitudes = (a, b, c, d, e, f)
                    dimension_solutions += 1
                    sign_sets = [(-1, 1) if value else (1,)
                                 for value in magnitudes]
                    for signs in itertools.product(*sign_sets):
                        signed = tuple(value*sign for value, sign
                                       in zip(magnitudes, signs))
                        pfaffian = (signed[0]*signed[5]
                                     - signed[1]*signed[4]
                                     + signed[2]*signed[3])
                        if abs(pfaffian) != 1:
                            continue
                        unimodular_signed += 1
                        if structural_rank(signed) == 60:
                            full_rank_signed += 1
                            if first_full is None:
                                first_full = signed

check("[DERIVED] weighted dimension equation has 57,563 magnitude solutions",
      dimension_solutions == 57563)
check("[DERIVED] 23,584 oriented signed designs are integrally unimodular",
      unimodular_signed == 23584,
      "raw signed count; sheet reversal is not quotiented")
check("[DERIVED] 3,592 designs also have full first-order structural rank",
      full_rank_signed == 3592)

EXPLICIT_DESIGN = (0, -1, 1, -8, 7, -3)
cap = sy.zeros(4)
for (i, j), value in zip(pairs, EXPLICIT_DESIGN):
    cap[i, j] = value
    cap[j, i] = -value
physical_k0_scale = sy.diag(1, 1, 2, 1)
physical_cap = physical_k0_scale * cap * physical_k0_scale
check("[DERIVED] explicit C+C+M2+M3 design closes every cheap gate",
      first_full == EXPLICIT_DESIGN
      and sum(abs(x)*w for x, w in zip(EXPLICIT_DESIGN, weights)) == 60
      and cap.det() == 1
      and structural_rank(EXPLICIT_DESIGN) == 60,
      f"signed multiplicities={EXPLICIT_DESIGN}, det(Cap_C)={cap.det()}")
check("[DERIVED] quaternionic K0 normalization remains nondegenerate",
      physical_cap.det() == 4,
      "H has no complex rank-one projection; nondegeneracy survives")

# ---------------------------------------------------------------------------
# 2. The tempting 60=12+48 chamber decomposition, and why it is not matter.
# ---------------------------------------------------------------------------

with contextlib.redirect_stdout(io.StringIO()):
    chamber = runpy.run_path("verify_oriented_chamber_double.py")

gamma = chamber["gamma"]
plus = np.flatnonzero(gamma == 1)
plus_position = {int(x): i for i, x in enumerate(plus)}
rotations = chamber["chamber_rotations"]


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


identity120 = tuple(range(120))


def permutation_order(p):
    power = identity120
    for order in range(1, 7):
        power = compose(p, power)
        if power == identity120:
            return order
    raise AssertionError("unexpected A5 order")


g5 = next(p for p in rotations if permutation_order(p) == 5)
R5 = np.zeros((60, 60), dtype=np.int64)
for chamber_index in plus:
    R5[plus_position[int(g5[int(chamber_index)])],
       plus_position[int(chamber_index)]] = 1
Pnum = np.zeros((60, 60), dtype=np.int64)
power = np.eye(60, dtype=np.int64)
for _ in range(5):
    Pnum += power
    power = power @ R5

D_fixed = chamber["D"].toarray().astype(np.int64)
J_fixed = chamber["J"].toarray().astype(np.int64)
S_fixed = (D_fixed @ J_fixed)[np.ix_(plus, plus)]
check("[DERIVED] C5 averaging cuts one chamber sheet as 12+48",
      np.array_equal(Pnum @ Pnum, 5*Pnum)
      and int(np.trace(Pnum)) == 60,
      "rank(P_C5)=Tr(P_C5)=12; complement rank=48=3*16")
check("[DERIVED] the 12 and 48 sectors are preserved by fixed chamber D,J,gamma",
      np.array_equal(S_fixed @ Pnum, Pnum @ S_fixed),
      "equivalently [S,P_C5]=0 on H+")

# Exact 2I character data and all Hopf harmonic inductions.
with contextlib.redirect_stdout(io.StringIO()):
    harmonic = runpy.run_path("verify_fibonacci_nonbinary_dynamics.py")

irreps = harmonic["irreps_2i"]
irrep_dims = harmonic["irrep_dims"]
induced_c10 = tuple(tuple(int(value) for value in row)
                    for row in harmonic["all_induced"])

# In the verifier's irrep order, diagonal M16 restricts as
# 2*rho1 + rho2 + 2*rho5 + rho9.
m16 = (2, 1, 0, 0, 2, 0, 0, 0, 1)
three_m16 = tuple(3*x for x in m16)
center_trace_m16 = sum(mult*int(character[1])
                       for mult, character in zip(m16, irreps))
check("[DERIVED NEGATIVE] scalar 48 complement is not three diagonal M16s",
      48 != 3*center_trace_m16,
      "central traces: scalar complement +48, three M16 modules 0")

# Four C10-induced line modules are the only direct sums of whole Hopf-orbit
# line modules with total dimension 48.  Exhaust them with repetition.
four_harmonic_hits = []
for harmonics in itertools.combinations_with_replacement(range(10), 4):
    total = tuple(sum(induced_c10[q][i] for q in harmonics)
                  for i in range(9))
    if total == three_m16:
        four_harmonic_hits.append(harmonics)
c10_matrix = sy.Matrix(induced_c10)
augmented_c10 = c10_matrix.col_join(sy.Matrix([m16]))
check("[DERIVED NEGATIVE] no four C10 harmonics reproduce three M16s",
      not four_harmonic_hits)
check("[DERIVED NEGATIVE] M16 is not even in the virtual C10 induction span",
      c10_matrix.rank() == 6 and augmented_c10.rank() == 7)

# Add every character of the other two derived stabilizers.  Power-class
# maps use the same class order as the exact 2I table:
# 1A,2A,4A,6A,3A,10A,5A,5B,10B.
power_classes = {
    4: (0, 2, 1, 2),
    6: (0, 3, 4, 1, 4, 3),
}
all_rows = list(induced_c10)
row_names = [f"C10:{q}" for q in range(10)]
for order, classes in power_classes.items():
    root = sy.exp(2*sy.pi*sy.I/order)
    for q in range(order):
        row = []
        for character in irreps:
            inner = sum(character[class_index] * root**(-q*k)
                        for k, class_index in enumerate(classes)) / order
            value = sy.simplify(sy.expand_complex(inner))
            assert value.is_Integer
            row.append(int(value))
        all_rows.append(tuple(row))
        row_names.append(f"C{order}:{q}")

all_induction_matrix = sy.Matrix(all_rows)
check("[DERIVED] C10+C4+C6 inductions span the full 2I representation ring",
      all_induction_matrix.rank() == 9)

# Sparse exact virtual identity found by the target-driven search:
# M16 = Ind_C10 chi0 - 2 Ind_C10 chi2 - Ind_C10 chi3
#       + Ind_C6 chi0 + Ind_C6 chi1.
coefficients = [0]*20
coefficients[0] = 1
coefficients[2] = -2
coefficients[3] = -1
coefficients[14] = 1
coefficients[15] = 1
virtual_sum = tuple(sum(coefficients[r]*all_rows[r][i]
                        for r in range(20)) for i in range(9))
check("[DERIVED identity; STRUCTURAL interpretation] an orbifold virtual M16 exists",
      virtual_sum == m16,
      "C10:0 - 2 C10:2 - C10:3 + C6:0 + C6:1")

# ---------------------------------------------------------------------------
# 3. A relabelling of chambers by 2I cannot repair the missing spin lift.
# ---------------------------------------------------------------------------

with contextlib.redirect_stdout(io.StringIO()):
    binary = runpy.run_path("verify_nonnormal_c10_selection.py")
orders = [binary["element_order"](g) for g in range(len(binary["group"]))]
unique_involution = orders.count(2) == 1
noncommutative = any(
    binary["mul"][a][b] != binary["mul"][b][a]
    for a in range(120) for b in range(120)
)
check("[DERIVED NEGATIVE] no connected cubic Cayley graph exists on 2I",
      unique_involution and noncommutative,
      "an inverse-closed 3-set is {-1,g,g^-1}, which generates an abelian subgroup")

# ---------------------------------------------------------------------------
# 4. Numerical detector: fixed-D inverse spectral compatibility.
# ---------------------------------------------------------------------------


def allowed_mask(signed):
    cells = cells_from_signed(signed)
    starts = np.cumsum([0] + [dimension for _, _, dimension in cells])
    mask = np.zeros((60, 60), dtype=bool)
    for q, (i, j, _) in enumerate(cells):
        for r, (k, l, _) in enumerate(cells):
            if q != r and (k == j or l == i):
                mask[starts[q]:starts[q+1], starts[r]:starts[r+1]] = True
    return mask, cells


mask, numerical_cells = allowed_mask(EXPLICIT_DESIGN)
target_eigenvalues = np.sort(np.linalg.eigvalsh(S_fixed.astype(float)))
rng = np.random.default_rng(0)
trial = rng.normal(size=(60, 60))
trial = (trial + trial.T)/2
trial *= mask
best_residual = float("inf")
best_orbit_point = None
for _ in range(120):
    _, eigenvectors = np.linalg.eigh(trial)
    orbit_point = (eigenvectors*target_eigenvalues) @ eigenvectors.T
    residual = np.linalg.norm(orbit_point[~mask]) / np.linalg.norm(orbit_point)
    if residual < best_residual:
        best_residual = residual
        best_orbit_point = orbit_point.copy()
    trial = orbit_point * mask

eigenvalue_residual = np.max(np.abs(
    np.linalg.eigvalsh(best_orbit_point) - target_eigenvalues
))
check("[PATTERN] real isospectral first-order detector reaches machine precision",
      best_residual < 1e-9 and eigenvalue_residual < 1e-11,
      f"forbidden relative residual={best_residual:.3e}")

# Connectedness detector.  In the J-paired basis, on a cell (i,j), A acts
# from factor i on H+ and from factor j on H-.  Full complex matrix units are
# a stronger screen than restricting the M2 factor to its quaternionic form.
starts = np.cumsum([0] + [dimension for _, _, dimension in numerical_cells])
plus_actions = []
minus_actions = []
for factor, size in enumerate(sizes):
    for a in range(size):
        for b in range(size):
            plus_action = np.zeros((60, 60))
            minus_action = np.zeros((60, 60))
            unit = np.zeros((size, size))
            unit[a, b] = 1
            for q, (i, j, dimension) in enumerate(numerical_cells):
                multiplicity = dimension // (sizes[i]*sizes[j])
                block = slice(starts[q], starts[q+1])
                if i == factor:
                    plus_action[block, block] = np.kron(
                        np.kron(unit, np.eye(sizes[j])), np.eye(multiplicity)
                    )
                if j == factor:
                    minus_action[block, block] = np.kron(
                        np.kron(np.eye(sizes[i]), unit), np.eye(multiplicity)
                    )
            plus_actions.append(plus_action)
            minus_actions.append(minus_action)

commutator_map = np.column_stack([
    (best_orbit_point*0 + best_orbit_point @ minus_action
     - plus_action @ best_orbit_point).reshape(-1)
    for plus_action, minus_action in zip(plus_actions, minus_actions)
])
commutator_singulars = np.linalg.svd(commutator_map, compute_uv=False)
numerical_kernel = len(commutator_singulars) - np.count_nonzero(
    commutator_singulars > 1e-8
)
check("[PATTERN] numerical fixed-D candidate is connected and fluctuating",
      numerical_kernel == 1 and commutator_singulars[0] > 1,
      (f"commutant kernel={numerical_kernel}; "
       f"smallest nonzero singular={commutator_singulars[-2]:.3e}"))

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT_ALGEBRA_TYPE=NOT_OBSTRUCTED_STRUCTURALLY")
print("VERDICT_48_SCALAR_COMPLEMENT=REFUTED_AS_THREE_M16")
print("VERDICT_C10_ONLY_MATTER_FUNCTOR=REFUTED")
print("VERDICT_ORBIFOLD_INDEX_ROUTE=OPEN")
print("VERDICT_FIXED_D_PHYSICAL_TYPE=NUMERICAL_PATTERN_ONLY")
if passed != tests:
    raise SystemExit(1)

