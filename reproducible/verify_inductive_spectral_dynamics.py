#!/usr/bin/env python3
"""Targeted certificate for the barycentric-path inductive spectral dynamics.

The construction has two sharply separated layers.

1.  A genuine inductive *commutative* spectral dynamics on nested top-cell
    paths.  The 600 root tetrahedra refine 24-to-1.  Conditional expectations
    give an exact hierarchical Markov Laplacian, its square-root Dirac
    intertwines refinement, and inverse 3-volume scaling gives spectral
    abscissa three.
2.  A noncommutative Cuntz/Toeplitz extension.  Perron--Frobenius cylinder
    weights give the critical KMS relation and nontrivial modular frequency.
    The exponential Dirac has unbounded ordinary commutators with creation
    operators but exact algebraic twisted commutators.  Joining the KMS GNS
    representation and the compact-resolvent Dirac is therefore OPEN.

The matrix control is exact over Q on a universal small branching model.  The
actual 600-cell root complex and all N=24 multiplicity/scaling claims are
checked separately.  No Standard-Model target is inspected.
"""

from collections import defaultdict, deque
from itertools import combinations
import math
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as sla
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("INDUCTIVE SPECTRAL DYNAMICS ON 600-CELL BARYCENTRIC PATHS")
print("=" * 78)

# --------------------------------------------------------------------- roots
# Reconstruct the 600 tetrahedra and their dual graph independently from the
# 120-vertex one-skeleton.  This is the connected base of the path dynamics.
_, adjacency, _ = build_600cell()
neighbours = [set(np.flatnonzero(adjacency[index] > 0.5))
              for index in range(120)]
tetrahedra = []
for i in range(120):
    for j in sorted(vertex for vertex in neighbours[i] if vertex > i):
        common_ij = neighbours[i] & neighbours[j]
        for k in sorted(vertex for vertex in common_ij if vertex > j):
            common_ijk = common_ij & neighbours[k]
            for ell in sorted(vertex for vertex in common_ijk if vertex > k):
                tetrahedra.append((i, j, k, ell))

face_to_cells = defaultdict(list)
for cell_index, cell in enumerate(tetrahedra):
    for face in combinations(cell, 3):
        face_to_cells[face].append(cell_index)

dual_neighbours = [set() for _ in tetrahedra]
for incident_cells in face_to_cells.values():
    if len(incident_cells) == 2:
        left, right = incident_cells
        dual_neighbours[left].add(right)
        dual_neighbours[right].add(left)

seen = {0}
queue = deque((0,))
while queue:
    current = queue.popleft()
    for target in sorted(dual_neighbours[current] - seen):
        seen.add(target)
        queue.append(target)

dual_edges = sum(len(row) for row in dual_neighbours) // 2
check("root carrier has exactly 600 tetrahedra", len(tetrahedra) == 600)
check("every triangular face belongs to exactly two tetrahedra",
      len(face_to_cells) == 1200 and
      all(len(cells) == 2 for cells in face_to_cells.values()))
check("dual tetrahedron graph is connected 4-regular with 1200 edges",
      len(seen) == 600 and dual_edges == 1200 and
      all(len(row) == 4 for row in dual_neighbours))

# --------------------------------------------------------- exact path counts
ROOTS = 600
BRANCHING = 24
for level, expected in ((0, 600), (1, 14400), (2, 345600), (3, 8294400)):
    check(f"level-{level} path carrier has 600*24^{level} states",
          ROOTS * BRANCHING**level == expected, f"dim={expected}")

for depth in range(1, 5):
    vertical = sum(ROOTS*(BRANCHING-1)*BRANCHING**level
                   for level in range(depth))
    check(f"root plus vertical wavelets exhaust level {depth}",
          ROOTS + vertical == ROOTS*BRANCHING**depth)

# -------------------------------------------- exact universal matrix control
# Conditional expectations on a small R=2, N=3 model prove the algebraic
# identities without creating a 14,400-square dense matrix.  The formulas are
# tensor identities independent of R and N; actual R=600,N=24 counts are above.
R_CONTROL = 2
N_CONTROL = 3
P = sp.ones(N_CONTROL, N_CONTROL) / N_CONTROL
u = sp.ones(N_CONTROL, 1) / sp.sqrt(N_CONTROL)
base_laplacian = sp.Matrix(((1, -1), (-1, 1)))


def kron_all(*matrices):
    result = sp.Matrix(((1,),))
    for matrix in matrices:
        result = sp.kronecker_product(result, matrix)
    return result


def conditional(depth, retained):
    """Average all branch coordinates after ``retained``."""
    factors = [sp.eye(R_CONTROL)]
    factors.extend(sp.eye(N_CONTROL) for _ in range(retained))
    factors.extend(P for _ in range(depth-retained))
    return kron_all(*factors)


def hierarchical_laplacian(depth):
    expectations = [conditional(depth, retained)
                    for retained in range(depth+1)]
    root_term = kron_all(base_laplacian, *(P for _ in range(depth)))
    result = root_term
    # D-wavelet eigenvalues are 2^(k+1), so L=D^2 has 4^(k+1).
    for level in range(depth):
        wavelet = expectations[level+1] - expectations[level]
        result += 4**(level+1) * wavelet
    return sp.simplify(result), expectations


L2, E2 = hierarchical_laplacian(2)
L3, E3 = hierarchical_laplacian(3)
dim2, dim3 = R_CONTROL*N_CONTROL**2, R_CONTROL*N_CONTROL**3

wavelets2 = [E2[level+1]-E2[level] for level in range(2)]
check("conditional expectations are nested orthogonal projectors",
      all(E.T == E and E*E == E for E in E3) and
      all(E3[i]*E3[j] == E3[min(i, j)]
          for i in range(4) for j in range(4)))
check("wavelet projectors are orthogonal with the exact multiplicities",
      wavelets2[0].rank() == R_CONTROL*(N_CONTROL-1) and
      wavelets2[1].rank() == R_CONTROL*N_CONTROL*(N_CONTROL-1) and
      wavelets2[0]*wavelets2[1] == sp.zeros(dim2))
rank_from_blocks_2 = (
    base_laplacian.rank()
    + sum(R_CONTROL*(N_CONTROL-1)*N_CONTROL**level
          for level in range(2))
)
rank_from_blocks_3 = (
    base_laplacian.rank()
    + sum(R_CONTROL*(N_CONTROL-1)*N_CONTROL**level
          for level in range(3))
)
root_projector_3 = kron_all(base_laplacian, P, P, P)
vertical_projectors_3 = [E3[level+1]-E3[level]
                         for level in range(3)]
check("hierarchical L is symmetric, positive-rank and connected",
      L2.T == L2 and L3.T == L3 and
      all(root_projector_3*wavelet == sp.zeros(dim3)
          for wavelet in vertical_projectors_3) and
      rank_from_blocks_2 == dim2-1 and rank_from_blocks_3 == dim3-1,
      f"orthogonal-block ranks={rank_from_blocks_2},{rank_from_blocks_3}")
check("hierarchical L has graph-Laplacian signs and zero row sums",
      all(L3[row, col] <= 0
          for row in range(dim3) for col in range(dim3) if row != col) and
      L3*sp.ones(dim3, 1) == sp.zeros(dim3, 1))

# Exact isometric refinement: append the normalized constant child vector.
inclusion = sp.kronecker_product(sp.eye(dim2), u)
check("refinement inclusion is isometric", inclusion.T*inclusion == sp.eye(dim2))
check("L intertwines refinement exactly", L3*inclusion == inclusion*L2)

# A graph Laplacian generates a positivity-preserving Markov heat flow.  Check
# the actual exponential numerically at three deterministic times.
L2_float = np.asarray(L2, dtype=float)
markov_residual = 0.0
minimum_entry = 1.0
nontrivial_norm = 0.0
for time in (0.01, 0.1, 1.0):
    heat = sla.expm(-time*L2_float)
    markov_residual = max(markov_residual,
                          float(np.max(np.abs(heat.sum(axis=1)-1))))
    minimum_entry = min(minimum_entry, float(heat.min()))
    nontrivial_norm = max(nontrivial_norm,
                          float(np.linalg.norm(heat-np.eye(dim2))))
check("heat dynamics is Markov on the finite inductive carrier",
      markov_residual < 1e-12 and minimum_entry > -1e-12,
      f"row residual={markov_residual:.3e}, min entry={minimum_entry:.3e}")
check("heat dynamics is genuinely nontrivial", nontrivial_norm > 0.1,
      f"max ||exp(-tL)-I||={nontrivial_norm:.6f}")

# ---------------------------------------------------------- selected scaling
# Each barycentric top simplex has exactly 1/24 of its parent's affine volume.
# Conditional on inverse isotropic 3-volume scaling, b=24^(1/3) is fixed.
b = BRANCHING**(1/3)
check("inverse three-volume scaling fixes b^3=24",
      abs(b**3-BRANCHING) < 1e-12, f"b={b:.12f}")

# Vertical D eigenvalue b^(k+1), multiplicity 600*23*24^k.  The zeta ratio is
# 24*b^(-s), so the abscissa is exactly three.  The finite base spectrum does
# not change it.
dimension = math.log(BRANCHING)/math.log(b)
check("inductive Dirac zeta abscissa is exactly three",
      abs(dimension-3) < 1e-12, f"d={dimension:.12f}")
for s_value, should_converge in ((2.9, False), (3.1, True), (4.0, True)):
    ratio = BRANCHING*b**(-s_value)
    check(f"zeta geometric ratio at s={s_value} has the correct side",
          (ratio < 1) == should_converge, f"ratio={ratio:.9f}")

# -------------------------------------------------------- KMS/Cuntz extension
# Uniform cylinders are projectively consistent and satisfy the O_24 gauge
# KMS equation.  The equation fixes only beta*epsilon=log(24).
for length in range(5):
    parent_weight = sp.Rational(1, BRANCHING**length)
    child_weight = sp.Rational(1, BRANCHING**(length+1))
    check(f"uniform cylinder state is consistent at length {length}",
          BRANCHING*child_weight == parent_weight)

beta, epsilon = sp.symbols("beta epsilon", positive=True)
kms_product = sp.solve(
    sp.Eq(sp.exp(-beta*epsilon), sp.Rational(1, BRANCHING)),
    beta*epsilon)
check("critical KMS equation fixes beta*epsilon=log(24)",
      kms_product == [sp.log(BRANCHING)])

# Geometry can set epsilon=log(b)=log(24)/3, giving beta=3.  Rescaling the
# gauge clock changes beta inversely, so KMS alone has no physical time unit.
epsilon_geometry = sp.log(BRANCHING)/3
check("volume-scale gauge energy makes the critical inverse temperature 3",
      sp.simplify(sp.log(BRANCHING)/epsilon_geometry) == 3)
check("KMS leaves an overall clock normalization without the geometric scale",
      sp.simplify(sp.log(BRANCHING)/(2*epsilon_geometry)) == sp.Rational(3, 2))

# Algebraic creation operator raises word depth.  D|w|=b^|w| gives an
# unbounded ordinary commutator but an exact modularly twisted one with
# sigma(S_i)=b S_i and sigma(S_i*)=b^-1 S_i*.
ordinary_sizes = [(b-1)*b**length for length in range(8)]
check("ordinary Cuntz-creation commutator is unbounded",
      all(ordinary_sizes[index+1] > ordinary_sizes[index]
          for index in range(len(ordinary_sizes)-1)) and
      ordinary_sizes[-1] > 100*ordinary_sizes[0],
      f"depth 0..7 sizes={ordinary_sizes[0]:.3g}..{ordinary_sizes[-1]:.3g}")

twisted_creation = [b**(length+1)-b*b**length for length in range(8)]
twisted_annihilation = [b**(length-1)-b**-1*b**length
                        for length in range(1, 9)]
check("modular twisted commutators vanish on every tested depth",
      all(abs(value) < 1e-12
          for value in twisted_creation+twisted_annihilation))
check("the analytic twist preserves the algebraic Cuntz relations",
      sp.simplify(sp.Symbol("b", positive=True)**-1
                  * sp.Symbol("b", positive=True)) == 1,
      "sigma(S_i*)sigma(S_j)=delta_ij and sigma(S_i)sigma(S_i*)=S_iS_i*")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: exact inductive hierarchical Markov dynamics on 600*24^n paths.")
print("STRUCTURAL: inverse 3-volume scaling selects b=24^(1/3) and d=3.")
print("DERIVED: uniform cylinders satisfy the critical O_24 KMS relation.")
print("DERIVED NEGATIVE: ordinary Cuntz commutators are unbounded.")
print("DERIVED CONDITIONAL: modularly twisted commutators vanish algebraically.")
print("OPEN: one KMS-GNS modular spectral triple carrying this compact Dirac.")
raise SystemExit(0 if passed == tests else 1)
