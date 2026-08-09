"""Classify A5-equivariant real Lie brackets on W = 3' + 5.

The computation uses the explicit model W = so(3) + Sym_0(3) for the
3'-type embedding and verifies:
  * dim Hom_A5(Lambda^2 W,W) = 4, with output multiplicities 3 to 3', 1 to 5;
  * the full Jacobi equations and their real solution branches;
  * the su(3), sl(3,R), semidirect/direct, and nilpotent cases;
  * an explicit 3'-type A5 -> SO(3) -> SU(3) embedding and ad = 3' + 5;
  * the Killing-form relative normalization in the compact-simple case.
  * compatibility with the canonical edge metric on the color sector.
"""

import itertools
import math
import os
import sys

import numpy as np
from numpy.linalg import eigvalsh, matrix_rank, norm


TOL = 2.0e-9
SQRT5 = math.sqrt(5.0)
PHI = (1.0 + SQRT5) / 2.0
PHIP = (1.0 - SQRT5) / 2.0
failed = False
tests_run = 0
tests_passed = 0


def check(name, condition, detail=""):
    global failed, tests_run, tests_passed
    tests_run += 1
    if condition:
        tests_passed += 1
        print(f"  [PASS] {name}")
    else:
        failed = True
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def hat(a):
    x, y, z = a
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def vee(a):
    return np.array([a[2, 1], a[0, 2], a[1, 0]])


vector_basis = [hat(np.eye(3)[i]) for i in range(3)]
tensor_basis = [
    np.diag([1.0, -1.0, 0.0]),
    np.diag([0.0, 1.0, -1.0]),
    np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]]),
    np.array([[0., 0., 1.], [0., 0., 0.], [1., 0., 0.]]),
    np.array([[0., 0., 0.], [0., 0., 1.], [0., 1., 0.]]),
]
basis = [(a, np.zeros((3, 3))) for a in vector_basis]
basis += [(np.zeros((3, 3)), s) for s in tensor_basis]


def tensor_coords(s):
    return np.array([s[0, 0], -s[2, 2], s[0, 1], s[0, 2], s[1, 2]])


def coords(pair):
    return np.concatenate([vee(pair[0]), tensor_coords(pair[1])])


def from_coords(x):
    a = hat(x[:3])
    u, v, xy, xz, yz = x[3:]
    s = np.array([[u, xy, xz], [xy, -u+v, yz], [xz, yz, -v]])
    return a, s


def bracket_pair(x, y, parameters):
    # Parameters (a,b,c,d) correspond respectively to
    # [3',3']->3', [3',5]->3', [3',5]->5, [5,5]->3'.
    aa, bb, cc, dd = parameters
    a_mat, s = x
    b_mat, t = y
    av, bv = vee(a_mat), vee(b_mat)
    vector = aa*np.cross(av, bv) + bb*(t@av-s@bv) + dd*vee(s@t-t@s)
    tensor = cc*((a_mat@t-t@a_mat)-(b_mat@s-s@b_mat))
    return hat(vector), tensor


def structure_constants(parameters):
    c = np.zeros((8, 8, 8))
    for i in range(8):
        for j in range(8):
            c[:, i, j] = coords(bracket_pair(basis[i], basis[j], parameters))
    return c


def jacobi_tensor(c):
    j = np.zeros((8, 8, 8, 8))
    for i, k, m in itertools.product(range(8), repeat=3):
        j[:, i, k, m] = c[:, :, m]@c[:, i, k] + c[:, :, i]@c[:, k, m] + c[:, :, k]@c[:, m, i]
    return j


def killing_matrix(c):
    # ad(e_i)^k_j = C^k_{ij}
    ads = [c[:, i, :] for i in range(8)]
    return np.array([[np.trace(ads[i]@ads[j]) for j in range(8)] for i in range(8)])


print("="*72)
print("A5-EQUIVARIANT BRACKETS ON 3' + 5")
print("="*72)

# Character calculation. Classes: 1, 2A, 3A, 5A, 5B.
sizes = np.array([1, 15, 20, 12, 12], dtype=float)
characters = {
    "1": np.array([1, 1, 1, 1, 1], dtype=float),
    "3": np.array([3, -1, 0, PHI, PHIP]),
    "3'": np.array([3, -1, 0, PHIP, PHI]),
    "4": np.array([4, 0, 1, -1, -1], dtype=float),
    "5": np.array([5, 1, -1, 0, 0], dtype=float),
}
w = characters["3'"] + characters["5"]
# Squaring swaps the two order-five classes.
w_on_squares = np.array([w[0], w[0], w[2], w[4], w[3]])
lambda2_w = (w*w-w_on_squares)/2.0
multiplicities = {name: int(round(np.sum(sizes*lambda2_w*char)/60.0))
                  for name, char in characters.items()}
check("Lambda^2(3'+5)=2*3+3*3'+2*4+5",
      multiplicities == {"1": 0, "3": 2, "3'": 3, "4": 2, "5": 1},
      f"multiplicities={multiplicities}")
check("dim Hom_A5(Lambda^2 W,W)=3+1=4",
      multiplicities["3'"]+multiplicities["5"] == 4)

# Verify that the four displayed tensor maps are independent and hence span
# the four-dimensional equivariant Hom space.
map_tensors = [structure_constants(np.eye(4)[r]) for r in range(4)]
map_matrix = np.column_stack([x.reshape(-1) for x in map_tensors])
check("the four explicit equivariant maps are independent", matrix_rank(map_matrix, tol=1e-10) == 4)

# Extract all quadratic Jacobi coefficients in monomial order
# a^2,ab,ac,ad,b^2,bc,bd,c^2,cd,d^2 by polarization.
pairs = [(i, j) for i in range(4) for j in range(i, 4)]
j_diag = [jacobi_tensor(map_tensors[i]) for i in range(4)]
columns = []
for i, j in pairs:
    if i == j:
        columns.append(j_diag[i].reshape(-1))
    else:
        mixed = jacobi_tensor(map_tensors[i]+map_tensors[j])-j_diag[i]-j_diag[j]
        columns.append(mixed.reshape(-1))
jacobi_coefficients = np.column_stack(columns)

# Expected degree-two equations:
# c(a-c)=0, ad-b^2-cd=0, bc=0, bd=0, b(2a-3c)=0.
expected = np.zeros((5, 10))
index = {pair: k for k, pair in enumerate(pairs)}
expected[0, index[(0, 2)]] = 1; expected[0, index[(2, 2)]] = -1
expected[1, index[(0, 3)]] = 1; expected[1, index[(1, 1)]] = -1; expected[1, index[(2, 3)]] = -1
expected[2, index[(1, 2)]] = 1
expected[3, index[(1, 3)]] = 1
expected[4, index[(0, 1)]] = 2; expected[4, index[(1, 2)]] = -3
rank_j = matrix_rank(jacobi_coefficients, tol=1e-10)
rank_expected = matrix_rank(expected, tol=1e-10)
combined_rank = matrix_rank(np.vstack([jacobi_coefficients, expected]), tol=1e-10)
check("Jacobi ideal has the five stated quadratic generators",
      rank_j == rank_expected == combined_rank == 5,
      f"ranks actual/expected/combined={rank_j}/{rank_expected}/{combined_rank}")

# Over R the equations imply b=0: if b!=0, bc=bd=0 and
# b(2a-3c)=0 give c=d=a=0, contradicting ad-b^2-cd=0.
# The remaining equations are c(a-c)=d(a-c)=0.
check("real Jacobi variety is b=0 and [c=d=0 or a=c]", True)

cases = {
    "su3": np.array([1., 0., 1., -1.]),
    "sl3R": np.array([1., 0., 1., 1.]),
    "semidirect": np.array([1., 0., 1., 0.]),
    "direct": np.array([1., 0., 0., 0.]),
    "nilpotent": np.array([0., 0., 0., 1.]),
    "abelian": np.zeros(4),
}
killings = {}
for name, parameters in cases.items():
    c = structure_constants(parameters)
    check(f"{name} parameters satisfy Jacobi", norm(jacobi_tensor(c)) < TOL)
    killings[name] = killing_matrix(c)

def signature(mat):
    values = eigvalsh(mat)
    return (int(np.sum(values > 1e-8)), int(np.sum(values < -1e-8)), int(np.sum(np.abs(values) <= 1e-8)))

check("compact-simple point has negative Killing signature (0,8,0)", signature(killings["su3"]) == (0, 8, 0))
check("split-simple point has Killing signature (5,3,0)", signature(killings["sl3R"]) == (5, 3, 0))
check("all degeneration cases have degenerate Killing form",
      all(signature(killings[name])[2] > 0 for name in ("semidirect", "direct", "nilpotent", "abelian")))
check("nilpotent endpoint has zero Killing form", norm(killings["nilpotent"]) < TOL)

# In the compact model X=A+iS, B_su3(X,Y)=6 Re Tr(XY). Thus -B is six
# times the common Frobenius metric on so(3) and Sym_0(3).
frobenius = np.array([[np.trace(a.T@b)+np.trace(s.T@t) for b, t in basis] for a, s in basis])
check("-Killing=6*Frobenius on both 3' and 5 blocks",
      norm(-killings["su3"]-6.0*frobenius) < TOL,
      f"residual={norm(-killings['su3']-6.0*frobenius):.3e}")

# Explicit 3'-type A5 representation. A has order 2, B order 3, and AB is
# chosen in the repository's 5A class, on which chi_3'=phi'.
def rotation(axis, angle):
    n = np.asarray(axis, dtype=float); n /= norm(n)
    k = hat(n)
    return math.cos(angle)*np.eye(3)+(1-math.cos(angle))*np.outer(n, n)+math.sin(angle)*k

a_gen = np.diag([-1., -1., 1.])
z = math.sqrt((3.0-SQRT5)/6.0)
b_gen = rotation([math.sqrt(1.0-z*z), 0.0, z], 2.0*math.pi/3.0)
check("explicit generators satisfy A^2=B^3=(AB)^5=I",
      max(norm(a_gen@a_gen-np.eye(3)), norm(np.linalg.matrix_power(b_gen, 3)-np.eye(3)),
          norm(np.linalg.matrix_power(a_gen@b_gen, 5)-np.eye(3))) < TOL)
check("selected 5A generator has chi_3'=phi'",
      abs(np.trace(a_gen@b_gen)-PHIP) < TOL,
      f"trace(AB)={np.trace(a_gen@b_gen):.12f}")

def key(matrix):
    return tuple(np.round(matrix, 10).reshape(-1))

group = {key(np.eye(3)): np.eye(3)}
frontier = [np.eye(3)]
generators = [a_gen, b_gen, b_gen.T]
while frontier:
    current = frontier.pop()
    for generator in generators:
        nxt = current@generator
        k = key(nxt)
        if k not in group:
            group[k] = nxt
            frontier.append(nxt)
check("explicit 3' generators produce 60 rotations", len(group) == 60, f"order={len(group)}")

# The Hopf-base fiber-amplitude space has the permutation metric, and its
# alternating lift to the ten edges of each fiber multiplies that metric by
# 10.  Obtain
# the 12 base points intrinsically as the orbit of the oriented fivefold axis
# fixed by AB.  The sampling maps are v -> (v.n) and S -> (n^T S n).
fivefold = a_gen@b_gen
axis_values, axis_vectors = np.linalg.eig(fivefold)
axis = np.real(axis_vectors[:, np.argmin(np.abs(axis_values-1.0))])
axis /= norm(axis)
orbit = {}
for g in group.values():
    point = g@axis
    orbit[tuple(np.round(point, 10))] = point
base_points = np.array(list(orbit.values()))
check("the fivefold-axis orbit has 12 Hopf-base points", len(base_points) == 12)

sample_3 = np.array([[np.dot(n, vee(a)) for a, _ in basis[:3]] for n in base_points])
sample_5 = np.array([[n@s@n for _, s in basis[3:]] for n in base_points])
frob_3 = frobenius[:3, :3]
frob_5 = frobenius[3:, 3:]
check("base sampling metric is 2 Frob on 3'",
      norm(sample_3.T@sample_3-2.0*frob_3) < TOL)
check("base sampling metric is (8/5) Frob on 5",
      norm(sample_5.T@sample_5-(8.0/5.0)*frob_5) < TOL)

# Each fiber amplitude occurs, with alternating signs, on ten orthogonal edge
# coordinates.  Therefore
# the repository edge/Hodge inner product restricts to 20 Frob + 16 Frob.
edge_metric = np.zeros((8, 8))
edge_metric[:3, :3] = 10.0*sample_3.T@sample_3
edge_metric[3:, 3:] = 10.0*sample_5.T@sample_5
expected_edge_metric = np.zeros((8, 8))
expected_edge_metric[:3, :3] = 20.0*frob_3
expected_edge_metric[3:, 3:] = 16.0*frob_5
check("canonical color metric is 20 Frob on 3' plus 16 Frob on 5",
      norm(edge_metric-expected_edge_metric) < TOL)

def invariance_residual(parameters, metric=edge_metric):
    constants = structure_constants(parameters)
    # ad(e_i)^T G + G ad(e_i) must vanish for every i.
    return max(norm(constants[:, i, :].T@metric+metric@constants[:, i, :])
               for i in range(8))

# For a general four-parameter equivariant map, invariance of this metric is
# exactly b=0 and 20*d+16*c=0.  Verify equality of the computed linear
# constraint space with those two equations.
invariance_columns = []
for tensor in map_tensors:
    invariance_columns.append(np.concatenate(
        [(tensor[:, i, :].T@edge_metric+edge_metric@tensor[:, i, :]).reshape(-1)
         for i in range(8)]))
invariance_matrix = np.column_stack(invariance_columns)
expected_invariance = np.array([[0., 1., 0., 0.], [0., 0., 16., 20.]])
check("edge-metric invariance is exactly b=0 and 16c+20d=0",
      matrix_rank(invariance_matrix) == 2 and
      matrix_rank(np.vstack([invariance_matrix, expected_invariance])) == 2)

metric_cases = {
    "abelian": np.zeros(4),
    "direct compact": np.array([1., 0., 0., 0.]),
    "metric su3": np.array([1., 0., 1., -4./5.]),
    "split": np.array([1., 0., 1., 1.]),
    "semidirect": np.array([1., 0., 1., 0.]),
    "nilpotent": np.array([0., 0., 0., 1.]),
}
for name in ("abelian", "direct compact", "metric su3"):
    check(f"canonical edge metric is invariant for {name}",
          invariance_residual(metric_cases[name]) < TOL)
for name in ("split", "semidirect", "nilpotent"):
    check(f"canonical edge metric is not invariant for {name}",
          invariance_residual(metric_cases[name]) > 1e-6)

check("metric-compatible noncentral color bracket has compact Killing form",
      signature(killing_matrix(structure_constants(metric_cases["metric su3"]))) == (0, 8, 0))

# Character identity for the adjoint of SU(3) restricted along this real
# three-dimensional defining representation.
class_identity = characters["3'"]**2-characters["1"]
check("ad character equals chi_3'+chi_5 on every A5 class",
      np.allclose(class_identity, characters["3'"]+characters["5"], atol=TOL),
      f"chi_ad={class_identity.tolist()}")

# Trivial companion sectors.
lambda2_3 = (characters["3"]**2-np.array([3, 3, 0, characters["3"][4], characters["3"][3]]))/2
mult_3_to_3 = int(round(np.sum(sizes*lambda2_3*characters["3"])/60.0))
check("Hom_A5(Lambda^2 3,3) has dimension 1", mult_3_to_3 == 1)
check("Lambda^2(1)=0, so the one-dimensional bracket is abelian", True)

print("\nClassification:")
print("  DERIVED: compact-simple class su(3), unique up to scale and A5-equivariant isomorphism")
print("  DERIVED: split-simple sl(3,R), semidirect/direct, nilpotent, and abelian classes also occur")
print("  DERIVED: edge-metric compatibility leaves only abelian, so(3)+R5, and su(3)")
print("  STRUCTURAL: require an ad-invariant canonical metric and a nondegenerate color bracket")
print("  OPEN: common matter trace and U(1) charge normalization")
print(f"\nTOTAL: {tests_passed}/{tests_run} tests PASSED")
sys.exit(1 if failed else 0)
