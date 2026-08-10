#!/usr/bin/env python3
"""Exact all-degree Whitney/Kaehler--Dirac barycentric induction certificate.

The carrier is one regular Euclidean tetrahedron and its full barycentric
subdivision.  Local affine naturality and assembly make every exact identity
valid tetrahedron-by-tetrahedron on a piecewise-flat simplicial 3-manifold,
including the boundary complex of the 600-cell.

The certificate distinguishes exact Galerkin compression from the stronger,
false operator-intertwining claim.  No refinement weights are fitted.
"""

from itertools import combinations, permutations
from math import factorial

import numpy as np
import scipy.linalg as sla
import sympy as sy


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
print("ALL-DEGREE WHITNEY / KAEHLER--DIRAC INDUCTIVE CERTIFICATE")
print("=" * 78)

# An unscaled regular tetrahedron keeps every affine integral rational.
reference_vertices = {
    0: sy.Matrix((1, 1, 1)),
    1: sy.Matrix((1, -1, -1)),
    2: sy.Matrix((-1, 1, -1)),
    3: sy.Matrix((-1, -1, 1)),
}
coarse_top = [(0, 1, 2, 3)]


def subset_members(mask):
    return tuple(index for index in range(4) if mask & (1 << index))


fine_barycentric = {
    mask: sy.Matrix(tuple(
        sy.Rational(int(index in subset_members(mask)),
                    len(subset_members(mask)))
        for index in range(4)
    ))
    for mask in range(1, 16)
}
fine_vertices = {
    mask: sum((fine_barycentric[mask][index]*reference_vertices[index]
               for index in range(4)), sy.zeros(3, 1))
    for mask in range(1, 16)
}
fine_top = []
for ordering in permutations(range(4)):
    cumulative = 0
    flag = []
    for vertex in ordering:
        cumulative |= 1 << vertex
        flag.append(cumulative)
    fine_top.append(tuple(flag))


def all_simplices(top_cells):
    return [
        sorted({tuple(face)
                for top in top_cells
                for face in combinations(top, degree+1)})
        for degree in range(4)
    ]


coarse = all_simplices(coarse_top)
fine = all_simplices(fine_top)
coarse_dims = tuple(map(len, coarse))
fine_dims = tuple(map(len, fine))
check("coarse tetrahedron f-vector is (4,6,4,1)",
      coarse_dims == (4, 6, 4, 1), str(coarse_dims))
check("barycentric tetrahedron f-vector is (15,50,60,24)",
      fine_dims == (15, 50, 60, 24), str(fine_dims))


def incidence(simplices):
    indices = [{cell: index for index, cell in enumerate(layer)}
               for layer in simplices]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(simplices[degree+1]),
                          len(simplices[degree]))
        for row, cell in enumerate(simplices[degree+1]):
            for omit in range(degree+2):
                face = cell[:omit] + cell[omit+1:]
                matrix[row, indices[degree][face]] = (-1)**omit
        differentials.append(matrix)
    return differentials


d_coarse = incidence(coarse)
d_fine = incidence(fine)
check("coarse and fine coboundaries square to zero exactly",
      all(d_coarse[k+1]*d_coarse[k]
          == sy.zeros(len(coarse[k+2]), len(coarse[k])) for k in range(2))
      and all(d_fine[k+1]*d_fine[k]
              == sy.zeros(len(fine[k+2]), len(fine[k])) for k in range(2)))


coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([[covector[index] for index in basis]
                          for covector in covectors]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    """Exact Whitney-k L2 Gram matrix on one affine tetrahedron."""
    affine = sy.Matrix.hstack(points[1]-points[0], points[2]-points[0],
                              points[3]-points[0])
    inverse = affine.inv()
    gradients = [-sum((sy.Matrix(inverse.row(row)).T for row in range(3)),
                      sy.zeros(3, 1))]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det())/6
    barycentric_second_moment = volume*(sy.ones(4, 4)+sy.eye(4))/20
    local_forms = list(combinations(range(4), degree+1))
    coefficient_matrices = []
    for form in local_forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree+1):
                covectors = [gradients[form[index]]
                             for index in range(degree+1)
                             if index != omitted]
                components = (factorial(degree)*(-1)**omitted
                              * wedge_components(covectors, degree))
                coefficients[:, form[omitted]] += components
        coefficient_matrices.append(coefficients)

    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficient_matrices):
        for col, right in enumerate(coefficient_matrices):
            mass[row, col] = sy.simplify(sum(
                (left[basis, :]*barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def assemble_mass(vertices, top_cells, simplices, degree):
    indices = {cell: index for index, cell in enumerate(simplices[degree])}
    result = sy.zeros(len(simplices[degree]), len(simplices[degree]))
    local_faces = list(combinations(range(4), degree+1))
    for top in top_cells:
        local = local_whitney_mass(tuple(vertices[vertex] for vertex in top),
                                   degree)
        for local_row, face_left in enumerate(local_faces):
            global_left = tuple(top[index] for index in face_left)
            row = indices[global_left]
            for local_col, face_right in enumerate(local_faces):
                global_right = tuple(top[index] for index in face_right)
                col = indices[global_right]
                result[row, col] += local[local_row, local_col]
    return sy.simplify(result)


mass_coarse = [assemble_mass(reference_vertices, coarse_top, coarse, degree)
               for degree in range(4)]
mass_fine = [assemble_mass(fine_vertices, fine_top, fine, degree)
             for degree in range(4)]
check("all eight Whitney mass matrices are exact symmetric Gram matrices",
      all(matrix == matrix.T for matrix in mass_coarse+mass_fine))

minimum_mass_eigenvalue = min(
    float(np.linalg.eigvalsh(np.asarray(matrix, dtype=float))[0])
    for matrix in mass_coarse+mass_fine
)
check("all Whitney mass matrices are positive definite [numerical audit]",
      minimum_mass_eigenvalue > 1e-10,
      f"minimum eigenvalue={minimum_mass_eigenvalue:.6e}")

# The coefficient of a coarse Whitney k-form on a fine k-simplex is its
# degree of freedom: the integral of that form over the fine simplex.  This
# integral is the determinant of the corresponding barycentric-coordinate
# minor.  It includes point evaluation at k=0.
inclusions = []
for degree in range(4):
    inclusion = sy.zeros(fine_dims[degree], coarse_dims[degree])
    for row, fine_simplex in enumerate(fine[degree]):
        barycentric_rows = [fine_barycentric[vertex]
                            for vertex in fine_simplex]
        for col, coarse_simplex in enumerate(coarse[degree]):
            minor = sy.Matrix([
                [barycentric[vertex] for vertex in coarse_simplex]
                for barycentric in barycentric_rows
            ])
            inclusion[row, col] = minor.det()
    inclusions.append(inclusion)

check("all four Whitney inclusions have full coarse rank",
      all(inclusions[degree].rank() == coarse_dims[degree]
          for degree in range(4)))

commuting_residuals = [
    d_fine[degree]*inclusions[degree]
    - inclusions[degree+1]*d_coarse[degree]
    for degree in range(3)
]
check("Whitney inclusion commutes with d in every degree exactly",
      all(residual == sy.zeros(*residual.shape)
          for residual in commuting_residuals))

isometry_residuals = [
    sy.simplify(inclusions[degree].T*mass_fine[degree]
                * inclusions[degree]-mass_coarse[degree])
    for degree in range(4)
]
check("Whitney inclusion is L2-isometric in every degree exactly",
      all(residual == sy.zeros(*residual.shape)
          for residual in isometry_residuals),
      "P_k^T M_f,k P_k=M_c,k for k=0,1,2,3")

hilbert_adjoints = [
    mass_coarse[degree].inv()*inclusions[degree].T*mass_fine[degree]
    for degree in range(4)
]
check("Hilbert adjoints are exact left inverses of all inclusions",
      all(hilbert_adjoints[degree]*inclusions[degree]
          == sy.eye(coarse_dims[degree]) for degree in range(4)))

# Coarse codifferentials use the exact Whitney metrics.  Fine inverses are not
# required to certify the compressed adjoint; the weak-form identity follows
# from dP=Pd and P^T M_f P=M_c and is checked directly.
delta_coarse = [
    mass_coarse[degree].inv()*d_coarse[degree].T*mass_coarse[degree+1]
    for degree in range(3)
]
compressed_d = [
    hilbert_adjoints[degree+1]*d_fine[degree]*inclusions[degree]
    for degree in range(3)
]
compressed_delta = [
    mass_coarse[degree].inv()*inclusions[degree].T
    * d_fine[degree].T*mass_fine[degree+1]*inclusions[degree+1]
    for degree in range(3)
]
check("d and its Whitney adjoint compress to their coarse operators exactly",
      all(compressed_d[degree] == d_coarse[degree]
          and compressed_delta[degree] == delta_coarse[degree]
          for degree in range(3)))

# Strong codifferential intertwining would make the inherited subspace
# reducing.  Multiplication by M_f avoids forming dense fine inverses:
# M_f(delta_f P-P delta_c) is the exact residual below.
weak_adjoint_leakage = [
    sy.simplify(d_fine[degree].T*mass_fine[degree+1]
                * inclusions[degree+1]
                - mass_fine[degree]*inclusions[degree]*delta_coarse[degree])
    for degree in range(3)
]
leakage_ranks = tuple(residual.rank() for residual in weak_adjoint_leakage)
check("strong codifferential leakage has exact ranks (3,3,0)",
      leakage_ranks == (3, 3, 0),
      f"exact weak-residual ranks={leakage_ranks}")
check("every adjoint leakage is invisible under coarse compression",
      all(inclusions[degree].T*weak_adjoint_leakage[degree]
          == sy.zeros(coarse_dims[degree], coarse_dims[degree+1])
          for degree in range(3)))

# Assemble the exact weak-form Dirac A=M D, the block metric, inclusion and
# form-parity grading.  These matrices avoid every inverse on the fine level.
def exact_weak_dirac(masses, differentials):
    offsets = np.cumsum((0,) + tuple(matrix.rows for matrix in masses))
    total = int(offsets[-1])
    metric = sy.diag(*masses)
    operator = sy.zeros(total, total)
    for degree, differential in enumerate(differentials):
        lo0, hi0 = offsets[degree], offsets[degree+1]
        lo1, hi1 = offsets[degree+1], offsets[degree+2]
        forward = masses[degree+1]*differential
        operator[lo1:hi1, lo0:hi0] = forward
        operator[lo0:hi0, lo1:hi1] = forward.T
    grading = sy.diag(*[
        sy.eye(mass.rows) if degree % 2 == 0 else -sy.eye(mass.rows)
        for degree, mass in enumerate(masses)
    ])
    return offsets, metric, operator, grading


offsets_c, exact_metric_c, exact_weak_c, grading_c = exact_weak_dirac(
    mass_coarse, d_coarse)
offsets_f, exact_metric_f, exact_weak_f, grading_f = exact_weak_dirac(
    mass_fine, d_fine)
inclusion_all = sy.diag(*inclusions)
check("weak Kahler--Dirac forms are exactly symmetric and form-odd",
      exact_weak_c == exact_weak_c.T
      and exact_weak_f == exact_weak_f.T
      and exact_weak_c*grading_c+grading_c*exact_weak_c
      == sy.zeros(sum(coarse_dims))
      and exact_weak_f*grading_f+grading_f*exact_weak_f
      == sy.zeros(sum(fine_dims)))
check("the all-degree metric and Dirac form compress exactly",
      inclusion_all.T*exact_metric_f*inclusion_all == exact_metric_c
      and inclusion_all.T*exact_weak_f*inclusion_all == exact_weak_c)
check("form parity itself is compatible with refinement",
      grading_f*inclusion_all == inclusion_all*grading_c
      and inclusion_all.T*exact_metric_f*grading_f*inclusion_all
      == exact_metric_c*grading_c)

# Locality is a support statement for the variational/weak generator.  Every
# nonzero kinetic matrix entry must connect simplices contained in a common
# top tetrahedron.  The strong coefficient operator M^-1 A need not retain
# this exact sparsity, a limitation recorded in the note.
fine_stars = []
for degree in range(4):
    fine_stars.append([
        {top_index for top_index, top in enumerate(fine_top)
         if set(simplex).issubset(top)}
        for simplex in fine[degree]
    ])

locality_violations = []
kinetic_nonzeros = 0
for degree in range(3):
    forward = mass_fine[degree+1]*d_fine[degree]
    for high in range(forward.rows):
        for low in range(forward.cols):
            if forward[high, low] != 0:
                kinetic_nonzeros += 1
                if not (fine_stars[degree+1][high]
                        & fine_stars[degree][low]):
                    locality_violations.append((degree, high, low))
check("the fine weak Dirac stencil is exactly simplex-star local",
      kinetic_nonzeros > 0 and not locality_violations,
      f"directed kinetic entries={kinetic_nonzeros}, violations=0")

# Assemble finite generalized Dirac eigenproblems A v=lambda M v, where
# A=M D is symmetric.  This is a numerical spectral audit of an otherwise
# exact rational construction.
def generalized_dirac(masses, differentials):
    offsets = np.cumsum((0,) + tuple(matrix.rows for matrix in masses))
    total = int(offsets[-1])
    metric = np.zeros((total, total))
    weak_dirac = np.zeros((total, total))
    for degree, mass in enumerate(masses):
        lo, hi = offsets[degree], offsets[degree+1]
        metric[lo:hi, lo:hi] = np.asarray(mass, dtype=float)
    for degree, differential in enumerate(differentials):
        lo0, hi0 = offsets[degree], offsets[degree+1]
        lo1, hi1 = offsets[degree+1], offsets[degree+2]
        forward = (np.asarray(masses[degree+1], dtype=float)
                   @ np.asarray(differential, dtype=float))
        weak_dirac[lo1:hi1, lo0:hi0] = forward
        weak_dirac[lo0:hi0, lo1:hi1] = forward.T
    return metric, weak_dirac


metric_c, weak_c = generalized_dirac(mass_coarse, d_coarse)
metric_f, weak_f = generalized_dirac(mass_fine, d_fine)
eigen_c = sla.eigh(weak_c, metric_c, eigvals_only=True)
eigen_f = sla.eigh(weak_f, metric_f, eigvals_only=True)
zero_c = int(np.count_nonzero(np.abs(eigen_c) < 1e-8))
zero_f = int(np.count_nonzero(np.abs(eigen_f) < 1e-8))
check("coarse and fine Kahler--Dirac kernels equal the tetrahedron Betti sum",
      zero_c == 1 and zero_f == 1,
      f"kernel dimensions coarse/fine={zero_c}/{zero_f}")
check("both finite Whitney--Dirac dynamics are spectrally nontrivial",
      np.max(np.abs(eigen_c)) > 1 and np.max(np.abs(eigen_f)) > 1,
      f"spectral radii={np.max(abs(eigen_c)):.6f}/{np.max(abs(eigen_f)):.6f}")
check("both spectra have the grading-forced plus/minus pairing",
      np.max(np.abs(eigen_c+eigen_c[::-1])) < 1e-9
      and np.max(np.abs(eigen_f+eigen_f[::-1])) < 1e-8)

# A single graded Hamiltonian combines propagation and a constant internal
# mass eigenvalue: H=cD+mu*gamma.  The exact anticommutation above proves
# H^2=c^2 D^2+mu^2 for arbitrary c,mu.  A fixed rational witness checks the
# resulting finite-time M-unitary evolution; it is not a physical fit.
speed_witness = 1.0
mass_witness = 2.0
grading_f_float = np.asarray(grading_f, dtype=float)
strong_d_f = sla.solve(metric_f, weak_f, assume_a="pos")
strong_h_f = speed_witness*strong_d_f + mass_witness*grading_f_float
mass_shell_residual = np.max(np.abs(
    strong_h_f@strong_h_f
    - speed_witness**2*(strong_d_f@strong_d_f)
    - mass_witness**2*np.eye(strong_h_f.shape[0])
))
check("graded Whitney Hamiltonian satisfies one mass shell numerically",
      mass_shell_residual < 1e-9,
      f"max residual={mass_shell_residual:.3e}")

time_witness = 0.137
evolution = sla.expm(-1j*time_witness*strong_h_f)
unitarity_residual = np.max(np.abs(
    evolution.conj().T@metric_f@evolution-metric_f
))
check("massive fine-level spectral evolution is metric-unitary and nontrivial",
      unitarity_residual < 1e-9
      and np.linalg.norm(evolution-np.eye(evolution.shape[0])) > 1,
      f"M-unitarity residual={unitarity_residual:.3e}")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: all-degree Whitney spaces nest isometrically and commute with d.")
print("DERIVED: the metric Kahler--Dirac compresses exactly at quadratic-form level.")
print("DERIVED: its weak stencil is simplex-star local and its evolution is nontrivial.")
print("DERIVED: top-degree codifferential intertwines exactly under refinement.")
print("DERIVED NEGATIVE: lower degrees have vertical adjoint leakage of ranks 3,3.")
print("DERIVED: both finite Galerkin Dirac dynamics are nontrivial with correct kernel.")
print("DERIVED CONDITIONAL: cD+gamma M obeys one mass shell for any supplied c,M.")
print("OPEN: certify controlled spectral convergence and causal continuum dynamics.")
raise SystemExit(0 if passed == tests else 1)
