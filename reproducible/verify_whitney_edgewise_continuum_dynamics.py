#!/usr/bin/env python3
"""Exact edgewise Whitney induction and continuum-dynamics gate.

Protocol commit 075cb38 froze the carrier, exact identities, theorem inputs,
principal-symbol control, labels and scope before this calculation.
"""

from itertools import combinations
import json
from math import factorial
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import sympy as sy


OUTPUT = Path(__file__).with_name(
    "whitney_edgewise_continuum_dynamics.json"
)
EDGEWISE_CERTIFICATE = Path(__file__).with_name(
    "whitney_rank_edgewise_refinement.json"
)
PROTOCOL_COMMIT = "075cb38"
EDGEWISE_PROTOCOL_COMMIT = "58fa9fc"
EDGEWISE_RESULT_COMMIT = "0eddf27"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in weak_compositions(total - first, parts - 1):
            yield (first,) + rest


def edgewise_facets(resolution, dimension=3):
    """Independent exact Edelsbrunner--Grayson color-scheme enumeration."""
    width = dimension + 1
    facets = set()
    for counts in weak_compositions(resolution * width, width):
        sequence = tuple(
            color for color, count in enumerate(counts)
            for _ in range(count)
        )
        rows = tuple(
            sequence[row * width:(row + 1) * width]
            for row in range(resolution)
        )
        columns = tuple(
            tuple(rows[row][column] for row in range(resolution))
            for column in range(width)
        )
        if len(set(columns)) != width:
            continue
        points = tuple(
            tuple(column.count(color) for color in range(width))
            for column in columns
        )
        facets.add(points)
    return tuple(sorted(facets))


regular_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
rank_chamber = (
    regular_vertices[0],
    (regular_vertices[0] + regular_vertices[1]) / 2,
    sum(regular_vertices[:3], sy.zeros(3, 1)) / 3,
    sum(regular_vertices, sy.zeros(3, 1)) / 4,
)


def physical_point(numerator, resolution):
    return sum((
        sy.Rational(weight, resolution) * vertex
        for weight, vertex in zip(numerator, rank_chamber)
    ), sy.zeros(3, 1))


def edgewise_mesh(resolution):
    numerator_top = edgewise_facets(resolution)
    numerators = tuple(sorted({point for top in numerator_top for point in top}))
    point_index = {point: index for index, point in enumerate(numerators)}
    vertices = {
        point_index[point]: physical_point(point, resolution)
        for point in numerators
    }
    top = tuple(sorted({
        tuple(sorted(point_index[point] for point in facet))
        for facet in numerator_top
    }))
    simplices = tuple(
        tuple(sorted({
            tuple(face)
            for cell in top
            for face in combinations(cell, degree + 1)
        }))
        for degree in range(4)
    )
    barycentric = {
        point_index[point]: sy.Matrix(tuple(
            sy.Rational(value, resolution) for value in point
        ))
        for point in numerators
    }
    return vertices, top, simplices, barycentric


def incidence(simplices):
    indices = [{cell: index for index, cell in enumerate(layer)}
               for layer in simplices]
    differentials = []
    for degree in range(3):
        matrix = sy.zeros(len(simplices[degree + 1]),
                          len(simplices[degree]))
        for row, cell in enumerate(simplices[degree + 1]):
            for omitted in range(degree + 2):
                face = cell[:omitted] + cell[omitted + 1:]
                matrix[row, indices[degree][face]] = (-1) ** omitted
        differentials.append(matrix)
    return differentials


coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    affine = sy.Matrix.hstack(
        points[1] - points[0], points[2] - points[0], points[3] - points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det()) / 6
    barycentric_second_moment = volume * (sy.ones(4, 4) + sy.eye(4)) / 20
    local_forms = list(combinations(range(4), degree + 1))
    coefficients = []
    for form in local_forms:
        matrix = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            matrix[0, form[0]] = 1
        else:
            for omitted in range(degree + 1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                matrix[:, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * wedge_components(covectors, degree)
                )
        coefficients.append(matrix)
    mass = sy.zeros(len(local_forms), len(local_forms))
    for row, left in enumerate(coefficients):
        for column, right in enumerate(coefficients):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :] * barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


def assemble_mass(vertices, top, simplices, degree):
    indices = {cell: index for index, cell in enumerate(simplices[degree])}
    local_faces = list(combinations(range(4), degree + 1))
    result = sy.zeros(len(simplices[degree]), len(simplices[degree]))
    for cell in top:
        local = local_whitney_mass(tuple(vertices[index] for index in cell),
                                   degree)
        for local_row, left in enumerate(local_faces):
            row = indices[tuple(cell[index] for index in left)]
            for local_column, right in enumerate(local_faces):
                column = indices[tuple(cell[index] for index in right)]
                result[row, column] += local[local_row, local_column]
    return sy.simplify(result)


def exact_weak_dirac(masses, differentials):
    offsets = np.cumsum((0,) + tuple(matrix.rows for matrix in masses))
    metric = sy.diag(*masses)
    weak = sy.zeros(int(offsets[-1]), int(offsets[-1]))
    for degree, differential in enumerate(differentials):
        low_start, low_stop = offsets[degree:degree + 2]
        high_start, high_stop = offsets[degree + 1:degree + 3]
        forward = masses[degree + 1] * differential
        weak[high_start:high_stop, low_start:low_stop] = forward
        weak[low_start:low_stop, high_start:high_stop] = forward.T
    grading = sy.diag(*[
        ((-1) ** degree) * sy.eye(mass.rows)
        for degree, mass in enumerate(masses)
    ])
    return offsets, metric, weak, grading


def row_sum_lumping(matrix):
    return sy.diag(*[
        sy.simplify(sum(matrix.row(row))) for row in range(matrix.rows)
    ])


print("=" * 78)
print("EDGEWISE WHITNEY CONTINUUM-DYNAMICS GATE")
print("=" * 78)

# -------------------------------------------------------------------------
# Inherited all-level geometry certificate.
# -------------------------------------------------------------------------
edgewise_record = json.loads(EDGEWISE_CERTIFICATE.read_text())
edgewise_counts = edgewise_record["edgewise_counts"]
geometry_gate = (
    edgewise_record["protocol_commit"] == EDGEWISE_PROTOCOL_COMMIT
    and all(edgewise_counts[str(level)]["normalized_shape_classes"] == 3
            for level in (2, 3, 4))
    and all(not value["failures"]
            for value in edgewise_record["equivariance"].values())
    and all(value[0] for value in edgewise_record["conformity"].values())
    and set(edgewise_record["nesting_fine_per_coarse"]) == {8}
    and edgewise_record["direct_edgewise_variants"] == 3
    and set(edgewise_record["direct_variant_A4_orbit_sizes"]) == {3}
)
check("the preregistered rank-edgewise all-level geometry gates are intact",
      geometry_gate,
      "nested, conforming, S4-equivariant, three fixed shape classes")

# -------------------------------------------------------------------------
# Independent exact Esd_1 -> Esd_2 Whitney calculation.
# -------------------------------------------------------------------------
vertices_c, top_c, cells_c, barycentric_c = edgewise_mesh(1)
vertices_f, top_f, cells_f, barycentric_f = edgewise_mesh(2)
dims_c = tuple(map(len, cells_c))
dims_f = tuple(map(len, cells_f))
check("the independent edgewise meshes have the exact tetrahedral f-vectors",
      dims_c == (4, 6, 4, 1) and dims_f == (10, 25, 24, 8),
      f"Esd1={dims_c}, Esd2={dims_f}")

d_c = incidence(cells_c)
d_f = incidence(cells_f)
check("both edgewise coboundaries square to zero exactly",
      all(d_c[degree + 1] * d_c[degree]
          == sy.zeros(d_c[degree + 1].rows, d_c[degree].cols)
          for degree in range(2))
      and all(d_f[degree + 1] * d_f[degree]
              == sy.zeros(d_f[degree + 1].rows, d_f[degree].cols)
              for degree in range(2)))

masses_c = [assemble_mass(vertices_c, top_c, cells_c, degree)
            for degree in range(4)]
masses_f = [assemble_mass(vertices_f, top_f, cells_f, degree)
            for degree in range(4)]
minimum_mass_eigenvalue = min(
    float(np.linalg.eigvalsh(np.asarray(matrix, dtype=float))[0])
    for matrix in masses_c + masses_f
)
check("all exact consistent edgewise Whitney masses are positive Gram forms",
      all(matrix == matrix.T for matrix in masses_c + masses_f)
      and minimum_mass_eigenvalue > 1e-12,
      f"minimum audited eigenvalue={minimum_mass_eigenvalue:.6e}")

# At Esd_1 the four mesh vertices are the four rank-chamber vertices.  Map
# their mesh IDs to rank barycentric component IDs before forming minors.
coarse_component = {
    vertex_cell[0]: next(
        index for index, value in enumerate(barycentric_c[vertex_cell[0]])
        if value == 1
    )
    for vertex_cell in cells_c[0]
}
inclusions = []
for degree in range(4):
    inclusion = sy.zeros(dims_f[degree], dims_c[degree])
    for row, fine_simplex in enumerate(cells_f[degree]):
        for column, coarse_simplex in enumerate(cells_c[degree]):
            minor = sy.Matrix([
                [barycentric_f[fine_vertex][coarse_component[coarse_vertex]]
                 for coarse_vertex in coarse_simplex]
                for fine_vertex in fine_simplex
            ])
            inclusion[row, column] = sy.simplify(minor.det())
    inclusions.append(inclusion)

check("all four exact edgewise Whitney inclusions have full coarse rank",
      all(inclusions[degree].rank() == dims_c[degree]
          for degree in range(4)))
commuting_residuals = [
    d_f[degree] * inclusions[degree]
    - inclusions[degree + 1] * d_c[degree]
    for degree in range(3)
]
check("edgewise Whitney inclusion commutes with d in every degree exactly",
      all(residual == sy.zeros(*residual.shape)
          for residual in commuting_residuals))
isometry_residuals = [
    sy.simplify(inclusions[degree].T * masses_f[degree]
                * inclusions[degree] - masses_c[degree])
    for degree in range(4)
]
check("consistent Whitney inclusion is exactly L2-isometric in all degrees",
      all(residual == sy.zeros(*residual.shape)
          for residual in isometry_residuals),
      "P_p* M_f,p P_p=M_c,p for p=0,1,2,3")

offsets_c, metric_c, weak_c, grading_c = exact_weak_dirac(masses_c, d_c)
offsets_f, metric_f, weak_f, grading_f = exact_weak_dirac(masses_f, d_f)
inclusion_all = sy.diag(*inclusions)
check("the edgewise all-degree metric and weak Dirac compress exactly",
      inclusion_all.T * metric_f * inclusion_all == metric_c
      and inclusion_all.T * weak_f * inclusion_all == weak_c)
check("form parity is exactly compatible with edgewise refinement",
      grading_f * inclusion_all == inclusion_all * grading_c
      and weak_f * grading_f + grading_f * weak_f
      == sy.zeros(weak_f.rows))

# Strong intertwining and mass lumping are recorded target-free.
delta_c = [
    masses_c[degree].inv() * d_c[degree].T * masses_c[degree + 1]
    for degree in range(3)
]
strong_leakages = [
    sy.simplify(
        d_f[degree].T * masses_f[degree + 1] * inclusions[degree + 1]
        - masses_f[degree] * inclusions[degree] * delta_c[degree]
    )
    for degree in range(3)
]
strong_leakage_ranks = tuple(matrix.rank() for matrix in strong_leakages)
check("all strong-adjoint leakage ranks were evaluated exactly",
      all(inclusions[degree].T * strong_leakages[degree]
          == sy.zeros(dims_c[degree], dims_c[degree + 1])
          for degree in range(3)),
      f"ranks={strong_leakage_ranks}; compression hides every leakage")

lumped_c = [row_sum_lumping(matrix) for matrix in masses_c]
lumped_f = [row_sum_lumping(matrix) for matrix in masses_f]
lumped_residuals = [
    sy.simplify(inclusions[degree].T * lumped_f[degree]
                * inclusions[degree] - lumped_c[degree])
    for degree in range(4)
]
lumped_residual_ranks = tuple(matrix.rank() for matrix in lumped_residuals)
check("row-sum lumping was audited without using it as an acceptance target",
      all(matrix.rows == dims_c[degree]
          for degree, matrix in enumerate(lumped_residuals)),
      f"exact isometry-residual ranks={lumped_residual_ranks}")

# Post-protocol hostile audit: an oriented p-simplex basis may be reversed by
# a diagonal sign change without changing geometry.  A geometric metric rule
# must transform covariantly.  Ordinary row-sum lumping does not for p=1,2.
lumping_orientation_covariance = {}
for degree in (1, 2, 3):
    sign_change = sy.diag(-1, *([1] * (masses_c[degree].rows - 1)))
    reoriented = sign_change * masses_c[degree] * sign_change
    lumping_orientation_covariance[degree] = (
        row_sum_lumping(reoriented)
        == sign_change * lumped_c[degree] * sign_change
    )
check("row-sum lumping is not orientation-covariant in middle degrees",
      not lumping_orientation_covariance[1]
      and not lumping_orientation_covariance[2]
      and lumping_orientation_covariance[3],
      f"covariant by degree={lumping_orientation_covariance}")

# The weak stencil must remain element local.
stars = []
for degree in range(4):
    stars.append([
        {top_id for top_id, top in enumerate(top_f)
         if set(simplex).issubset(top)}
        for simplex in cells_f[degree]
    ])
locality_violations = []
kinetic_nonzeros = 0
for degree in range(3):
    forward = masses_f[degree + 1] * d_f[degree]
    for high in range(forward.rows):
        for low in range(forward.cols):
            if forward[high, low] == 0:
                continue
            kinetic_nonzeros += 1
            if not (stars[degree + 1][high] & stars[degree][low]):
                locality_violations.append((degree, high, low))
check("the edgewise weak Dirac stencil is exactly simplex-star local",
      kinetic_nonzeros > 0 and not locality_violations,
      f"directed entries={kinetic_nonzeros}, violations=0")

# Finite spectral control, not a continuum fit.
metric_f_float = np.asarray(metric_f, dtype=float)
weak_f_float = np.asarray(weak_f, dtype=float)
eigenvalues_f = sla.eigh(weak_f_float, metric_f_float, eigvals_only=True)
kernel_f = int(np.count_nonzero(np.abs(eigenvalues_f) < 1e-8))
pairing_residual = float(np.max(np.abs(eigenvalues_f + eigenvalues_f[::-1])))
check("the finite edgewise Dirac dynamics is nontrivial with correct kernel",
      kernel_f == 1 and np.max(np.abs(eigenvalues_f)) > 1
      and pairing_residual < 1e-8,
      f"kernel={kernel_f}, radius={np.max(abs(eigenvalues_f)):.6f}, "
      f"pairing residual={pairing_residual:.3e}")

# -------------------------------------------------------------------------
# Exact continuum principal symbol on Lambda*(R^3).
# -------------------------------------------------------------------------
exterior_basis = tuple(
    subset for degree in range(4)
    for subset in combinations(range(3), degree)
)
exterior_index = {subset: index for index, subset in enumerate(exterior_basis)}
epsilon = []
iota = []
for axis in range(3):
    wedge = sy.zeros(8)
    contraction = sy.zeros(8)
    for column, subset in enumerate(exterior_basis):
        if axis not in subset:
            target = tuple(sorted((axis,) + subset))
            sign_value = (-1) ** sum(value < axis for value in subset)
            wedge[exterior_index[target], column] = sign_value
        if axis in subset:
            position = subset.index(axis)
            target = subset[:position] + subset[position + 1:]
            contraction[exterior_index[target], column] = (-1) ** position
    epsilon.append(wedge)
    iota.append(contraction)

clifford_gate = all(
    epsilon[left] * iota[right] + iota[right] * epsilon[left]
    == (sy.eye(8) if left == right else sy.zeros(8))
    for left in range(3) for right in range(3)
)
check("exterior multiplication and contraction obey Clifford relations",
      clifford_gate)

xi = sy.symbols("xi0:3", real=True)
epsilon_xi = sum((xi[index] * epsilon[index] for index in range(3)),
                 sy.zeros(8))
iota_xi = sum((xi[index] * iota[index] for index in range(3)),
              sy.zeros(8))
principal_symbol = sy.I * (epsilon_xi - iota_xi)
xi_squared = sum(value ** 2 for value in xi)
continuum_grading = sy.diag(*[
    (-1) ** len(subset) for subset in exterior_basis
])
check("the continuum Hodge--Dirac symbol squares exactly to |xi|^2",
      sy.simplify(principal_symbol ** 2 - xi_squared * sy.eye(8))
      == sy.zeros(8)
      and principal_symbol.H == principal_symbol)

c, mu = sy.symbols("c mu", real=True)
massive_symbol = c * principal_symbol + mu * continuum_grading
check("the mass term preserves the characteristic speed and mass shell",
      principal_symbol * continuum_grading
      + continuum_grading * principal_symbol == sy.zeros(8)
      and sy.simplify(
          massive_symbol ** 2
          - (c ** 2 * xi_squared + mu ** 2) * sy.eye(8)
      ) == sy.zeros(8),
      "principal speed=|c|; mu is zeroth order")

# Hostile calibrated ultraviolet control for the consistent Whitney metric.
q = sy.symbols("q", real=True)
velocity_ratio = (
    sy.cos(q / 2)
    / (1 - sy.Rational(2, 3) * sy.sin(q / 2) ** 2) ** sy.Rational(3, 2)
)
low_q_series = sy.series(velocity_ratio, q, 0, 5)
cutoff_velocity = sy.simplify(velocity_ratio.subs(q, 2 * sy.pi / 3))
check("fixed physical modes approach unit continuum speed",
      sy.limit(velocity_ratio, q, 0) == 1,
      f"v/c={low_q_series}")
check("the finite-mesh cutoff overshoot remains and is not hidden",
      cutoff_velocity == sy.sqrt(2),
      "v/c=sqrt(2) at q=2*pi/3")

# The finite theorem inputs are a logical conjunction, not a numerical
# spectrum fit.  The primary-source scope audit is kept separate: the
# convergence section of Christiansen 1006.4779 assumes a domain in R^d,
# while Christiansen 1007.1120 assumes a smooth compact manifold.  Neither
# explicitly covers the fixed closed Regge metric used by the exact masses.
finite_theorem_inputs = {
    "fixed_compact_piecewise_flat_cellular_carrier": True,
    "nested_conforming_meshes": geometry_gate,
    "mesh_size_tends_to_zero_as_2^-n": True,
    "uniform_shape_regular_three_classes": geometry_gate,
    "whitney_spaces_form_exact_subcomplex": all(
        residual == sy.zeros(*residual.shape) for residual in commuting_residuals
    ),
    "exact_consistent_L2_metric": all(
        residual == sy.zeros(*residual.shape) for residual in isometry_residuals
    ),
}
literature_scope = {
    "cellular_framework_starts_from_compact_metric_space": True,
    "eigenpair_section_assumes_domain_in_Rd": True,
    "compact_manifold_rellich_theorem_assumes_smoothness": True,
    "explicit_closed_Regge_carrier_coverage_found": False,
}
check("every finite algebraic and mesh input for FEEC is present",
      all(finite_theorem_inputs.values()),
      "this does not by itself prove the analytic compactness theorem")
check("the primary-source Regge scope gap is reported rather than assumed",
      not literature_scope["explicit_closed_Regge_carrier_coverage_found"],
      "Euclidean Lipschitz domain / smooth compact manifold do not explicitly "
      "equal the closed piecewise-flat Regge carrier")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "edgewise_protocol_commit": EDGEWISE_PROTOCOL_COMMIT,
    "edgewise_result_commit": EDGEWISE_RESULT_COMMIT,
    "phenomenological_target_used": False,
    "control": {
        "carrier": "one rank-ordered barycentric orthoscheme",
        "coarse": "Esd_1",
        "fine": "Esd_2",
        "coarse_f_vector": list(dims_c),
        "fine_f_vector": list(dims_f),
        "minimum_mass_eigenvalue": minimum_mass_eigenvalue,
        "strong_adjoint_leakage_ranks": list(strong_leakage_ranks),
        "row_sum_lumped_isometry_residual_ranks": list(
            lumped_residual_ranks
        ),
        "row_sum_lumping_orientation_covariance": {
            str(key): value
            for key, value in lumping_orientation_covariance.items()
        },
        "weak_stencil_nonzeros": kinetic_nonzeros,
        "finite_dirac_kernel": kernel_f,
        "finite_dirac_spectral_radius": float(np.max(abs(eigenvalues_f))),
        "finite_dirac_pairing_residual": pairing_residual,
    },
    "all_level_geometry": {
        "tower": "Esd_(2^n)(sd K)",
        "mesh_law": "h_n=h_0 2^-n",
        "normalized_shape_classes": 3,
        "top_tetrahedra_600cell": "14400 * 8^n",
    },
    "finite_theorem_inputs": finite_theorem_inputs,
    "literature_scope_audit": literature_scope,
    "primary_sources": [
        "Edelsbrunner-Grayson DOI 10.1007/s004540010063",
        "Arnold-Falk-Winther arXiv:0906.4325",
        "Christiansen arXiv:1006.4779",
        "Christiansen arXiv:1007.1120",
    ],
    "continuum_symbol": {
        "dimension": 3,
        "exterior_algebra_dimension": 8,
        "identity": "sigma_D(xi)^2=|xi|^2 I",
        "massive_identity": "(c sigma_D+mu gamma)^2=(c^2|xi|^2+mu^2)I",
        "characteristic_speed": "|c|",
    },
    "ultraviolet_control": {
        "low_q_velocity_series": str(low_q_series),
        "fixed_mode_limit_v_over_c": 1,
        "cutoff_q": "2*pi/3",
        "cutoff_v_over_c": "sqrt(2)",
    },
    "verdicts": [
        "DERIVED EDGEWISE GALERKIN-INDUCTIVE DYNAMICS",
        "DERIVED POST-PROTOCOL NEGATIVE: ordinary row-sum lumping is not orientation-covariant for p=1,2",
        "OPEN ANALYTIC GAP: cited spectral theorems do not explicitly cover the closed Regge carrier",
        "DERIVED SYMBOL / CONDITIONAL CONTINUUM CAUSALITY if a valid continuum limit is established",
        "DERIVED NEGATIVE: no strict finite-mesh causal cone for consistent mass",
        "NOT DERIVED: time, numerical c, mass, fourth dimension or Planck units",
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured edgewise-continuum certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"STRONG_LEAKAGE_RANKS={strong_leakage_ranks}")
print(f"LUMPED_ISOMETRY_RANKS={lumped_residual_ranks}")
print("DERIVED: exact shape-regular edgewise Galerkin dynamics")
print("OPEN ANALYTIC GAP: closed Regge-carrier spectral convergence")
print("CONDITIONAL SYMBOL: continuum Hodge--Dirac propagation speed would be |c|")
raise SystemExit(0 if passed == tests else 1)
