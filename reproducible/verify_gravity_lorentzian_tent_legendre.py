#!/usr/bin/env python3
"""Complex-action and pre/post Legendre audit for the Lorentzian tent.

Protocol commit: 9b10ed4.  The mixed-Hessian rank was not computed before
registration.  The calculation uses the plus complex-angle branch of
Borissova--Dittrich and covers the 12 old/12 new cone-edge star sector with
the 30 common link edges held fixed.
"""

import cmath
import itertools
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_lorentzian_tent_legendre.json"
REGULARITY_RESULT = HERE / "gravity_lorentzian_tent_regular_evolution.json"
PROTOCOL_COMMIT = "9b10ed4"
PROTOCOL_CORRECTION_COMMIT = "bd4eaa5"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


graph = nx.icosahedral_graph()
triangles = sorted(
    tuple(sorted(clique))
    for clique in nx.enumerate_all_cliques(graph)
    if len(clique) == 3
)
OLD = 12
NEW = 13
RHO_INDEX = 24
simplexes = [(OLD, NEW, *triangle) for triangle in triangles]
distances = nx.single_source_shortest_path_length(graph, 0)
shells = [
    sorted(node for node in graph if distances[node] == shell)
    for shell in range(4)
]


def projected_theta(rho, q1, q2, q3):
    """Previously certified real angle at internal hinge q1."""
    c1 = (1-rho-q1)/2
    c2 = (1-rho-q2)/2
    c3 = (1-rho-q3)/2
    denominator = 4*(c1*c1+rho)
    p11 = (4*c1*c1-4*c1*c2+4*c2*c2+3*rho)/denominator
    p22 = (4*c1*c1-4*c1*c3+4*c3*c3+3*rho)/denominator
    p12 = (
        2*c1*c1-2*c1*c2-2*c1*c3+4*c2*c3+rho
    )/denominator
    return math.acos(np.clip(p12/math.sqrt(p11*p22), -1.0, 1.0))


def internal_weight(rho, q):
    return (1+q+rho)/(4*math.sqrt(4*rho+(1-q-rho)**2))


def simple_pole_equation(rho, q_values):
    deficits = {node: 2*math.pi for node in graph}
    for triangle in triangles:
        for node in triangle:
            other = [vertex for vertex in triangle if vertex != node]
            deficits[node] -= projected_theta(
                rho, q_values[node], q_values[other[0]], q_values[other[1]]
            )
    return sum(
        deficits[node]*internal_weight(rho, q_values[node]) for node in graph
    )


def edge_square_and_variable(vertex_a, vertex_b, variables):
    """Return signed squared length, variable index and ds/d(variable)."""
    a, b = sorted((vertex_a, vertex_b))
    if (a, b) == (OLD, NEW):
        return -variables[RHO_INDEX], RHO_INDEX, -1.0
    if b == OLD and a < 12:
        return variables[a], a, 1.0
    if b == NEW and a < 12:
        return variables[12+a], 12+a, 1.0
    if a < 12 and b < 12 and graph.has_edge(a, b):
        return 1.0, None, 0.0
    raise ValueError(f"non-edge requested: {(vertex_a, vertex_b)}")


def simplex_distance_matrix(simplex, variables):
    matrix = np.zeros((5, 5), dtype=float)
    for i in range(5):
        for j in range(i+1, 5):
            value, _, _ = edge_square_and_variable(simplex[i], simplex[j], variables)
            matrix[i, j] = matrix[j, i] = value
    return matrix


def log_minus(value):
    """Principal logarithm with the negative-real cut approached from below."""
    scale = max(1.0, abs(value))
    if abs(value.imag) < 2e-13*scale:
        value = complex(value.real, 0.0)
        if value.real < 0:
            return complex(math.log(-value.real), -math.pi)
    return cmath.log(value)


def signed_volume_square(squared, local_vertices):
    """Signed n-volume squared of a simplex selected from a distance matrix."""
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return 1.0
    base = vertices[0]
    others = vertices[1:]
    gram = np.array([
        [
            (squared[base, i]+squared[base, j]-squared[i, j])/2
            for j in others
        ]
        for i in others
    ])
    return float(np.linalg.det(gram)/(math.factorial(dimension)**2))


def simplex_angles(simplex, variables):
    """All ten theta+ angles and diagnostic normal data."""
    squared = simplex_distance_matrix(simplex, variables)
    gram = np.array([
        [
            (squared[0, i]+squared[0, j]-squared[i, j])/2
            for j in range(1, 5)
        ]
        for i in range(1, 5)
    ])
    inverse = np.linalg.inv(gram)
    normal_gram = np.empty((5, 5), dtype=float)
    normal_gram[1:, 1:] = inverse
    normal_gram[0, 1:] = -inverse.sum(axis=0)
    normal_gram[1:, 0] = -inverse.sum(axis=1)
    normal_gram[0, 0] = inverse.sum()
    simplex_volume_square = signed_volume_square(squared, range(5))

    angles = {}
    arguments = {}
    for a, b in itertools.combinations(range(5), 2):
        hinge_vertices = [index for index in range(5) if index not in (a, b)]
        facet_a = [index for index in range(5) if index != a]
        facet_b = [index for index in range(5) if index != b]
        hinge_volume_square = signed_volume_square(squared, hinge_vertices)
        facet_a_volume_square = signed_volume_square(squared, facet_a)
        facet_b_volume_square = signed_volume_square(squared, facet_b)
        # A.38 fixes the complex cosine branch.  The familiar normalized
        # outward-normal formula differs by a sign when the two facet normals
        # have opposite causal type, so using it unqualified breaks
        # Schlaefli at spacelike boundary hinges.
        gram_derivative = np.zeros((4, 4))
        opposite_edge = {a, b}
        for i in range(1, 5):
            for j in range(1, 5):
                hit_0i = {0, i} == opposite_edge
                hit_0j = {0, j} == opposite_edge
                hit_ij = i != j and {i, j} == opposite_edge
                gram_derivative[i-1, j-1] = (
                    float(hit_0i)+float(hit_0j)-float(hit_ij)
                )/2
        volume_derivative = simplex_volume_square*np.trace(
            inverse@gram_derivative
        )
        cosine = 16*volume_derivative/(
            cmath.sqrt(complex(facet_a_volume_square))
            * cmath.sqrt(complex(facet_b_volume_square))
        )
        # A.39 of Borissova--Dittrich fixes the sign/branch of sine.  A bare
        # principal sqrt(1-cos^2) gets spacelike-hinge orientations wrong and
        # violates the complex Schlaefli identity.
        sine = -(4/3)*(
            cmath.sqrt(complex(hinge_volume_square))
            * cmath.sqrt(complex(simplex_volume_square))
        )/(
            cmath.sqrt(complex(facet_a_volume_square))
            * cmath.sqrt(complex(facet_b_volume_square))
        )
        if abs(cosine*cosine+sine*sine-1) > 2e-10:
            raise RuntimeError("complex sine/cosine identity failed")
        argument = cosine+1j*sine
        angle = -1j*log_minus(argument)
        hinge = tuple(sorted(
            simplex[index] for index in range(5) if index not in (a, b)
        ))
        angles[hinge] = angle
        arguments[hinge] = argument
    return angles, arguments, np.linalg.eigvalsh(gram)


def cosine_branch_audit(simplex, variables):
    """Compare A.38 with the tempting normalized-normal cosine."""
    squared = simplex_distance_matrix(simplex, variables)
    gram = np.array([
        [
            (squared[0, i]+squared[0, j]-squared[i, j])/2
            for j in range(1, 5)
        ]
        for i in range(1, 5)
    ])
    inverse = np.linalg.inv(gram)
    normal_gram = np.empty((5, 5), dtype=float)
    normal_gram[1:, 1:] = inverse
    normal_gram[0, 1:] = -inverse.sum(axis=0)
    normal_gram[1:, 0] = -inverse.sum(axis=1)
    normal_gram[0, 0] = inverse.sum()
    simplex_volume_square = signed_volume_square(squared, range(5))
    records = []
    for a, b in itertools.combinations(range(5), 2):
        facet_a = [index for index in range(5) if index != a]
        facet_b = [index for index in range(5) if index != b]
        va = signed_volume_square(squared, facet_a)
        vb = signed_volume_square(squared, facet_b)
        gram_derivative = np.zeros((4, 4))
        opposite_edge = {a, b}
        for i in range(1, 5):
            for j in range(1, 5):
                gram_derivative[i-1, j-1] = (
                    float({0, i} == opposite_edge)
                    + float({0, j} == opposite_edge)
                    - float(i != j and {i, j} == opposite_edge)
                )/2
        volume_derivative = simplex_volume_square*np.trace(
            inverse@gram_derivative
        )
        volume_cosine = 16*volume_derivative/(
            cmath.sqrt(complex(va))*cmath.sqrt(complex(vb))
        )
        normal_cosine = -normal_gram[a, b]/(
            cmath.sqrt(complex(normal_gram[a, a]))
            * cmath.sqrt(complex(normal_gram[b, b]))
        )
        same_causal_type = normal_gram[a, a]*normal_gram[b, b] > 0
        records.append((same_causal_type, volume_cosine, normal_cosine))
    return records


def hinge_type(hinge):
    vertices = set(hinge)
    if OLD in vertices and NEW in vertices:
        return "internal"
    if OLD in vertices:
        return "old_boundary"
    if NEW in vertices:
        return "new_boundary"
    return "link_corner"


def hinge_k(hinge, link_k=0):
    kind = hinge_type(hinge)
    if kind == "internal":
        return 2
    if kind == "link_corner":
        return link_k
    return 1


def triangle_area_square_and_gradient(hinge, variables):
    """Signed area square V_h and derivatives in all 25 variables."""
    edge_pairs = [(hinge[0], hinge[1]), (hinge[0], hinge[2]), (hinge[1], hinge[2])]
    edge_data = [edge_square_and_variable(a, b, variables) for a, b in edge_pairs]
    x, y, z = [entry[0] for entry in edge_data]
    area_square = (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16
    partials = [(y+z-x)/8, (x+z-y)/8, (x+y-z)/8]
    gradient = np.zeros(25, dtype=float)
    for partial, (_, index, signed_jacobian) in zip(partials, edge_data):
        if index is not None:
            gradient[index] += partial*signed_jacobian
    return area_square, gradient


def evaluate_action_and_gradient(variables, link_k=0, diagnostics=False):
    """Full complex action and Schlaefli-reduced gradient."""
    angle_incidence = {}
    minimum_argument = math.inf
    simplex_inertias = []
    angle_records = []
    for simplex in simplexes:
        angles, arguments, eigenvalues = simplex_angles(simplex, variables)
        simplex_inertias.append((eigenvalues < -1e-10).sum())
        for hinge, angle in angles.items():
            angle_incidence.setdefault(hinge, []).append(angle)
            minimum_argument = min(minimum_argument, abs(arguments[hinge]))
            if diagnostics:
                angle_records.append((simplex, hinge, angle))

    curvatures = {
        hinge: math.pi*hinge_k(hinge, link_k)+sum(angles)
        for hinge, angles in angle_incidence.items()
    }
    i_action = 0j
    gradient = np.zeros(25, dtype=complex)
    for hinge, curvature in curvatures.items():
        area_square, area_square_gradient = triangle_area_square_and_gradient(
            hinge, variables
        )
        root_area = cmath.sqrt(complex(area_square))
        i_action += root_area*curvature
        gradient += -1j*curvature*area_square_gradient/(2*root_area)
    action = -1j*i_action
    data = {
        "angle_incidence": angle_incidence,
        "curvatures": curvatures,
        "minimum_angle_argument_modulus": minimum_argument,
        "simplex_negative_eigenvalue_counts": simplex_inertias,
        "angle_records": angle_records,
    }
    return action, gradient, data


print("="*78)
print("LORENTZIAN TENT COMPLEX ACTION AND PRE/POST LEGENDRE MAP")
print("="*78)

# Re-isolate the previously certified shell root using the independent real
# internal-angle expression.
def q_map_for_x(x):
    shell_values = [x, 1.5, 0.8, 1.5]
    return {node: shell_values[distances[node]] for node in graph}


x_root = brentq(
    lambda x: simple_pole_equation(0.25, q_map_for_x(x)),
    11/25,
    9/20,
    xtol=5e-15,
    rtol=1e-14,
)
variables = np.empty(25, dtype=float)
variables[:12] = 1.0
variables[12:24] = [q_map_for_x(x_root)[node] for node in range(12)]
variables[RHO_INDEX] = 0.25

action, gradient, diagnostic = evaluate_action_and_gradient(
    variables, diagnostics=True
)
incidence_counts = [len(values) for values in diagnostic["angle_incidence"].values()]
type_counts = {
    kind: sum(
        1 for hinge in diagnostic["angle_incidence"] if hinge_type(hinge) == kind
    )
    for kind in ("internal", "old_boundary", "new_boundary", "link_corner")
}
check(
    "the complete 92-hinge tent census and Lorentzian simplex signatures are rebuilt",
    len(diagnostic["angle_incidence"]) == 92
    and sorted(type_counts.values()) == [12, 20, 30, 30]
    and incidence_counts.count(5) == 12
    and incidence_counts.count(2) == 60
    and incidence_counts.count(1) == 20
    and set(diagnostic["simplex_negative_eigenvalue_counts"]) == {1},
    f"types={type_counts}; incidences={{5:{incidence_counts.count(5)}, "
    f"2:{incidence_counts.count(2)}, 1:{incidence_counts.count(1)}}}",
)

# Compare every internal theta+ against the independently certified projected
# angle formula.  Identify which q belongs to the hinge and which two to the
# remaining link vertices of the simplex.
internal_angle_residuals = []
for simplex, hinge, angle in diagnostic["angle_records"]:
    if hinge_type(hinge) != "internal":
        continue
    node = next(vertex for vertex in hinge if vertex < 12)
    other = [vertex for vertex in simplex if vertex < 12 and vertex != node]
    expected = projected_theta(
        variables[RHO_INDEX],
        variables[12+node],
        variables[12+other[0]],
        variables[12+other[1]],
    )
    internal_angle_residuals.append(angle+expected)
check(
    "the complex plus branch reproduces every prior real internal angle",
    max(abs(value) for value in internal_angle_residuals) < 2e-14,
    f"max |theta+ + theta_real|={max(abs(value) for value in internal_angle_residuals):.3e}",
)

branch_records = [
    record for simplex in simplexes
    for record in cosine_branch_audit(simplex, variables)
]
same_type_residual = max(
    abs(volume-normal)
    for same, volume, normal in branch_records if same
)
opposite_type_residual = max(
    abs(volume+normal)
    for same, volume, normal in branch_records if not same
)
check(
    "A.38 agrees for same-type normals and flips the mixed-causal normal cosine",
    len(branch_records) == 200
    and same_type_residual < 3e-14
    and opposite_type_residual < 3e-14
    and any(not record[0] for record in branch_records),
    f"200 incidences; same residual={same_type_residual:.3e}; "
    f"opposite residual={opposite_type_residual:.3e}",
)

curvature_by_type = {
    kind: [
        value for hinge, value in diagnostic["curvatures"].items()
        if hinge_type(hinge) == kind
    ]
    for kind in type_counts
}
check(
    "all timelike curvatures are real and all spacelike boundary curvatures imaginary",
    max(abs(value.imag) for value in curvature_by_type["internal"]) < 3e-14
    and max(
        abs(value.real)
        for kind in ("old_boundary", "new_boundary", "link_corner")
        for value in curvature_by_type[kind]
    ) < 3e-14,
)
check(
    "the complete Lorentzian action and all Schlaefli derivatives are real",
    abs(action.imag) < 2e-10
    and np.max(np.abs(gradient.imag)) < 2e-10
    and diagnostic["minimum_angle_argument_modulus"] > 0.1,
    f"S={action}; max Im(grad)={np.max(np.abs(gradient.imag)):.3e}; "
    f"min |log argument|={diagnostic['minimum_angle_argument_modulus']:.6g}",
)

# Post-failure protocol correction: check Schlaefli simplex by simplex at the
# witness and two deterministic nearby causal configurations.  This directly
# differentiates every complex angle and does not use the reduced gradient.
perturbation_one = variables.copy()
perturbation_two = variables.copy()
perturbation_one[:12] += 4e-4*np.array([(-1)**index for index in range(12)])
perturbation_one[12:24] += 3e-4*np.array([(index % 3)-1 for index in range(12)])
perturbation_one[RHO_INDEX] += 5e-4
perturbation_two[:12] += 3e-4*np.array([(index % 4)-1.5 for index in range(12)])
perturbation_two[12:24] += 2e-4*np.array([1.5-(index % 4) for index in range(12)])
perturbation_two[RHO_INDEX] -= 4e-4
schlaefli_configurations = [variables, perturbation_one, perturbation_two]
schlaefli_residuals = []
schlaefli_step = 2e-6
perturbed_inertias = []
for configuration in schlaefli_configurations:
    for simplex in simplexes:
        _, _, eigenvalues = simplex_angles(simplex, configuration)
        perturbed_inertias.append(int((eigenvalues < -1e-10).sum()))
        relevant_indices = set()
        for local_a, local_b in itertools.combinations(range(5), 2):
            _, index, _ = edge_square_and_variable(
                simplex[local_a], simplex[local_b], configuration
            )
            if index is not None:
                relevant_indices.add(index)
        for index in relevant_indices:
            plus = configuration.copy()
            minus = configuration.copy()
            plus[index] += schlaefli_step
            minus[index] -= schlaefli_step
            plus_angles = simplex_angles(simplex, plus)[0]
            minus_angles = simplex_angles(simplex, minus)[0]
            residual = 0j
            for hinge in plus_angles:
                area_square = triangle_area_square_and_gradient(
                    hinge, configuration
                )[0]
                angle_derivative = (
                    plus_angles[hinge]-minus_angles[hinge]
                )/(2*schlaefli_step)
                residual += cmath.sqrt(complex(area_square))*angle_derivative
            schlaefli_residuals.append(residual)
check(
    "direct angle differences satisfy complex Schlaefli in every tested simplex",
    set(perturbed_inertias) == {1}
    and max(abs(value) for value in schlaefli_residuals) < 2e-7,
    f"{len(schlaefli_residuals)} directional identities at 3 configurations; "
    f"max residual={max(abs(value) for value in schlaefli_residuals):.3e}",
)

# The common-link corner k changes only a constant because all of its edge
# lengths are frozen.  This is checked directly rather than asserted.
action_link_one, gradient_link_one, _ = evaluate_action_and_gradient(
    variables, link_k=1
)
check(
    "the common-link corner k choice leaves every canonical derivative unchanged",
    np.max(np.abs(gradient_link_one-gradient)) < 2e-13,
    f"max gradient change={np.max(np.abs(gradient_link_one-gradient)):.3e}; "
    f"action shift={action_link_one-action}",
)

simple_E = simple_pole_equation(variables[RHO_INDEX], q_map_for_x(x_root))
check(
    "the full action pole derivative reproduces the prior internal equation",
    abs(gradient[RHO_INDEX].real-simple_E) < 2e-12
    and abs(simple_E) < 2e-12,
    f"S_rho={gradient[RHO_INDEX].real:.3e}; E={simple_E:.3e}",
)

# Centered differences of the action attack the Schlaefli-reduced gradient.
action_step = 1e-6
finite_action_gradient = np.zeros(25)
for index in range(25):
    plus = variables.copy()
    minus = variables.copy()
    plus[index] += action_step
    minus[index] -= action_step
    plus_action = evaluate_action_and_gradient(plus)[0]
    minus_action = evaluate_action_and_gradient(minus)[0]
    finite_action_gradient[index] = ((plus_action-minus_action)/(2*action_step)).real
gradient_real = gradient.real
gradient_error = np.max(
    np.abs(finite_action_gradient-gradient_real)
    / np.maximum(1.0, np.abs(gradient_real))
)
check(
    "centered differences of the full complex action reproduce all 25 derivatives",
    gradient_error < 2e-6,
    f"maximum declared relative error={gradient_error:.3e}",
)

pre_momenta = -gradient_real[:12]
post_momenta = gradient_real[12:24]


def shell_orbit_values(values):
    return [float(np.mean(values[shell])) for shell in shells]


def shell_spread(values):
    return max(float(np.ptp(values[shell])) for shell in shells)


check(
    "pre- and post-momenta are real and constant on all four stabilizer shells",
    np.all(np.isfinite(pre_momenta))
    and np.all(np.isfinite(post_momenta))
    and shell_spread(pre_momenta) < 2e-13
    and shell_spread(post_momenta) < 2e-13,
    f"P- shells={shell_orbit_values(pre_momenta)}; "
    f"P+ shells={shell_orbit_values(post_momenta)}",
)


def hessian_from_gradient(step):
    hessian = np.empty((25, 25), dtype=float)
    max_imaginary = 0.0
    for index in range(25):
        plus = variables.copy()
        minus = variables.copy()
        plus[index] += step
        minus[index] -= step
        plus_gradient = evaluate_action_and_gradient(plus)[1]
        minus_gradient = evaluate_action_and_gradient(minus)[1]
        column = (plus_gradient-minus_gradient)/(2*step)
        max_imaginary = max(max_imaginary, float(np.max(np.abs(column.imag))))
        hessian[:, index] = column.real
    return hessian, max_imaginary


hessian_steps = [2e-5, 1e-5, 5e-6]
hessians = []
hessian_imaginary = []
for step in hessian_steps:
    hessian, max_imaginary = hessian_from_gradient(step)
    hessians.append(hessian)
    hessian_imaginary.append(max_imaginary)

hessian_scale = max(1.0, float(np.max(np.abs(hessians[-1]))))
symmetry_residuals = [
    float(np.max(np.abs(hessian-hessian.T))/hessian_scale)
    for hessian in hessians
]
step_discrepancies = [
    float(np.max(np.abs(hessian-hessians[-1]))/hessian_scale)
    for hessian in hessians[:-1]
]
check(
    "the complete 25x25 Hessian is real, symmetric and step-stable",
    max(hessian_imaginary) < 2e-8
    and max(symmetry_residuals) < 2e-5
    and max(step_discrepancies) < 2e-4,
    f"symmetry={symmetry_residuals}; step discrepancies={step_discrepancies}; "
    f"max Im={max(hessian_imaginary):.3e}",
)

with REGULARITY_RESULT.open() as stream:
    regularity_certificate = json.load(stream)
certified_text = regularity_certificate["regularity"]["E_rho_at_root"]
certified_pole_hessian = float(certified_text.split("+/-")[0].lstrip("["))
pole_hessian_values = [hessian[RHO_INDEX, RHO_INDEX] for hessian in hessians]
pole_hessian_error = max(
    abs(value-certified_pole_hessian)/abs(certified_pole_hessian)
    for value in pole_hessian_values
)
check(
    "the full Hessian pole entry reproduces the independent Arb certificate",
    pole_hessian_error < 2e-6,
    f"H_rhorho={pole_hessian_values}; Arb midpoint={certified_pole_hessian:.15g}; "
    f"max rel={pole_hessian_error:.3e}",
)

# Eliminate the internal pole and form the on-shell old/new mixed block for
# every preregistered differentiation step.
mixed_blocks = []
singular_values = []
ranks = []
relative_singular_values = []
for hessian in hessians:
    symmetric = (hessian+hessian.T)/2
    boundary = symmetric[:24, :24]
    boundary_pole = symmetric[:24, RHO_INDEX]
    effective = boundary-np.outer(boundary_pole, boundary_pole)/symmetric[RHO_INDEX, RHO_INDEX]
    mixed = effective[:12, 12:24]
    values = np.linalg.svd(mixed, compute_uv=False)
    relative = values/values[0]
    mixed_blocks.append(mixed)
    singular_values.append(values)
    relative_singular_values.append(relative)
    ranks.append(int(np.count_nonzero(relative > 1e-7)))

reference_singular = singular_values[-1]
singular_stability = max(
    float(np.max(np.abs(values-reference_singular)/np.maximum(reference_singular, 1e-15)))
    for values in singular_values[:-1]
)
smallest_ratio = float(relative_singular_values[-1][-1])
near_threshold = 1e-8 <= smallest_ratio <= 1e-6
full_rank = ranks == [12, 12, 12] and not near_threshold
rank_deficient = len(set(ranks)) == 1 and ranks[0] < 12 and not near_threshold
check(
    "the on-shell mixed singular spectrum satisfies the preregistered stability rule",
    singular_stability < 5e-3
    and len(set(ranks)) == 1
    and not near_threshold,
    f"ranks={ranks}; s_min/s_max={smallest_ratio:.6e}; "
    f"max singular-value drift={singular_stability:.3e}",
)

# Independent post-correction confirmation: obtain only the raw blocks needed
# for W from second differences of the full action, never calling the
# Schlaefli gradient.  The correction protocol freezes h=2e-5.
direct_step = 2e-5


def action_real_at(configuration):
    return float(evaluate_action_and_gradient(configuration)[0].real)


def direct_mixed_second(index_a, index_b):
    plus_plus = variables.copy()
    plus_minus = variables.copy()
    minus_plus = variables.copy()
    minus_minus = variables.copy()
    plus_plus[index_a] += direct_step
    plus_plus[index_b] += direct_step
    plus_minus[index_a] += direct_step
    plus_minus[index_b] -= direct_step
    minus_plus[index_a] -= direct_step
    minus_plus[index_b] += direct_step
    minus_minus[index_a] -= direct_step
    minus_minus[index_b] -= direct_step
    return (
        action_real_at(plus_plus)-action_real_at(plus_minus)
        -action_real_at(minus_plus)+action_real_at(minus_minus)
    )/(4*direct_step*direct_step)


direct_pq = np.empty((12, 12))
for old_index in range(12):
    for new_node in range(12):
        direct_pq[old_index, new_node] = direct_mixed_second(
            old_index, 12+new_node
        )
direct_p_rho = np.array([
    direct_mixed_second(index, RHO_INDEX) for index in range(12)
])
direct_rho_q = np.array([
    direct_mixed_second(RHO_INDEX, 12+index) for index in range(12)
])
rho_plus = variables.copy()
rho_minus = variables.copy()
rho_plus[RHO_INDEX] += direct_step
rho_minus[RHO_INDEX] -= direct_step
direct_rho_rho = (
    action_real_at(rho_plus)-2*action.real+action_real_at(rho_minus)
)/(direct_step*direct_step)
direct_mixed = direct_pq-np.outer(
    direct_p_rho, direct_rho_q
)/direct_rho_rho
direct_block_error = float(
    np.linalg.norm(direct_mixed-mixed_blocks[0])
    /np.linalg.norm(mixed_blocks[0])
)
direct_singular = np.linalg.svd(direct_mixed, compute_uv=False)
direct_relative = direct_singular/direct_singular[0]
direct_rank = int(np.count_nonzero(direct_relative > 1e-7))
direct_singular_error = float(np.max(
    np.abs(direct_singular-singular_values[0])
    /np.maximum(singular_values[0], 1e-15)
))
independent_confirmation = (
    direct_block_error < 2e-4
    and direct_rank == 12
    and direct_singular_error < 5e-3
)
check(
    "direct full-action second differences reproduce the on-shell mixed block",
    direct_block_error < 2e-4,
    f"relative Frobenius error={direct_block_error:.3e}; "
    f"direct S_rhorho={direct_rho_rho:.12g}",
)
check(
    "the independent full-action singular spectrum confirms rank twelve",
    direct_rank == 12 and direct_singular_error < 5e-3,
    f"rank={direct_rank}; s_min/s_max={direct_relative[-1]:.6e}; "
    f"max singular error={direct_singular_error:.3e}",
)

check(
    "the star-sector on-shell Legendre-map rank has a decisive verdict",
    (full_rank and independent_confirmation) or rank_deficient,
    "full rank, independently confirmed" if full_rank and independent_confirmation else (
        f"stable rank deficiency {ranks[-1]}" if rank_deficient else "inconclusive"
    ),
)

if full_rank and independent_confirmation:
    verdict = "DERIVED REGULAR LOCAL LEGENDRE MAP"
elif rank_deficient:
    verdict = "DERIVED CONSTRAINED/DEGENERATE STAR MAP"
else:
    verdict = "INCONCLUSIVE"

result = {
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "scope": {
        "carrier": "one icosahedral 4D tent move",
        "variables": "12 old + 12 new squared cone lengths; fixed link; one eliminated pole",
        "complex_angle_branch": "Borissova-Dittrich plus branch",
        "volume_coefficient": 0,
    },
    "witness": {
        "x_root": x_root,
        "rho": 0.25,
        "action": action.real,
        "action_imaginary_residual": action.imag,
        "pole_gradient": gradient_real[RHO_INDEX],
    },
    "hinges": {
        "total": len(diagnostic["angle_incidence"]),
        "type_counts": type_counts,
        "incidence_counts": {
            "five": incidence_counts.count(5),
            "two": incidence_counts.count(2),
            "one": incidence_counts.count(1),
        },
        "minimum_log_argument_modulus": diagnostic["minimum_angle_argument_modulus"],
    },
    "momenta_squared_length_coordinates": {
        "pre_by_node": pre_momenta.tolist(),
        "post_by_node": post_momenta.tolist(),
        "pre_shell_orbits": shell_orbit_values(pre_momenta),
        "post_shell_orbits": shell_orbit_values(post_momenta),
    },
    "hessian": {
        "steps": hessian_steps,
        "symmetry_residuals": symmetry_residuals,
        "step_discrepancies": step_discrepancies,
        "pole_entries": pole_hessian_values,
        "pole_arb_midpoint": certified_pole_hessian,
    },
    "on_shell_mixed_block": {
        "ranks_at_relative_threshold_1e-7": ranks,
        "singular_values": [values.tolist() for values in singular_values],
        "relative_singular_values_finest": relative_singular_values[-1].tolist(),
        "smallest_over_largest": smallest_ratio,
        "maximum_relative_step_drift": singular_stability,
        "direct_action_rank": direct_rank,
        "direct_action_smallest_over_largest": float(direct_relative[-1]),
        "direct_action_block_relative_error": direct_block_error,
        "direct_action_singular_relative_error": direct_singular_error,
    },
    "verdict": {
        "star_sector_legendre_map": verdict,
        "full_global_evolution": "OPEN",
        "constraint_matching": "OPEN",
        "physical_clock": "OPEN",
    },
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(result, indent=2)+"\n")

print("-"*78)
print(f"Result: {passed}/{tests} checks passed")
print(f"Verdict: {verdict}")
print(f"Wrote {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
