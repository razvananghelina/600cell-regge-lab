#!/usr/bin/env python3
"""Exact audit of the Hopf projector fields as a round-S3 spin-two carrier.

Protocol commit 4004c25 froze the hypotheses, formulas, scope and falsifiers.
The expected TT/Casimir result was known before this verifier was written;
this is a structural derivation audit, not a blind target discovery.
"""

from itertools import product
import json
from pathlib import Path

import sympy as sy


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "hopf_spin2_tensor_carrier.json"
PROTOCOL_COMMIT = "4004c25"

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def exact_zero(value):
    return sy.simplify(sy.expand(value)) == 0


def exact_zero_matrix(matrix):
    return all(exact_zero(value) for value in matrix)


def frobenius(left, right):
    return sy.simplify(sy.trace(left.T * right))


def projector(vector):
    vector = sy.Matrix(vector)
    return sy.simplify(vector * vector.T / vector.dot(vector))


def matrix_key(matrix):
    return tuple(sy.radsimp(sy.simplify(value)) for value in matrix)


def vertex_index(vertices, candidate):
    for index, vertex in enumerate(vertices):
        if exact_zero_matrix(candidate - vertex):
            return index
    return None


def compose_permutations(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def generated_permutation_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = list(generators)
    while frontier:
        permutation = frontier.pop()
        if permutation in group:
            continue
        group.add(permutation)
        for generator in generators:
            frontier.append(compose_permutations(permutation, generator))
            frontier.append(compose_permutations(generator, permutation))
    return group


print("Hopf symmetric-tensor carrier on the unit round S3")

sqrt5 = sy.sqrt(5)
phi = (1 + sqrt5) / 2

# The fixed twelve icosahedral vertices give the six unoriented fivefold axes.
vertices = []
for first, second in product((1, -1), repeat=2):
    vertices.extend((
        sy.Matrix((0, first, second * phi)),
        sy.Matrix((first, second * phi, 0)),
        sy.Matrix((first * phi, 0, second)),
    ))

projector_by_key = {}
for vertex in vertices:
    candidate = projector(vertex)
    projector_by_key.setdefault(matrix_key(candidate), candidate)
projectors = tuple(projector_by_key.values())
centered = tuple(sy.simplify(item - sy.eye(3) / 3) for item in projectors)

check(
    "twelve vertices give exactly six rank-one unoriented projectors",
    len(vertices) == 12 and len(projectors) == 6
    and all(
        exact_zero_matrix(item * item - item)
        and exact_zero(sy.trace(item) - 1)
        for item in projectors
    ),
)

center_sum = sy.zeros(3)
for tensor in centered:
    center_sum += tensor
gram = sy.Matrix([
    [frobenius(left, right) for right in centered]
    for left in centered
])
expected_gram = sy.Matrix(
    6, 6,
    lambda row, column: (
        sy.Rational(2, 3) if row == column else -sy.Rational(2, 15)
    ),
)
flat_columns = [sy.Matrix(tensor).reshape(9, 1) for tensor in centered]
span_rank = sy.Matrix.hstack(*flat_columns).rank()
check(
    "the centered projectors form the exact rank-five regular simplex",
    center_sum == sy.zeros(3)
    and gram == expected_gram
    and span_rank == 5,
    "span=Sym^2_0(R^3), norm^2=2/3, cross=-2/15",
)

sqrt2 = sy.sqrt(2)
sqrt6 = sy.sqrt(6)
tensor_basis = (
    sy.diag(1, -1, 0) / sqrt2,
    sy.diag(1, 1, -2) / sqrt6,
    sy.Matrix(((0, 1, 0), (1, 0, 0), (0, 0, 0))) / sqrt2,
    sy.Matrix(((0, 0, 1), (0, 0, 0), (1, 0, 0))) / sqrt2,
    sy.Matrix(((0, 0, 0), (0, 0, 1), (0, 1, 0))) / sqrt2,
)
coordinates = tuple(sy.Matrix([
    frobenius(tensor, basis_tensor) for basis_tensor in tensor_basis
]) for tensor in centered)
frame = sy.zeros(5)
for coordinate in coordinates:
    frame += coordinate * coordinate.T

h_symbols = sy.symbols("h0:5", real=True)
generic_tensor = sy.Matrix((
    (h_symbols[0], h_symbols[2], h_symbols[3]),
    (h_symbols[2], h_symbols[1], h_symbols[4]),
    (h_symbols[3], h_symbols[4], -h_symbols[0] - h_symbols[1]),
))
reconstructed = sy.zeros(3)
for tensor in centered:
    reconstructed += sy.Rational(5, 4) * frobenius(
        generic_tensor, tensor
    ) * tensor
check(
    "the six tensors give the exact (4/5) tight frame and reconstruction",
    exact_zero_matrix(frame - sy.Rational(4, 5) * sy.eye(5))
    and exact_zero_matrix(reconstructed - generic_tensor),
    "H=(5/4) sum_i Tr(H T_i) T_i",
)

# Quaternionic left/right multiplication maps Im(H) orthogonally onto T_qS3.
w, xq, yq, zq = sy.symbols("w x y z", real=True)
q = sy.Matrix((w, xq, yq, zq))
radius_squared = sy.expand(q.dot(q))
left_frame = sy.Matrix((
    (-xq, -yq, -zq),
    (w, -zq, yq),
    (zq, w, -xq),
    (-yq, xq, w),
))
right_frame = sy.Matrix((
    (-xq, -yq, -zq),
    (w, zq, -yq),
    (-zq, w, xq),
    (yq, -xq, w),
))
check(
    "left and right quaternion frames are tangent orthogonal frames",
    exact_zero_matrix(left_frame.T * q)
    and exact_zero_matrix(right_frame.T * q)
    and exact_zero_matrix(
        left_frame.T * left_frame - radius_squared * sy.eye(3)
    )
    and exact_zero_matrix(
        right_frame.T * right_frame - radius_squared * sy.eye(3)
    ),
    "for |q|=1 both frames are orthonormal bases of T_qS3",
)

lift_controls = []
for frame_matrix in (left_frame, right_frame):
    for tensor in centered:
        lifted = sy.simplify(frame_matrix * tensor * frame_matrix.T)
        lift_controls.append(
            exact_zero_matrix(lifted - lifted.T)
            and exact_zero(sy.trace(lifted))
            and exact_zero_matrix(lifted * q)
        )
check(
    "all twelve handed lifts are symmetric, tracefree and tangent",
    all(lift_controls),
)

# Levi-Civita connection in a left-invariant unit frame on the unit S3.
connection = []
for derivative_index in range(3):
    matrix = sy.zeros(3)
    for input_index in range(3):
        for output_index in range(3):
            matrix[output_index, input_index] = sy.LeviCivita(
                derivative_index, input_index, output_index
            )
    connection.append(matrix)

connection_brackets_hold = True
for left_index in range(3):
    for right_index in range(3):
        expected = sy.zeros(3)
        for output_index in range(3):
            expected += sy.LeviCivita(
                left_index, right_index, output_index
            ) * connection[output_index]
        connection_brackets_hold &= exact_zero_matrix(
            connection[left_index] * connection[right_index]
            - connection[right_index] * connection[left_index]
            - expected
        )
check(
    "connection matrices realize the exact so(3) algebra",
    connection_brackets_hold,
)

covariant_derivatives = [
    sy.simplify(matrix * generic_tensor - generic_tensor * matrix)
    for matrix in connection
]
divergence = sy.Matrix([
    sy.simplify(sum(
        covariant_derivatives[derivative_index][derivative_index, column]
        for derivative_index in range(3)
    ))
    for column in range(3)
])
check(
    "every constant-frame symmetric tensor has zero divergence",
    exact_zero_matrix(divergence),
    "antisymmetric connection contracts to zero against H_ab=H_ba",
)

rough_laplacian = sy.zeros(3)
for derivative_index in range(3):
    rough_laplacian -= (
        connection[derivative_index] * covariant_derivatives[derivative_index]
        - covariant_derivatives[derivative_index] * connection[derivative_index]
    )
check(
    "the connection Laplacian is exactly 6 on Sym^2_0(R^3)",
    exact_zero_matrix(rough_laplacian - 6 * generic_tensor),
    "nabla*nabla H=6H for unit round S3",
)

actual_tensor_checks = []
for tensor in centered:
    derivatives = [
        matrix * tensor - tensor * matrix for matrix in connection
    ]
    tensor_divergence = sy.Matrix([
        sum(derivatives[index][index, column] for index in range(3))
        for column in range(3)
    ])
    tensor_laplacian = sy.zeros(3)
    for index in range(3):
        tensor_laplacian -= (
            connection[index] * derivatives[index]
            - derivatives[index] * connection[index]
        )
    actual_tensor_checks.append(
        exact_zero(sy.trace(tensor))
        and exact_zero_matrix(tensor_divergence)
        and exact_zero_matrix(tensor_laplacian - 6 * tensor)
    )
check(
    "the six Hopf tensors span five homogeneous TT eigenfields",
    all(actual_tensor_checks) and span_rank == 5,
)

# The sign reversal of the right-invariant connection leaves divergence and
# the double commutator unchanged.
right_derivatives = [
    -matrix * generic_tensor + generic_tensor * matrix
    for matrix in connection
]
right_divergence = sy.Matrix([
    sy.simplify(sum(
        right_derivatives[index][index, column] for index in range(3)
    ))
    for column in range(3)
])
right_laplacian = sy.zeros(3)
for index in range(3):
    right_laplacian -= (
        (-connection[index]) * right_derivatives[index]
        - right_derivatives[index] * (-connection[index])
    )
check(
    "the opposite-handed homogeneous space is also TT with eigenvalue 6",
    exact_zero_matrix(right_divergence)
    and exact_zero_matrix(right_laplacian - 6 * generic_tensor),
)

# Exact icosahedral rotation generators.  They permute the twelve vertices,
# generate all 60 rotations, and fix no tracefree symmetric tensor.
rotation_three = sy.Matrix(((0, 0, 1), (1, 0, 0), (0, 1, 0)))
rotation_five = sy.Matrix((
    ((sqrt5 - 1) / 4, -(sqrt5 + 1) / 4, sy.Rational(1, 2)),
    ((sqrt5 + 1) / 4, sy.Rational(1, 2), (sqrt5 - 1) / 4),
    (-sy.Rational(1, 2), (sqrt5 - 1) / 4, (sqrt5 + 1) / 4),
))
rotation_controls = (
    exact_zero_matrix(rotation_three.T * rotation_three - sy.eye(3))
    and exact_zero_matrix(rotation_five.T * rotation_five - sy.eye(3))
    and exact_zero(sy.det(rotation_three) - 1)
    and exact_zero(sy.det(rotation_five) - 1)
    and exact_zero_matrix(rotation_three**3 - sy.eye(3))
    and exact_zero_matrix(rotation_five**5 - sy.eye(3))
)

permutations = []
for rotation in (rotation_three, rotation_five):
    permutation = tuple(
        vertex_index(vertices, sy.simplify(rotation * vertex))
        for vertex in vertices
    )
    permutations.append(permutation)
icosahedral_group = generated_permutation_group(permutations)
check(
    "the exact order-three/order-five rotations generate the 60-element A5 action",
    rotation_controls
    and all(index is not None for permutation in permutations for index in permutation)
    and len(icosahedral_group) == 60,
    f"generated permutation-group order={len(icosahedral_group)}",
)

q_three = {
    w: sy.Rational(1, 2),
    xq: sy.Rational(1, 2),
    yq: sy.Rational(1, 2),
    zq: sy.Rational(1, 2),
}
q_five = {
    w: phi / 2,
    xq: 0,
    yq: -1 / (2 * phi),
    zq: -sy.Rational(1, 2),
}
relative_three = sy.simplify(
    (left_frame.T * right_frame).subs(q_three)
)
relative_five = sy.simplify(
    (left_frame.T * right_frame).subs(q_five)
)
check(
    "relative left/right frames realize the exact A5 generators",
    exact_zero_matrix(relative_three - rotation_three.T)
    and exact_zero_matrix(relative_five - rotation_five)
    and exact_zero(sum(value**2 for value in q_three.values()) - 1)
    and exact_zero(sum(value**2 for value in q_five.values()) - 1),
    "L_q^T R_q gives the adjoint icosahedral rotations",
)

invariance_equations = []
for rotation in (rotation_three, rotation_five):
    invariance_equations.extend(list(
        sy.simplify(rotation * generic_tensor * rotation.T - generic_tensor)
    ))
invariance_matrix, _ = sy.linear_eq_to_matrix(
    invariance_equations, h_symbols
)
check(
    "A5 fixes no nonzero tracefree symmetric tensor",
    invariance_matrix.rank() == 5,
    "the invariant subspace of Sym^2_0 is zero",
)
check(
    "left- and right-invariant tracefree homogeneous tensor spaces intersect trivially",
    invariance_matrix.rank() == 5,
    "a field constant in both frames would be A5-invariant",
)

# Intrinsic spin content.  x=cos(theta), and chi_l is the SO(3) spin-l
# character restricted to a rotation through theta.
character_variable = sy.symbols("X", real=True)


def spin_character(spin):
    return sy.expand(
        1 + 2 * sum(
            sy.chebyshevt(weight, character_variable)
            for weight in range(1, spin + 1)
        )
    )


characters = [spin_character(spin) for spin in range(4)]
exterior_character = sy.expand(2 * characters[0] + 2 * characters[1])
check(
    "the geometric exterior fibre is 2 V0 + 2 V1 and contains no intrinsic V2",
    exact_zero(exterior_character - (4 + 4 * character_variable))
    and not exact_zero(exterior_character - characters[2]),
    "Lambda*(R^3)=2 scalar + 2 vector; orbital l=2 does not change fibre spin",
)

double_angle_vector_character = sy.expand(
    1 + 2 * (2 * character_variable**2 - 1)
)
symmetric_square_character = sy.expand(
    (characters[1]**2 + double_angle_vector_character) / 2
)
check(
    "the symmetric square is V0 + V2 and its tracefree part is spin two",
    exact_zero(
        symmetric_square_character - characters[0] - characters[2]
    )
    and sy.Poly(characters[2], character_variable).degree() == 2,
)

# Scope/provenance guards from the already registered authoritative outputs.
regge_transfer = json.loads(
    (HERE / "whitney_regge_continuum_transfer.json").read_text()
)
gravity_hessian = json.loads(
    (HERE / "gravity_box4_full_hessian.json").read_text()
)
action_origin = json.loads(
    (HERE / "hopf_hessian_action_origin.json").read_text()
)
projector_certificate = json.loads(
    (HERE / "hopf_projector_cubic.json").read_text()
)
check(
    "the positive tensor result is not silently transferred to the fixed Regge metric",
    "exact flat metric is pushed forward"
    in regge_transfer["carrier"]["map_role"]
    and regge_transfer["carrier"]["map_role"]
    != "unit round metric",
)
check(
    "the existing edge-weight Hessian supplies no gauge zero mode",
    gravity_hessian["full_hessian"]["exact_inertia_from_psd_and_rank"]
    == {"positive": 720, "zero": 0, "negative": 0},
)
check(
    "the existing certified action still has no selected coupling to this carrier",
    "no label-Hessian coupling" in action_origin["verdict"]
    and "remain open" in projector_certificate["interpretation_boundary"],
)

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "provenance": (
        "post-recognition structural audit; expected TT/Casimir result was "
        "known before implementation"
    ),
    "phenomenological_target_used": False,
    "complete_hypotheses": {
        "metric": "unit round S3 with sectional curvature +1",
        "frames": ["L_q(v)=qv", "R_q(v)=vq"],
        "axes": "the six certified unoriented C10/fivefold axes",
        "connection": "[e_a,e_b]=2 epsilon_abc e_c; nabla_a e_b=epsilon_abc e_c",
        "not_assumed": [
            "Lorentzian time", "Einstein action", "metric fluctuation",
            "diffeomorphism quotient", "stress-energy source",
        ],
    },
    "tensor_carrier": {
        "local_bundle": "Sym^2_0(T*S3)",
        "real_fibre_dimension": 5,
        "six_frame_rank": span_rank,
        "frame_operator": "4/5 I_5",
        "reconstruction": "H=(5/4) sum_i Tr(H T_i) T_i",
        "tracefree": True,
        "divergence_free_for_constant_frame_coefficients": True,
        "connection_laplacian_eigenvalue_unit_round_S3": 6,
        "left_homogeneous_dimension": 5,
        "right_homogeneous_dimension": 5,
        "left_right_homogeneous_intersection_dimension": 0,
    },
    "representation_content": {
        "kahler_dirac_exterior_fibre": "2 V0 + 2 V1",
        "intrinsic_spin2_in_exterior_fibre": False,
        "symmetric_square_of_vector": "V0 + V2",
        "tracefree_symmetric_tensor": "V2",
        "orbital_warning": (
            "l=2 scalar/vector harmonics do not change intrinsic fibre spin"
        ),
    },
    "scope_boundary": {
        "round_result_transferred_to_fixed_regge": False,
        "existing_box4_edge_hessian_zero_modes": 0,
        "selected_action_or_source_coupling": False,
        "physical_graviton_derived": False,
    },
    "verdicts": [
        {
            "label": "DERIVED ROUND-S3 KINEMATICS",
            "claim": (
                "the six Hopf projector fields furnish the full symmetric-"
                "tracefree tangent fibre and five homogeneous TT eigenfields "
                "in each handed frame"
            ),
        },
        {
            "label": "DERIVED CORRECTION",
            "claim": (
                "the repository does contain a canonical spin-two tensor "
                "carrier; the old coexact one-form argument was the wrong carrier"
            ),
        },
        {
            "label": "OPEN PHYSICS",
            "claim": (
                "no selected gravitational action, gauge quotient, Lorentzian "
                "propagation, universal source coupling or Planck scale follows"
            ),
        },
    ],
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("machine-readable carrier certificate was written", OUTPUT.exists())

print(f"\nRESULT: {passed}/{tests} checks passed")
print("DERIVED: Hopf projectors give a genuine round-S3 TT tensor carrier.")
print("CORRECTION: this is the right spin-two type; coexact one-forms were not.")
print("OPEN: action, gauge symmetry, propagation, source coupling and scale.")
raise SystemExit(0 if passed == tests else 1)
