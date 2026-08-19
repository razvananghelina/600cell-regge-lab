#!/usr/bin/env python3
"""Exact affine-holonomy obstruction to a global 600-cell frustum flex."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_global_flex_holonomy.json"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_global_flex_holonomy_prior_art.md"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_global_flex_holonomy_protocol.md"
GLUING_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"
GLUING_SOURCE = HERE / "verify_gravity_600cell_two_frustum_face_gluing.py"
GLUING_JSON = HERE / "gravity_600cell_two_frustum_face_gluing.json"
GLUING_ADV_SOURCE = HERE / "verify_gravity_600cell_two_frustum_face_gluing_adversarial.py"
GLUING_ADV_JSON = HERE / "gravity_600cell_two_frustum_face_gluing_adversarial.json"
FINITE_SOURCE = HERE / "verify_finite_regge_a2_hessian.py"
FINITE_JSON = HERE / "finite_regge_a2_hessian.json"

PROTOCOL_COMMIT = "19de8f0"
EXPECTED_HASHES = {
    "prior_art": "e5477823bc765d83cf812d393282ff8376c502d2967617226f42f1707474d056",
    "protocol": "671cfcd02d902a8cc95969619c7ae9bdb3279efd4704ea210f00b0b337be66b1",
    "gluing_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
    "gluing_source": "52636ae59bd4e4568df175e32b7c3aeae4fbfbc3d475d255131b6db671c41ae7",
    "gluing_json": "0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e",
    "gluing_adversarial_source": "b7a1f63e193aad50783929c8448ce99c18f1b50dc8e5ea27e3ed1102ec9dfa26",
    "gluing_adversarial_json": "0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542",
    "finite_source": "4419c409c66488a246fbfdd2ff8ba265cf835884635738cfcf1bef8eab9ae5b2",
    "finite_json": "9f78212270f2dd2b3f73a2dd914f497a0914fc1f114c0483b1ea65b12036a025",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
TETRA = tuple(sp.Matrix(point) for point in (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
))
EDGE_PAIRS = ((0, 1), (0, 2))
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
ORIGIN_SHIFT_3 = sp.Matrix((2, -1, 3))
ORIGIN_SHIFT_4 = sp.Matrix((2, -1, 3, 0))

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def kernel(matrix):
    vectors = matrix.nullspace()
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(matrix.cols, 0)


def same_space(left, right):
    return bool(
        left.rows == right.rows
        and left.rank() == right.rank()
        and left.row_join(right).rank() == left.rank()
    )


def parity(permutation):
    return sum(permutation[i] > permutation[j]
               for i in range(len(permutation))
               for j in range(i + 1, len(permutation))) % 2


def exact_600cell_incidence():
    root5 = sp.sqrt(5)
    phi = (1 + root5) / 2
    vertices = set()
    for axis in range(4):
        for sign in (-1, 1):
            point = [sp.Integer(0)] * 4
            point[axis] = sp.Integer(sign)
            vertices.add(tuple(point))
    for signs in product((-1, 1), repeat=4):
        vertices.add(tuple(sp.Rational(sign, 2) for sign in signs))
    base = (phi / 2, sp.Rational(1, 2), 1 / (2 * phi), sp.Integer(0))
    for permutation in permutations(range(4)):
        if parity(permutation):
            continue
        permuted = tuple(base[permutation[index]] for index in range(4))
        nonzero = tuple(index for index, value in enumerate(permuted)
                        if value != 0)
        for signs in product((-1, 1), repeat=3):
            point = list(permuted)
            for index, sign in zip(nonzero, signs):
                point[index] = sp.expand(sign * point[index])
            vertices.add(tuple(point))
    vertices = tuple(vertices)
    target = phi / 2
    adjacency = [set() for _ in vertices]
    edges = []
    for left, right in combinations(range(len(vertices)), 2):
        dot = sp.expand(sum(vertices[left][axis] * vertices[right][axis]
                            for axis in range(4)))
        if sp.simplify(dot - target) == 0:
            adjacency[left].add(right)
            adjacency[right].add(left)
            edges.append((left, right))
    triangles = []
    for left, right in edges:
        for third in adjacency[left] & adjacency[right]:
            if right < third:
                triangles.append((left, right, third))
    tetrahedra = []
    for a, b, c in triangles:
        for d in adjacency[a] & adjacency[b] & adjacency[c]:
            if c < d:
                tetrahedra.append((a, b, c, d))
    edge_incidence = Counter()
    for tetrahedron in tetrahedra:
        for edge in combinations(tetrahedron, 2):
            edge_incidence[tuple(sorted(edge))] += 1
    face_incidence = Counter()
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            face_incidence[tuple(sorted(face))] += 1
    return {
        "vertices": vertices,
        "edges": tuple(edges),
        "triangles": tuple(triangles),
        "tetrahedra": tuple(tetrahedra),
        "edge_incidence": edge_incidence,
        "face_incidence": face_incidence,
    }


def lorentz_basis():
    result = []
    for a, b in ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)):
        generator = sp.zeros(4)
        generator[a, b] = 1
        generator[b, a] = -ETA[a, a] / ETA[b, b]
        result.append(generator)
    return tuple(result)


LORENTZ = lorentz_basis()
LORENTZ_COORDINATES = sp.Matrix.hstack(
    *(generator.reshape(16, 1) for generator in LORENTZ)
)


def lorentz_coordinates(matrix):
    solution, free = LORENTZ_COORDINATES.gauss_jordan_solve(
        matrix.reshape(16, 1)
    )
    if free.rows:
        raise RuntimeError("ambiguous Lorentz coordinates")
    return solution


def cross_matrix(unit):
    x, y, z = unit
    return sp.Matrix(((0, -z, y), (z, 0, -x), (-y, x, 0)))


def axis_rotation(first, second, cosine, sine):
    direction = second - first
    unit = direction / sp.sqrt((direction.T * direction)[0])
    rotation = sp.simplify(
        cosine * sp.eye(3)
        + (1 - cosine) * (unit * unit.T)
        + sine * cross_matrix(unit)
    )
    translation = sp.simplify(first - rotation * first)
    return unit, rotation, translation


def poincare_adjoint(rotation3, translation3):
    linear = sp.eye(4)
    linear[:3, :3] = rotation3
    translation = sp.Matrix((translation3[0], translation3[1],
                             translation3[2], 0))
    result = sp.zeros(10)
    for column in range(10):
        if column < 6:
            A = LORENTZ[column]
            b = sp.zeros(4, 1)
        else:
            A = sp.zeros(4)
            b = sp.eye(4)[:, column - 6]
        transformed_A = sp.simplify(linear * A * linear.inv())
        transformed_b = sp.simplify(
            linear * b - transformed_A * translation
        )
        result[:6, column] = lorentz_coordinates(transformed_A)
        result[6:10, column] = transformed_b
    return result


def local_kernel(scale, lapse):
    if scale == 1:
        result = sp.zeros(10, 6)
        result[:3, :3] = sp.eye(3)
        result[6:9, 3:6] = sp.eye(3)
        return result
    result = sp.zeros(10, 6)
    result[:6, :6] = sp.eye(6)
    for column, generator in enumerate(LORENTZ):
        result[6:10, column] = (
            sp.Rational(lapse, scale - 1) * generator * NORMAL
        )
    return result


def origin_parameter_map(shift):
    result = sp.eye(10)
    for column, generator in enumerate(LORENTZ):
        result[6:10, column] = -generator * shift
    return result


def fixed_dimension(adjoint_matrices, subspace):
    equations = sp.Matrix.vstack(*(
        (adjoint - sp.eye(10)) * subspace
        for adjoint in adjoint_matrices
    ))
    return kernel(equations).cols, equations.rank()


paths = {
    "prior_art": PRIOR_ART,
    "protocol": PROTOCOL,
    "gluing_result": GLUING_RESULT,
    "gluing_source": GLUING_SOURCE,
    "gluing_json": GLUING_JSON,
    "gluing_adversarial_source": GLUING_ADV_SOURCE,
    "gluing_adversarial_json": GLUING_ADV_JSON,
    "finite_source": FINITE_SOURCE,
    "finite_json": FINITE_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all global-holonomy inputs have exact frozen provenance",
      provenance_ok, str(hashes))

gluing = json.loads(GLUING_JSON.read_text())
gluing_adversarial = json.loads(GLUING_ADV_JSON.read_text())
finite = json.loads(FINITE_JSON.read_text())
upstream_ok = bool(
    gluing["outcome"] == "TWO_FRUSTUM_DIAGONAL_ONLY"
    and gluing["passed"] == gluing["tests"] == 9
    and gluing_adversarial["outcome"]
    == "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY"
    and gluing_adversarial["passed"] == gluing_adversarial["tests"] == 11
    and finite["passed"] == finite["tests"] == 18
    and finite["f_vector"] == [120, 720, 1200, 600]
)
check("the accepted gluing and complete finite-Regge controls persist",
      upstream_ok)

complex_data = exact_600cell_incidence()
f_vector = tuple(len(complex_data[name]) for name in
                 ("vertices", "edges", "triangles", "tetrahedra"))
incidence_ok = bool(
    f_vector == (120, 720, 1200, 600)
    and complex_data["edge_incidence"] == Counter({edge: 5
                                                    for edge in complex_data["edges"]})
    and set(complex_data["face_incidence"].values()) == {2}
    and len(complex_data["face_incidence"]) == 1200
)
check("exact golden-field reconstruction gives five tetrahedra per edge",
      incidence_ok, f"f={f_vector}")

edge_vector = TETRA[1] - TETRA[0]
left = TETRA[2] - TETRA[0]
right = TETRA[3] - TETRA[0]
left_perp = sp.simplify(
    left - edge_vector * (left.T * edge_vector)[0]
    / (edge_vector.T * edge_vector)[0]
)
right_perp = sp.simplify(
    right - edge_vector * (right.T * edge_vector)[0]
    / (edge_vector.T * edge_vector)[0]
)
dihedral_cosine = sp.simplify(
    (left_perp.T * right_perp)[0]
    / sp.sqrt((left_perp.T * left_perp)[0]
              * (right_perp.T * right_perp)[0])
)
c = dihedral_cosine
chebyshev_five = sp.expand(16 * c**5 - 20 * c**3 + 5 * c)
u_four = sp.expand(16 * c**4 - 12 * c**2 + 1)
sin_theta = sp.sqrt(1 - c**2)
sin_five_theta = sp.simplify(sin_theta * u_four)
cosine_delta = sp.simplify(chebyshev_five)
sine_delta = sp.simplify(-sin_five_theta)
angle_ok = bool(
    dihedral_cosine == sp.Rational(1, 3)
    and cosine_delta == sp.Rational(241, 243)
    and sine_delta == 22 * sp.sqrt(2) / 243
    and sp.simplify(cosine_delta**2 + sine_delta**2) == 1
    and sine_delta > 0
)
check("the regular dihedral angle gives the exact fivefold deficit rotation",
      angle_ok,
      f"cos(theta)={dihedral_cosine}, cos(delta)={cosine_delta}, "
      f"sin(delta)={sine_delta}")

holonomies = []
shifted_holonomies = []
axis_controls = True
for left_index, right_index in EDGE_PAIRS:
    point = TETRA[left_index]
    endpoint = TETRA[right_index]
    unit, rotation, translation = axis_rotation(
        point, endpoint, cosine_delta, sine_delta
    )
    shifted_point = point + ORIGIN_SHIFT_3
    shifted_endpoint = endpoint + ORIGIN_SHIFT_3
    _, shifted_rotation, shifted_translation = axis_rotation(
        shifted_point, shifted_endpoint, cosine_delta, sine_delta
    )
    local_axis = bool(
        all(entry.is_rational for entry in rotation)
        and sp.simplify(rotation.T * rotation) == sp.eye(3)
        and rotation.det() == 1
        and sp.simplify(rotation * point + translation - point)
        == sp.zeros(3, 1)
        and sp.simplify(rotation * endpoint + translation - endpoint)
        == sp.zeros(3, 1)
        and shifted_rotation == rotation
        and sp.simplify(
            shifted_translation
            - (translation + ORIGIN_SHIFT_3 - rotation * ORIGIN_SHIFT_3)
        ) == sp.zeros(3, 1)
    )
    axis_controls &= local_axis
    holonomies.append(poincare_adjoint(rotation, translation))
    shifted_holonomies.append(
        poincare_adjoint(shifted_rotation, shifted_translation)
    )

nonparallel = bool(
    (TETRA[1] - TETRA[0]).cross(TETRA[2] - TETRA[0]) != sp.zeros(3, 1)
    and holonomies[0] * holonomies[1]
    != holonomies[1] * holonomies[0]
)
axis_controls &= nonparallel
check("two nonparallel affine edge rotations are exact and noncommuting",
      axis_controls)

full_common = kernel(sp.Matrix.vstack(
    holonomies[0] - sp.eye(10), holonomies[1] - sp.eye(10)
))
time_translation = sp.zeros(10, 1)
time_translation[9, 0] = 1
full_control_ok = bool(
    full_common.shape == (10, 1)
    and same_space(full_common, time_translation)
)
check("the two-edge full-Poincare fixed space is exactly time translation",
      full_control_ok)

records = []
zero_deficit_ok = True
single_edge_ok = True
two_edge_zero = True
orientation_ok = True
origin_ok = True

shift_parameter = origin_parameter_map(ORIGIN_SHIFT_4)
for scale, lapse in REPRESENTATIVES:
    local = local_kernel(scale, lapse)
    zero_dimension, _ = fixed_dimension((sp.eye(10),), local)
    one_dimensions = [fixed_dimension((adjoint,), local)[0]
                      for adjoint in holonomies]
    two_dimension, two_rank = fixed_dimension(tuple(holonomies), local)
    inverse_one = [fixed_dimension((adjoint.inv(),), local)[0]
                   for adjoint in holonomies]
    inverse_two, _ = fixed_dimension(
        tuple(adjoint.inv() for adjoint in holonomies), local
    )
    shifted_local = shift_parameter * local
    shifted_one = [fixed_dimension((adjoint,), shifted_local)[0]
                   for adjoint in shifted_holonomies]
    shifted_two, _ = fixed_dimension(tuple(shifted_holonomies), shifted_local)

    expected_one = 2 if scale == 1 else 1
    local_zero = zero_dimension == 6
    local_single = one_dimensions == [expected_one, expected_one]
    local_two = two_dimension == 0 and two_rank == 6
    local_orientation = bool(
        inverse_one == one_dimensions and inverse_two == two_dimension
    )
    local_origin = bool(
        shifted_one == one_dimensions and shifted_two == two_dimension
    )
    zero_deficit_ok &= local_zero
    single_edge_ok &= local_single
    two_edge_zero &= local_two
    orientation_ok &= local_orientation
    origin_ok &= local_origin

    records.append({
        "scale": scale,
        "lapse": lapse,
        "zero_deficit_fixed_dimension": zero_dimension,
        "single_edge_fixed_dimensions": one_dimensions,
        "two_edge_common_fixed_dimension": two_dimension,
        "two_edge_equation_rank_on_local_kernel": two_rank,
        "inverse_single_edge_fixed_dimensions": inverse_one,
        "inverse_two_edge_common_fixed_dimension": inverse_two,
        "shifted_single_edge_fixed_dimensions": shifted_one,
        "shifted_two_edge_common_fixed_dimension": shifted_two,
    })

check("zero deficit leaves all six local flex directions", zero_deficit_ok)
check("one hinge leaves the preregistered nonzero fixed dimensions",
      single_edge_ok)
check("two nonparallel hinges kill the local seed on every stratum",
      two_edge_zero)
check("reversing every dual loop leaves the fixed dimensions unchanged",
      orientation_ok)
check("shifting the complete development leaves all decisions unchanged",
      origin_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok and angle_ok
    and axis_controls and full_control_ok and zero_deficit_ok
    and single_edge_ok and orientation_ok and origin_ok
)
killed = bool(controls_ok and two_edge_zero)
survives = bool(
    controls_ok and any(record["two_edge_common_fixed_dimension"] > 0
                        for record in records)
)

if not controls_ok:
    outcome = "GLOBAL_FLEX_HOLONOMY_CONTROL_FAILED"
elif killed:
    outcome = "GLOBAL_FLEX_SEED_KILLED_BY_HOLONOMY"
elif survives:
    outcome = "GLOBAL_FLEX_SEED_SURVIVES_HOLONOMY"
else:
    outcome = "GLOBAL_FLEX_HOLONOMY_OPEN"

allowed = {
    "GLOBAL_FLEX_HOLONOMY_CONTROL_FAILED",
    "GLOBAL_FLEX_SEED_KILLED_BY_HOLONOMY",
    "GLOBAL_FLEX_SEED_SURVIVES_HOLONOMY",
    "GLOBAL_FLEX_HOLONOMY_OPEN",
}
check("the preregistered holonomy hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "independent_600cell_f_vector": list(f_vector),
    "tetrahedra_per_edge_multiset": {
        str(key): value for key, value in
        Counter(complex_data["edge_incidence"].values()).items()
    },
    "dihedral_cosine": str(dihedral_cosine),
    "deficit_cosine": str(cosine_delta),
    "deficit_sine": str(sine_delta),
    "full_poincare_two_edge_fixed_dimension": full_common.cols,
    "records": records,
    "classification": {
        "propagated_local_flex_seed": (
            "KILLED BY TWO NONPARALLEL REGGE HOLONOMIES"
            if killed else "OPEN"
        ),
        "global_infinite_small_cellular_flex": (
            "CANDIDATE REFUTATION" if killed else "OPEN"
        ),
        "global_infinitesimal_rigidity": (
            "CANDIDATE, REQUIRES COMPLETE-DUAL-COMPLEX AUDIT"
            if killed else "OPEN"
        ),
        "finite_uniqueness_action_hessian_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("deficit: cos=", cosine_delta, "sin=", sine_delta)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"one={record['single_edge_fixed_dimensions']}, "
        f"two={record['two_edge_common_fixed_dimension']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
