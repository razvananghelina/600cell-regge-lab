#!/usr/bin/env python3
"""Complete dual-complex audit of the global 600-cell flex holonomy."""

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_global_flex_holonomy_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_global_flex_holonomy_adversarial_protocol.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_global_flex_holonomy_protocol.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_global_flex_holonomy.py"
PRIMARY_JSON = HERE / "gravity_600cell_global_flex_holonomy.json"
GLUING_RESULT = ROOT / "docs/gravity/gravity_600cell_two_frustum_face_gluing_result.md"
GLUING_ADV_JSON = HERE / "gravity_600cell_two_frustum_face_gluing_adversarial.json"

PROTOCOL_COMMIT = "5ae2648"
EXPECTED_HASHES = {
    "protocol": "9404da1f2a9ca5b1d7cf0038f81870aff8916a83f8014d7787ed14fa3915c325",
    "primary_protocol": "671cfcd02d902a8cc95969619c7ae9bdb3279efd4704ea210f00b0b337be66b1",
    "primary_source": "9e4c13cf944283fbe473c318853ac951701abe6ac7147c78f525a1de071d7120",
    "primary_json": "6852c4f0da3f747f178a697647bc0326a9668858ef414d0078668f2030875acf",
    "gluing_result": "b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3",
    "gluing_adversarial_json": "0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542",
}

ETA = sp.diag(1, 1, 1, -1)
NORMAL = sp.Matrix((0, 0, 0, 1))
CANONICAL = tuple(sp.Matrix(point) for point in (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
))
REPRESENTATIVES = ((1, 5), (2, 5), (3, 11))
TRACE_TARGET = sp.Rational(725, 243)

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


# Exact Q(phi), represented as a+b*phi with phi^2=phi+1.
def golden(a=0, b=0):
    return (Fraction(a), Fraction(b))


def gadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def gmul(left, right):
    a, b = left
    c, d = right
    return (a * c + b * d, a * d + b * c + b * d)


def gscale(value, scalar):
    scalar = Fraction(scalar)
    return (scalar * value[0], scalar * value[1])


def gdot(left, right):
    result = golden()
    for a, b in zip(left, right):
        result = gadd(result, gmul(a, b))
    return result


def permutation_parity(permutation):
    return sum(permutation[i] > permutation[j]
               for i in range(len(permutation))
               for j in range(i + 1, len(permutation))) % 2


def build_exact_complex():
    zero = golden()
    one = golden(1)
    half = golden(Fraction(1, 2))
    phi_half = golden(0, Fraction(1, 2))
    inverse_phi_half = golden(Fraction(-1, 2), Fraction(1, 2))
    vertices = set()
    for axis in range(4):
        for sign in (-1, 1):
            point = [zero] * 4
            point[axis] = gscale(one, sign)
            vertices.add(tuple(point))
    for signs in product((-1, 1), repeat=4):
        vertices.add(tuple(gscale(half, sign) for sign in signs))
    base = (phi_half, half, inverse_phi_half, zero)
    for permutation in permutations(range(4)):
        if permutation_parity(permutation):
            continue
        permuted = tuple(base[permutation[index]] for index in range(4))
        nonzero = tuple(index for index, value in enumerate(permuted)
                        if value != zero)
        for signs in product((-1, 1), repeat=3):
            point = list(permuted)
            for index, sign in zip(nonzero, signs):
                point[index] = gscale(point[index], sign)
            vertices.add(tuple(point))
    vertices = tuple(sorted(vertices))

    adjacency = [set() for _ in vertices]
    edges = []
    for left, right in combinations(range(len(vertices)), 2):
        if gdot(vertices[left], vertices[right]) == phi_half:
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
    tetrahedra = tuple(sorted(tetrahedra))

    face_to_tetrahedra = {}
    edge_to_tetrahedra = {}
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            face_to_tetrahedra.setdefault(tuple(sorted(face)), []).append(tetrahedron)
        for edge in combinations(tetrahedron, 2):
            edge_to_tetrahedra.setdefault(tuple(sorted(edge)), []).append(tetrahedron)
    dual = {tetrahedron: set() for tetrahedron in tetrahedra}
    for incident in face_to_tetrahedra.values():
        if len(incident) == 2:
            left, right = incident
            dual[left].add(right)
            dual[right].add(left)
    return {
        "vertices": vertices,
        "edges": tuple(edges),
        "triangles": tuple(triangles),
        "tetrahedra": tetrahedra,
        "face_to_tetrahedra": face_to_tetrahedra,
        "edge_to_tetrahedra": edge_to_tetrahedra,
        "dual": dual,
    }


def coordinates(tetrahedron):
    return {vertex: CANONICAL[index]
            for index, vertex in enumerate(tetrahedron)}


def reflect_across_face(point, face_points):
    first, second, third = face_points
    normal = (second - first).cross(third - first)
    return sp.simplify(
        point - 2 * normal * (normal.T * (point - first))[0]
        / (normal.T * normal)[0]
    )


def affine_map(domain, codomain):
    domain_h = sp.Matrix.hstack(*(point.col_join(sp.ones(1, 1))
                                  for point in domain))
    codomain_h = sp.Matrix.hstack(*(point.col_join(sp.ones(1, 1))
                                    for point in codomain))
    return sp.simplify(codomain_h * domain_h.inv())


def face_transition(source, target):
    source_coordinates = coordinates(source)
    target_coordinates = coordinates(target)
    shared = tuple(sorted(set(source) & set(target)))
    if len(shared) != 3:
        raise RuntimeError("face transition without a shared triangle")
    source_apex = next(vertex for vertex in source if vertex not in shared)
    target_apex = next(vertex for vertex in target if vertex not in shared)
    reflected = reflect_across_face(
        source_coordinates[source_apex],
        tuple(source_coordinates[vertex] for vertex in shared),
    )
    domain_vertices = tuple(target)
    domain = tuple(target_coordinates[vertex] for vertex in domain_vertices)
    codomain = tuple(
        source_coordinates[vertex] if vertex in shared else reflected
        for vertex in domain_vertices
    )
    transition = affine_map(domain, codomain)
    return transition


def edge_cycle(edge, incident, start, reverse=False):
    star = set(incident)
    neighbours = {tetrahedron: [] for tetrahedron in star}
    for left, right in combinations(star, 2):
        intersection = set(left) & set(right)
        if len(intersection) == 3 and set(edge).issubset(intersection):
            neighbours[left].append(right)
            neighbours[right].append(left)
    if any(len(values) != 2 for values in neighbours.values()):
        raise RuntimeError("edge star is not a cycle")
    choices = sorted(neighbours[start])
    current = choices[1 if reverse else 0]
    path = [start, current]
    previous = start
    while current != start:
        candidates = [value for value in neighbours[current] if value != previous]
        if len(candidates) != 1:
            raise RuntimeError("ambiguous edge cycle traversal")
        following = candidates[0]
        previous, current = current, following
        path.append(current)
        if len(path) > len(star) + 1:
            raise RuntimeError("edge cycle did not close")
    if len(path) != len(star) + 1:
        raise RuntimeError("incomplete edge cycle")
    return tuple(path)


def loop_product(path, transitions):
    result = sp.eye(4)
    for source, target in zip(path[:-1], path[1:]):
        result = sp.simplify(result * transitions[(source, target)])
    return result


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


def poincare_adjoint(spatial_affine):
    spatial_linear = spatial_affine[:3, :3]
    spatial_translation = spatial_affine[:3, 3]
    linear = sp.eye(4)
    linear[:3, :3] = spatial_linear
    translation = sp.Matrix((spatial_translation[0], spatial_translation[1],
                             spatial_translation[2], 0))
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


def fixed_dimension(adjoint_matrices, subspace):
    equations = sp.Matrix.vstack(*(
        (adjoint - sp.eye(10)) * subspace
        for adjoint in adjoint_matrices
    ))
    return kernel(equations).cols


paths = {
    "protocol": PROTOCOL,
    "primary_protocol": PRIMARY_PROTOCOL,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
    "gluing_result": GLUING_RESULT,
    "gluing_adversarial_json": GLUING_ADV_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all complete-dual-complex inputs have frozen provenance",
      provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
gluing = json.loads(GLUING_ADV_JSON.read_text())
upstream_ok = bool(
    primary["outcome"] == "GLOBAL_FLEX_SEED_KILLED_BY_HOLONOMY"
    and primary["passed"] == primary["tests"] == 12
    and all(record["two_edge_common_fixed_dimension"] == 0
            for record in primary["records"])
    and gluing["outcome"] == "ADVERSARIAL_TWO_FRUSTUM_DIAGONAL_ONLY"
    and gluing["passed"] == gluing["tests"] == 11
)
check("the primary holonomy and independent gluing results persist", upstream_ok)

complex_data = build_exact_complex()
tetrahedra = complex_data["tetrahedra"]
f_vector = (
    len(complex_data["vertices"]),
    len(complex_data["edges"]),
    len(complex_data["triangles"]),
    len(tetrahedra),
)
edge_incidence = Counter(map(len, complex_data["edge_to_tetrahedra"].values()))
face_incidence = Counter(map(len, complex_data["face_to_tetrahedra"].values()))
start = tetrahedra[0]
visited = {start}
queue = deque((start,))
while queue:
    current = queue.popleft()
    for neighbour in complex_data["dual"][current]:
        if neighbour not in visited:
            visited.add(neighbour)
            queue.append(neighbour)
incidence_ok = bool(
    f_vector == (120, 720, 1200, 600)
    and edge_incidence == Counter({5: 720})
    and face_incidence == Counter({2: 1200})
    and len(visited) == 600
)
check("the independent exact carrier and complete dual graph are correct",
      incidence_ok,
      f"f={f_vector}, edge incidence={dict(edge_incidence)}, "
      f"dual component={len(visited)}")

transitions = {}
transition_isometries = True
for face, incident in complex_data["face_to_tetrahedra"].items():
    left, right = incident
    for source, target in ((left, right), (right, left)):
        transition = face_transition(source, target)
        transitions[(source, target)] = transition
        linear = transition[:3, :3]
        transition_isometries &= bool(
            transition[3, :] == sp.Matrix([[0, 0, 0, 1]])
            and sp.simplify(linear.T * linear) == sp.eye(3)
            and abs(linear.det()) == 1
        )
    transition_isometries &= bool(
        sp.simplify(transitions[(left, right)]
                    * transitions[(right, left)]) == sp.eye(4)
    )
check("all 2400 directed face transitions are exact inverse isometries",
      transition_isometries and len(transitions) == 2400)

loop_records = []
loop_matrices = {}
all_loop_geometry = True
trace_counter = Counter()
for edge in complex_data["edges"]:
    incident = complex_data["edge_to_tetrahedra"][edge]
    base = min(incident)
    forward_path = edge_cycle(edge, incident, base, reverse=False)
    reverse_path = edge_cycle(edge, incident, base, reverse=True)
    forward = loop_product(forward_path, transitions)
    reverse = loop_product(reverse_path, transitions)
    loop_matrices[(edge, base, False)] = forward
    loop_matrices[(edge, base, True)] = reverse
    linear = forward[:3, :3]
    base_coordinates = coordinates(base)
    endpoint_control = all(
        sp.simplify(
            forward * base_coordinates[vertex].col_join(sp.ones(1, 1))
            - base_coordinates[vertex].col_join(sp.ones(1, 1))
        ) == sp.zeros(4, 1)
        for vertex in edge
    )
    local_ok = bool(
        len(forward_path) == len(reverse_path) == 6
        and sp.simplify(linear.T * linear) == sp.eye(3)
        and linear.det() == 1
        and sp.trace(linear) == TRACE_TARGET
        and endpoint_control
        and sp.simplify(forward * reverse) == sp.eye(4)
    )
    all_loop_geometry &= local_ok
    trace_counter[str(sp.trace(linear))] += 1
    loop_records.append({
        "edge": list(edge),
        "base_tetrahedron": list(base),
        "forward_path": [list(tetrahedron) for tetrahedron in forward_path],
        "trace": str(sp.trace(linear)),
        "determinant": str(linear.det()),
        "fixes_edge": endpoint_control,
        "reverse_is_inverse": sp.simplify(forward * reverse) == sp.eye(4),
    })

check("all 720 actual edge-star loops have the exact Regge holonomy",
      all_loop_geometry and len(loop_records) == 720,
      f"trace multiset={dict(trace_counter)}")

base_tetrahedron = tetrahedra[0]
base_edges = tuple(combinations(base_tetrahedron, 2))
base_forward = []
base_reverse = []
for edge in base_edges:
    incident = complex_data["edge_to_tetrahedra"][tuple(sorted(edge))]
    forward_path = edge_cycle(edge, incident, base_tetrahedron, reverse=False)
    reverse_path = edge_cycle(edge, incident, base_tetrahedron, reverse=True)
    base_forward.append(poincare_adjoint(loop_product(forward_path, transitions)))
    base_reverse.append(poincare_adjoint(loop_product(reverse_path, transitions)))

time_translation = sp.zeros(10, 1)
time_translation[9, 0] = 1
full_forward = kernel(sp.Matrix.vstack(*(
    adjoint - sp.eye(10) for adjoint in base_forward
)))
full_reverse = kernel(sp.Matrix.vstack(*(
    adjoint - sp.eye(10) for adjoint in base_reverse
)))
full_control = bool(
    full_forward.shape == full_reverse.shape == (10, 1)
    and same_space(full_forward, time_translation)
    and same_space(full_reverse, time_translation)
)
check("all six actual base-edge loops fix only time translation in Poincare",
      full_control)

records = []
single_controls = True
all_six_zero = True
reversal_controls = True
for scale, lapse in REPRESENTATIVES:
    local = local_kernel(scale, lapse)
    forward_single = [fixed_dimension((adjoint,), local)
                      for adjoint in base_forward]
    reverse_single = [fixed_dimension((adjoint,), local)
                      for adjoint in base_reverse]
    forward_common = fixed_dimension(tuple(base_forward), local)
    reverse_common = fixed_dimension(tuple(base_reverse), local)
    expected = 2 if scale == 1 else 1
    single_controls &= bool(
        forward_single == reverse_single == [expected] * 6
    )
    all_six_zero &= forward_common == 0
    reversal_controls &= reverse_common == forward_common
    records.append({
        "scale": scale,
        "lapse": lapse,
        "forward_single_edge_fixed_dimensions": forward_single,
        "reverse_single_edge_fixed_dimensions": reverse_single,
        "forward_six_edge_common_fixed_dimension": forward_common,
        "reverse_six_edge_common_fixed_dimension": reverse_common,
    })

check("each actual base-edge loop has the predicted one-edge fixed space",
      single_controls)
check("all six actual base-edge loops kill every local seed",
      all_six_zero)
check("reversing all six actual loops changes no local fixed dimension",
      reversal_controls)

swap = (1, 0, 2, 3)
relabel_map = affine_map(CANONICAL, tuple(CANONICAL[index] for index in swap))
relabel_adjoint = poincare_adjoint(relabel_map)
relabeled_forward = [
    sp.simplify(relabel_adjoint * adjoint * relabel_adjoint.inv())
    for adjoint in base_forward
]
relabeled_full = kernel(sp.Matrix.vstack(*(
    adjoint - sp.eye(10) for adjoint in relabeled_forward
)))
relabel_ok = bool(
    same_space(relabeled_full, relabel_adjoint * time_translation)
    and all(
        fixed_dimension(tuple(relabeled_forward),
                        relabel_adjoint * local_kernel(scale, lapse)) == 0
        for scale, lapse in REPRESENTATIVES
    )
)
check("an odd base-tetrahedron relabelling preserves the closure verdict",
      relabel_ok)

controls_ok = bool(
    provenance_ok and upstream_ok and incidence_ok
    and transition_isometries and all_loop_geometry
    and full_control and single_controls and reversal_controls and relabel_ok
)
killed = bool(controls_ok and all_six_zero)
survives = bool(
    controls_ok and any(record["forward_six_edge_common_fixed_dimension"] > 0
                        for record in records)
)

if not controls_ok:
    outcome = "ADVERSARIAL_GLOBAL_HOLONOMY_CONTROL_FAILED"
elif killed:
    outcome = "ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED"
elif survives:
    outcome = "ADVERSARIAL_GLOBAL_FLEX_SURVIVES"
else:
    outcome = "ADVERSARIAL_GLOBAL_HOLONOMY_OPEN"

allowed = {
    "ADVERSARIAL_GLOBAL_HOLONOMY_CONTROL_FAILED",
    "ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED",
    "ADVERSARIAL_GLOBAL_FLEX_SURVIVES",
    "ADVERSARIAL_GLOBAL_HOLONOMY_OPEN",
}
check("the adversarial hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "f_vector": list(f_vector),
    "directed_face_transitions": len(transitions),
    "dual_component_tetrahedra": len(visited),
    "edge_star_loop_census": len(loop_records),
    "loop_trace_multiset": dict(trace_counter),
    "full_poincare_six_edge_fixed_dimension": full_forward.cols,
    "base_tetrahedron": list(base_tetrahedron),
    "base_edges": [list(edge) for edge in base_edges],
    "records": records,
    "classification": {
        "complete_dual_complex_flex_closure": (
            "NO NONZERO LOCAL SEED" if killed else "OPEN"
        ),
        "global_infinitesimal_cellular_flex": (
            "REFUTED DERIVED EXACT" if killed else "OPEN"
        ),
        "global_infinitesimal_rigidity": (
            "ADVERSARIALLY CORROBORATED FOR THE COMPLETE LOCAL FLEX FAMILY"
            if killed else "OPEN"
        ),
        "finite_reconstruction_action_hessian_or_dynamics": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
print("all-edge trace census:", dict(trace_counter))
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"single={record['forward_single_edge_fixed_dimensions']}, "
        f"all-six={record['forward_six_edge_common_fixed_dimension']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)
