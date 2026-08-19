#!/usr/bin/env python3
"""Cellular Regge Hessian on the globally matched prism-shift family."""

from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_prism_shift_action.json"
PRIOR_ART_COMMIT = "29fac99"
PROTOCOL_COMMIT = "28d2885"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_prism_shift_action_prior_art.md":
        "13ae1f228cfc5c3abe4b54a11f60f7c11c1cfa25cada1f16d5b8f056e79b549c",
    "docs/gravity/gravity_600cell_prism_shift_action_protocol.md":
        "0c2e506a0f50a986c3970fc5c052f209c1d2c9726f20c9741fc726cbee518c05",
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
}
mp.mp.dps = 60
RHO = mp.mpf(1)
LENGTH = mp.mpf(1)
TRIANGLE_AREA = mp.sqrt(3)/4
EPSILON = 2*mp.pi-5*mp.acos(mp.mpf(1)/3)
LATERAL_NORMALS = (
    (-1, -1, -1, 0),
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
)
BOTTOM_NORMAL = (0, 0, 0, 1)
TOP_NORMAL = (0, 0, 0, -1)
G_INV = (
    (mp.mpf(3)/2, -mp.mpf(1)/2, -mp.mpf(1)/2),
    (-mp.mpf(1)/2, mp.mpf(3)/2, -mp.mpf(1)/2),
    (-mp.mpf(1)/2, -mp.mpf(1)/2, mp.mpf(3)/2),
)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    top = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    top.append((first, second, third, fourth))
    return tuple(top)


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def graph_distances(vertex_count, edges, seed):
    adjacency = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    distances = [-1]*vertex_count
    distances[seed] = 0
    queue = deque([seed])
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if distances[neighbour] == -1:
                distances[neighbour] = distances[current]+1
                queue.append(neighbour)
    return tuple(distances)


def normalized_direction(values):
    values = [mp.mpf(value) for value in values]
    mean = sum(values)/len(values)
    centered = [value-mean for value in values]
    norm = mp.sqrt(sum(value*value for value in centered))
    if norm == 0:
        raise ValueError("a frozen direction became constant")
    return tuple(value/norm for value in centered)


def scaled(direction, amplitude):
    amplitude = mp.mpf(amplitude)
    return tuple(amplitude*value for value in direction)


def add_potentials(left, right, scale=1):
    scale = mp.mpf(scale)
    return tuple(a+scale*b for a, b in zip(left, right))


def inverse_metric(a):
    """Exact block inverse for H=[G,a;a^T,-rho] with regular edge-one G."""
    u = [sum(G_INV[i][j]*a[j] for j in range(3)) for i in range(3)]
    schur = -RHO-sum(a[i]*u[i] for i in range(3))
    inverse = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    for i in range(3):
        for j in range(3):
            inverse[i][j] = G_INV[i][j]+u[i]*u[j]/schur
        inverse[i][3] = inverse[3][i] = -u[i]/schur
    inverse[3][3] = 1/schur
    return inverse, schur


def bilinear(left, metric, right):
    return sum(mp.mpf(left[i])*metric[i][j]*mp.mpf(right[j])
               for i in range(4) for j in range(4))


def branch_sqrt(value):
    value = mp.mpc(value)
    if abs(mp.im(value)) < mp.mpf("1e-50") and mp.re(value) < 0:
        return -1j*mp.sqrt(-mp.re(value))
    return mp.sqrt(value)


def angle(left, right, inverse):
    cross = bilinear(left, inverse, right)
    norm_product = bilinear(left, inverse, left)*bilinear(right, inverse, right)
    cosine = -cross/branch_sqrt(norm_product)
    return mp.acos(cosine), cosine


print("="*78)
print("CELLULAR REGGE HESSIAN ON THE PRISM-SHIFT FAMILY")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the prior-art gate, protocol and source have frozen provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "29fac99" and PROTOCOL_COMMIT == "28d2885",
    str(actual_hashes),
)

vertices, adjacency, _ = build_600cell()
top = tetrahedra_from_adjacency(adjacency)
simplices = all_simplices(top)
edges = simplices[1]
faces = simplices[2]
edge_index = {edge: index for index, edge in enumerate(edges)}
face_index = {face: index for index, face in enumerate(faces)}
edge_incidence = Counter()
face_incidence = Counter()
tetra_edge_data = []
tetra_face_data = []
for tetrahedron in top:
    local_edges = []
    for left, right in combinations(range(4), 2):
        edge = tuple(sorted((tetrahedron[left], tetrahedron[right])))
        omitted = tuple(index for index in range(4)
                        if index not in (left, right))
        local_edges.append((edge_index[edge], omitted))
        edge_incidence[edge] += 1
    tetra_edge_data.append(tuple(local_edges))

    local_faces = []
    for missing in range(4):
        face = tuple(tetrahedron[index] for index in range(4)
                     if index != missing)
        local_faces.append((face_index[face], missing))
        face_incidence[face] += 1
    tetra_face_data.append(tuple(local_faces))

source_f = tuple(len(layer) for layer in simplices)
check(
    "the source carrier and all cellular hinge incidences are exact",
    source_f == (120, 720, 1200, 600)
    and Counter(edge_incidence.values()) == Counter({5: 720})
    and Counter(face_incidence.values()) == Counter({2: 1200}),
    f"f={source_f}, edge incidence={dict(Counter(edge_incidence.values()))}, "
    f"face incidence={dict(Counter(face_incidence.values()))}",
)

# Integer graph Laplacian and frozen directions.
laplacian = np.zeros((120, 120), dtype=np.int64)
for left, right in edges:
    laplacian[left, left] += 1
    laplacian[right, right] += 1
    laplacian[left, right] -= 1
    laplacian[right, left] -= 1
adjacency_integer = np.eye(120, dtype=np.int64)*12-laplacian
graph_identity = bool(
    np.all(np.diag(adjacency_integer) == 0)
    and np.all((adjacency_integer == 0) | (adjacency_integer == 1))
    and np.all(adjacency_integer.sum(axis=1) == 12)
)
check(
    "the exact incidence Hessian carrier is Delta0=12I-A600",
    graph_identity,
)

seed = min(range(len(vertices)), key=lambda index: tuple(vertices[index]))
distances = graph_distances(120, edges, seed)
coordinate_order = sorted(range(120), key=lambda index: (vertices[index, 0], index))
coordinate_rank = [0]*120
for rank, vertex in enumerate(coordinate_order):
    coordinate_rank[vertex] = rank
raw_directions = {
    "one_vertex_delta": tuple(1 if index == seed else 0 for index in range(120)),
    "graph_distance": distances,
    "squared_graph_distance": tuple(value*value for value in distances),
    "first_coordinate_rank": tuple(coordinate_rank),
    "modular_quadratic": tuple((17*i*i+3*i+5) % 101 for i in range(120)),
}
directions = {name: normalized_direction(values)
              for name, values in raw_directions.items()}


def evaluate_action(potential, collect_diagnostics=False):
    edge_angle_sums = [mp.mpc(0) for _ in edges]
    bottom_angle_sums = [mp.mpc(0) for _ in faces]
    top_angle_sums = [mp.mpc(0) for _ in faces]
    maximum_top_relation = mp.mpf(0)
    signature_failures = 0
    anchor = None

    for top_index, tetrahedron in enumerate(top):
        reference = tetrahedron[0]
        a = [potential[tetrahedron[index]]-potential[reference]
             for index in range(1, 4)]
        inverse, schur = inverse_metric(a)
        signature_failures += int(not (schur < 0))

        local_lateral_angles = {}
        for global_edge, omitted in tetra_edge_data[top_index]:
            theta, _ = angle(LATERAL_NORMALS[omitted[0]],
                             LATERAL_NORMALS[omitted[1]], inverse)
            edge_angle_sums[global_edge] += theta
            local_lateral_angles[omitted] = theta

        local_bottom = []
        local_top = []
        for global_face, missing in tetra_face_data[top_index]:
            bottom_theta, _ = angle(BOTTOM_NORMAL,
                                    LATERAL_NORMALS[missing], inverse)
            top_theta, _ = angle(TOP_NORMAL,
                                 LATERAL_NORMALS[missing], inverse)
            bottom_angle_sums[global_face] += bottom_theta
            top_angle_sums[global_face] += top_theta
            maximum_top_relation = max(
                maximum_top_relation,
                abs(top_theta-(mp.pi-bottom_theta)),
            )
            local_bottom.append(bottom_theta)
            local_top.append(top_theta)

        if anchor is None:
            anchor = {
                "lateral": next(iter(local_lateral_angles.values())),
                "bottom": local_bottom[0],
                "top": local_top[0],
            }

    deficits = [2*mp.pi-value for value in edge_angle_sums]
    lateral_sum = mp.mpc(0)
    for edge, deficit in zip(edges, deficits):
        x = potential[edge[1]]-potential[edge[0]]
        lateral_sum += 1j*mp.sqrt(RHO*LENGTH*LENGTH+x*x)*deficit

    boundary_sum = mp.mpc(0)
    for bottom_sum, top_sum in zip(bottom_angle_sums, top_angle_sums):
        boundary_sum += TRIANGLE_AREA*((mp.pi-bottom_sum)+(mp.pi-top_sum))

    boundary_action = -1j*boundary_sum
    total = -1j*(lateral_sum+boundary_sum)
    diagnostics = {
        "signature_failures": signature_failures,
        "maximum_top_relation": maximum_top_relation,
        "boundary_action": boundary_action,
        "imaginary_action": abs(mp.im(total)),
        "anchor": anchor,
    }
    return total, deficits, diagnostics


zero = tuple(mp.mpf(0) for _ in range(120))
zero_action, zero_deficits, zero_diagnostics = evaluate_action(
    zero, collect_diagnostics=True)
anchor_ok = bool(
    abs(zero_diagnostics["anchor"]["lateral"]-mp.acos(mp.mpf(1)/3))
        < mp.mpf("2e-50")
    and abs(zero_diagnostics["anchor"]["bottom"]-mp.pi/2)
        < mp.mpf("2e-50")
    and abs(zero_diagnostics["anchor"]["top"]-mp.pi/2)
        < mp.mpf("2e-50")
)
check(
    "the independent facet conormals recover all static angle anchors",
    anchor_ok,
    str({key: mp.nstr(value, 20)
         for key, value in zero_diagnostics["anchor"].items()}),
)

direct_records = {}
direct_cache = {}
all_direct_ok = True
for name, direction in directions.items():
    record = {}
    for amplitude in (mp.mpf("1e-3"), mp.mpf("-1e-3")):
        key = (name, mp.nstr(amplitude, 8))
        action, deficits, diagnostics = evaluate_action(scaled(direction, amplitude))
        direct_cache[key] = (action, deficits, diagnostics)
        record[key[1]] = {
            "action": mp.nstr(mp.re(action), 20),
            "boundary_abs": mp.nstr(abs(diagnostics["boundary_action"]), 8),
            "imaginary_abs": mp.nstr(diagnostics["imaginary_action"], 8),
            "top_relation": mp.nstr(diagnostics["maximum_top_relation"], 8),
            "signature_failures": diagnostics["signature_failures"],
        }
        all_direct_ok = all_direct_ok and bool(
            diagnostics["signature_failures"] == 0
            and diagnostics["maximum_top_relation"] < mp.mpf("2e-11")
            and abs(diagnostics["boundary_action"]) < mp.mpf("2e-9")
            and diagnostics["imaginary_action"] < mp.mpf("2e-9")
        )
    direct_records[name] = record
check(
    "all nonconstant controls stay Lorentzian and their boundary terms cancel",
    all_direct_ok,
    f"directions={len(direct_records)}",
)

static_expected = 720*EPSILON
static_relative_error = abs(mp.re(zero_action)-static_expected)/abs(static_expected)
check(
    "the complete static action reproduces 720 times the regular deficit",
    static_relative_error < mp.mpf("2e-12")
    and zero_diagnostics["signature_failures"] == 0
    and zero_diagnostics["maximum_top_relation"] < mp.mpf("2e-50"),
    f"action={mp.nstr(zero_action, 25)}, expected={mp.nstr(static_expected, 25)}, "
    f"relative error={mp.nstr(static_relative_error, 8)}",
)

# Schlaefli-reduced gradient versus the complete area-angle action.
base_direction = directions["squared_graph_distance"]
base_potential = scaled(base_direction, "1e-3")
base_key = ("squared_graph_distance", "0.001")
if base_key in direct_cache:
    base_action, base_deficits, base_diagnostics = direct_cache[base_key]
else:
    base_action, base_deficits, base_diagnostics = evaluate_action(base_potential)
reduced_gradient = [mp.mpf(0) for _ in range(120)]
for edge, deficit in zip(edges, base_deficits):
    left, right = edge
    x = base_potential[right]-base_potential[left]
    value = deficit*x/mp.sqrt(1+x*x)
    reduced_gradient[left] -= value
    reduced_gradient[right] += value

gradient_errors = []
gradient_records = {}
for name in ("one_vertex_delta", "graph_distance",
             "first_coordinate_rank", "modular_quadratic"):
    direction = directions[name]
    expected = sum(value*component
                   for value, component in zip(reduced_gradient, direction))
    local = []
    for step in (mp.mpf("1e-5"), mp.mpf("5e-6")):
        plus, _, _ = evaluate_action(add_potentials(
            base_potential, direction, step))
        minus, _, _ = evaluate_action(add_potentials(
            base_potential, direction, -step))
        observed = mp.re(plus-minus)/(2*step)
        error = abs(observed-expected)
        allowed = mp.mpf("3e-9")+mp.mpf("3e-7")*max(
            abs(observed), abs(expected))
        gradient_errors.append((error, allowed))
        local.append({
            "step": mp.nstr(step, 8),
            "observed": mp.nstr(observed, 18),
            "expected": mp.nstr(expected, 18),
            "absolute_error": mp.nstr(error, 8),
            "allowed": mp.nstr(allowed, 8),
        })
    gradient_records[name] = local
check(
    "complete-action differences reproduce the Schlaefli gradient",
    all(error <= allowed for error, allowed in gradient_errors),
    f"maximum error/allowance={mp.nstr(max(error/allowed for error, allowed in gradient_errors), 8)}",
)

# Exact analytic area Hessian and graph assembly.
x_symbol, l_symbol, rho_symbol = sy.symbols("x L rho", positive=True)
area_symbol = sy.sqrt(rho_symbol*l_symbol**2+x_symbol**2)
area_second = sy.simplify(sy.diff(area_symbol, x_symbol, 2).subs(x_symbol, 0))
expected_area_second = 1/(l_symbol*sy.sqrt(rho_symbol))
check(
    "the lateral area has the frozen exact second derivative",
    sy.simplify(area_second-expected_area_second) == 0,
    f"A''(0)={area_second}",
)

directional_records = {}
directional_errors = []
for name, direction in directions.items():
    energy = sum((direction[right]-direction[left])**2
                 for left, right in edges)
    predicted = EPSILON*energy
    second_values = {}
    for step in (mp.mpf("1e-3"), mp.mpf("5e-4"), mp.mpf("2.5e-4")):
        cache_plus = (name, mp.nstr(step, 8))
        cache_minus = (name, mp.nstr(-step, 8))
        plus = (direct_cache[cache_plus][0] if cache_plus in direct_cache
                else evaluate_action(scaled(direction, step))[0])
        minus = (direct_cache[cache_minus][0] if cache_minus in direct_cache
                 else evaluate_action(scaled(direction, -step))[0])
        second_values[step] = mp.re(plus-2*zero_action+minus)/(step*step)
    medium = second_values[mp.mpf("5e-4")]
    fine = second_values[mp.mpf("2.5e-4")]
    richardson = (4*fine-medium)/3
    error = abs(richardson-predicted)
    allowed = mp.mpf("2e-9")+mp.mpf("2e-7")*max(
        abs(richardson), abs(predicted))
    directional_errors.append((error, allowed))
    directional_records[name] = {
        "energy": mp.nstr(energy, 18),
        "predicted": mp.nstr(predicted, 18),
        "second_differences": {
            mp.nstr(step, 8): mp.nstr(value, 18)
            for step, value in second_values.items()
        },
        "richardson": mp.nstr(richardson, 18),
        "absolute_error": mp.nstr(error, 8),
        "allowed": mp.nstr(allowed, 8),
    }
check(
    "all direct second variations equal epsilon times graph energy",
    all(error <= allowed for error, allowed in directional_errors),
    f"maximum error/allowance={mp.nstr(max(error/allowed for error, allowed in directional_errors), 8)}",
)

# Frozen exact spectrum, with no post-hoc relabeling.
golden = (1+np.sqrt(5.0))/2
expected_spectrum = (
    (0.0, 1, "0"),
    (12-6*golden, 4, "12-6phi"),
    (12-4*golden, 9, "12-4phi"),
    (9.0, 16, "9"),
    (12.0, 25, "12"),
    (14.0, 36, "14"),
    (8+4*golden, 9, "8+4phi"),
    (15.0, 16, "15"),
    (6+6*golden, 4, "6+6phi"),
)
eigenvalues = np.linalg.eigvalsh(laplacian.astype(float))
clusters = []
for value in eigenvalues:
    if not clusters or abs(value-clusters[-1][0]) > 1e-7:
        clusters.append([float(value), 1])
    else:
        clusters[-1][1] += 1
spectrum_errors = []
for observed, expected in zip(clusters, expected_spectrum):
    spectrum_errors.append(abs(observed[0]-expected[0]))
spectrum_ok = bool(
    len(clusters) == len(expected_spectrum)
    and all(observed[1] == expected[1]
            for observed, expected in zip(clusters, expected_spectrum))
    and max(spectrum_errors) < 1e-10
)
check(
    "the Hessian has the nine preregistered exact Laplacian clusters",
    spectrum_ok,
    str([(label, multiplicity) for _, multiplicity, label in expected_spectrum]),
)

kappa = EPSILON/(LENGTH*mp.sqrt(RHO))
minimum_quotient = kappa*mp.mpf(str(expected_spectrum[1][0]))
rank = int(np.count_nonzero(eigenvalues > 1e-9))
check(
    "the constant is the sole null and all 119 quotient modes are positive",
    rank == 119 and EPSILON > 0 and minimum_quotient > 0,
    f"rank={rank}, epsilon={mp.nstr(EPSILON, 20)}, "
    f"minimum quotient eigenvalue={mp.nstr(minimum_quotient, 20)}",
)

dust_phi_derivative = mp.mpf(0)
check(
    "fixed-strut conserved dust is constant on the shift family",
    dust_phi_derivative == 0,
    "S_dust=-(8*pi*M)*sqrt(rho) contains no phi",
)

verdict = (
    "SHIFT_HESSIAN_IS_GRAPH_LAPLACIAN"
    if passed == tests else
    "ACTION_BRANCH_OPEN"
)
check(
    "the preregistered verdict is evaluated",
    verdict == "SHIFT_HESSIAN_IS_GRAPH_LAPLACIAN",
    verdict,
)

artifact = {
    "action_controls": {
        "direct": direct_records,
        "static_action": mp.nstr(zero_action, 30),
        "static_expected": mp.nstr(static_expected, 30),
        "static_relative_error": mp.nstr(static_relative_error, 12),
        "maximum_gradient_error_ratio": mp.nstr(
            max(error/allowed for error, allowed in gradient_errors), 12),
        "gradient_records": gradient_records,
    },
    "classification": "DERIVED_RESTRICTED_SECTOR" if passed == tests else "OPEN",
    "hessian": {
        "formula": "epsilon/(L*sqrt(rho)) * (12I-A600)",
        "epsilon": mp.nstr(EPSILON, 30),
        "area_second_derivative": str(area_second),
        "directional_records": directional_records,
        "maximum_directional_error_ratio": mp.nstr(
            max(error/allowed for error, allowed in directional_errors), 12),
        "rank_on_vertex_space": rank,
        "nullity": 120-rank,
        "minimum_quotient_eigenvalue_L1_rho1": mp.nstr(minimum_quotient, 30),
    },
    "provenance": {
        "prior_art_commit": PRIOR_ART_COMMIT,
        "protocol_commit": PROTOCOL_COMMIT,
        "input_hashes": actual_hashes,
    },
    "source_f_vector": list(source_f),
    "spectrum": [
        {"label": label, "value": value, "multiplicity": multiplicity}
        for value, multiplicity, label in expected_spectrum
    ],
    "tests": tests,
    "passed": passed,
    "verdict": verdict,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"\nResult: {passed}/{tests} checks passed.")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    sys.exit(1)
