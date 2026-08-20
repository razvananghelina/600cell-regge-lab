#!/usr/bin/env python3
"""All-schedule H4 stationary-fill census on P(sd K_600).

Prior-art commit: 5518fa7.
Protocol commit: 3d36c54.
No root, Hessian, spectrum or physical target is computed.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_refined_h4_stationary_fill.json"
PRIOR_ART_COMMIT = "5518fa7"
PROTOCOL_COMMIT = "3d36c54"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_gravity_global_regge_orbits.py":
        "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "reproducible/gravity_600cell_refined_canonical_map_feasibility.json":
        "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py":
        "36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32",
    "reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json":
        "2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_acceleration_blind.py":
        "496ee770ad06cf4a7f0bca79153042cca7e2821c179bdb5b27f1fdb9f393ba2b",
    "reproducible/gravity_600cell_projected_rank_edgewise_local_dust.json":
        "53463e5271301ae41eb26564875d26991ddea8024a9e09ae3c302d428ad39779",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_local_dust.py":
        "e1064380d8580e458ebfcb990285181b2ca5b092f7738a6729b0feec528df2e0",
    "reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json":
        "0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5",
}
PAIR4 = tuple(combinations(range(4), 2))
LOCAL_TRIANGLES = np.asarray(tuple(combinations(range(5), 3)), dtype=np.int8)
TAU_TEXT = "0.0102"
PRIMARY_DPS = 100
SECONDARY_DPS = 140
FD_STEP_TEXTS = ("1e-15", "5e-16")
FD_GATE_TEXT = "1e-24"
EXPECTED_F = (2640, 17040, 28800, 14400)
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def mp_text(value, digits=80):
    return mp.nstr(value, digits)


def complex_record(value, digits=80):
    return {
        "real": mp_text(mp.re(value), digits),
        "imag": mp_text(mp.im(value), digits),
        "absolute": mp_text(abs(value), digits),
    }


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(row > 0.5)) for row in adjacency]
    result = []
    for first in range(len(adjacency)):
        for second in sorted(v for v in neighbours[first] if v > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(v for v in common_two if v > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(v for v in common_three if v > third):
                    result.append((first, second, third, fourth))
    return tuple(result)


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def barycentric_chambers(coarse_top):
    coarse_cells = all_simplices(coarse_top)
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: index for index, cell in enumerate(vertex_cells)}
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flag = (
                (ordering[0],),
                tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])),
                tetrahedron,
            )
            top.append(tuple(cell_index[cell] for cell in flag))
    colours = np.asarray([len(cell)-1 for cell in vertex_cells], dtype=np.int8)
    return vertex_cells, np.asarray(top, dtype=np.int32), colours


def pack_state_rows(rows):
    result = np.zeros(len(rows), dtype=np.int64)
    for column in range(rows.shape[1]):
        result = 8*result+rows[:, column].astype(np.int64)
    return result


def unpack_states(value, length):
    values = [0]*length
    for index in range(length-1, -1, -1):
        values[index] = int(value % 8)
        value //= 8
    return tuple((item % 4, item // 4) for item in values)


def state_label(states):
    return "|".join(f"r{rank}t{layer}" for rank, layer in states)


def staircase_slab(top, colours, order):
    rank = np.empty(4, dtype=np.int8)
    for position, colour in enumerate(order):
        rank[colour] = position
    tetra_colours = colours[top]
    permutation = np.argsort(rank[tetra_colours], axis=1)
    ordered = np.take_along_axis(top, permutation, axis=1)
    vertex_count = len(colours)
    rows = []
    for pivot in range(4):
        bottom = ordered[:, :pivot+1]
        top_part = ordered[:, pivot:]+vertex_count
        rows.append(np.concatenate((bottom, top_part), axis=1))
    return np.concatenate(rows).astype(np.int32)


def schedule_combinatorics(top, colours, order):
    slab = staircase_slab(top, colours, order)
    unique_slab = np.unique(np.sort(slab, axis=1), axis=0)
    vertex_count = len(colours)
    extended_states = np.concatenate((colours, colours+4))

    simplex_states = np.sort(extended_states[slab], axis=1)
    simplex_codes = pack_state_rows(simplex_states)
    simplex_counts = Counter(map(int, simplex_codes))

    triangle_rows = slab[:, LOCAL_TRIANGLES].reshape(-1, 3)
    triangle_rows = np.sort(triangle_rows, axis=1)
    occurrence_states = np.sort(extended_states[triangle_rows], axis=1)
    occurrence_triangle_codes = pack_state_rows(occurrence_states)
    occurrence_simplex_codes = np.repeat(simplex_codes, len(LOCAL_TRIANGLES))
    occurrence_pair_codes = occurrence_simplex_codes*512+occurrence_triangle_codes
    pair_codes = np.unique(occurrence_pair_codes)
    pair_positions = np.searchsorted(pair_codes, occurrence_pair_codes)

    unique_triangles, inverse, incidence = np.unique(
        triangle_rows, axis=0, return_inverse=True, return_counts=True
    )
    unique_triangle_states = np.sort(
        extended_states[unique_triangles], axis=1
    )
    unique_triangle_codes = pack_state_rows(unique_triangle_states)

    signatures = np.zeros((len(unique_triangles), len(pair_codes)), dtype=np.int16)
    np.add.at(signatures, (inverse, pair_positions), 1)

    triangle_records = []
    mixed_types = 0
    for triangle_code in sorted(set(map(int, unique_triangle_codes))):
        indices = np.flatnonzero(unique_triangle_codes == triangle_code)
        distinct_signatures = np.unique(signatures[indices], axis=0)
        if len(distinct_signatures) != 1:
            mixed_types += 1
            continue
        signature = distinct_signatures[0]
        contributions = []
        for pair_position in np.flatnonzero(signature):
            pair_code = int(pair_codes[pair_position])
            simplex_code, pair_triangle_code = divmod(pair_code, 512)
            contributions.append({
                "simplex": unpack_states(simplex_code, 5),
                "triangle": unpack_states(pair_triangle_code, 3),
                "multiplicity": int(signature[pair_position]),
            })
        triangle_records.append({
            "states": unpack_states(triangle_code, 3),
            "count": len(indices),
            "incidence": int(incidence[indices[0]]),
            "incidence_values": sorted(set(map(int, incidence[indices]))),
            "contributions": contributions,
        })

    simplex_records = [
        {"states": unpack_states(code, 5), "count": count}
        for code, count in sorted(simplex_counts.items())
    ]
    boundary_triangle_count = sum(
        record["count"] for record in triangle_records
        if len({layer for _, layer in record["states"]}) == 1
    )
    return {
        "order": tuple(map(int, order)),
        "pentachora": len(slab),
        "distinct_pentachora": len(unique_slab),
        "triangles": len(unique_triangles),
        "boundary_triangles": boundary_triangle_count,
        "mixed_triangle_types": mixed_types,
        "simplex_types": simplex_records,
        "triangle_types": triangle_records,
    }


def signed_volume_square(squared, local_vertices):
    vertices = list(local_vertices)
    dimension = len(vertices)-1
    if dimension == 0:
        return mp.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = mp.matrix([
        [
            (squared[base][left]+squared[base][right]
             - squared[left][right])/2
            for right in others
        ]
        for left in others
    ])
    return mp.det(gram)/(mp.factorial(dimension)**2)


def log_minus(value):
    scale = max(mp.mpf(1), abs(value))
    branch_floor = mp.power(10, -mp.mp.dps+20)*scale
    if abs(mp.im(value)) < branch_floor:
        real = mp.re(value)
        if real < 0:
            return mp.log(-real)-mp.j*mp.pi
        return mp.log(real)
    return mp.log(value)


def angle_record(squared):
    gram = mp.matrix([
        [
            (squared[0][left]+squared[0][right]-squared[left][right])/2
            for right in range(1, 5)
        ]
        for left in range(1, 5)
    ])
    inverse = gram**-1
    simplex_volume_square = signed_volume_square(squared, range(5))
    facet_volume_squares = {
        omitted: signed_volume_square(
            squared, [vertex for vertex in range(5) if vertex != omitted]
        )
        for omitted in range(5)
    }
    angles = {}
    maximum_identity_residual = mp.mpf(0)
    minimum_argument = mp.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = signed_volume_square(squared, hinge)
        derivative = mp.matrix(4, 4)
        opposite = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                derivative[left-1, right-1] = (
                    int({0, left} == opposite)
                    + int({0, right} == opposite)
                    - int(left != right and {left, right} == opposite)
                )/2
        product = inverse*derivative
        volume_derivative = simplex_volume_square*sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            mp.sqrt(mp.mpc(facet_volume_squares[omitted_a]))
            * mp.sqrt(mp.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16*volume_derivative/denominator
        sine = -mp.mpf(4)/3*(
            mp.sqrt(mp.mpc(hinge_volume_square))
            * mp.sqrt(mp.mpc(simplex_volume_square))
        )/denominator
        maximum_identity_residual = max(
            maximum_identity_residual, abs(cosine*cosine+sine*sine-1)
        )
        argument = cosine+mp.j*sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge] = -mp.j*log_minus(argument)
    return angles, maximum_identity_residual, minimum_argument


def exact_geometry(dps):
    with mp.workdps(dps):
        phi = (1+mp.sqrt(5))/2
        c_value = phi/2

        def norm_square(size):
            return size*(1+(size-1)*c_value)

        unit_squares = {}
        for left, right in PAIR4:
            a, b = left+1, right+1
            dot = a*(1+(b-1)*c_value)/mp.sqrt(
                norm_square(a)*norm_square(b)
            )
            unit_squares[left, right] = 2-2*dot

        squared = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
        for left, right in PAIR4:
            squared[left][right] = squared[right][left] = unit_squares[left, right]
        gram = mp.matrix([
            [
                (squared[3][left]+squared[3][right]-squared[left][right])/2
                for right in range(3)
            ]
            for left in range(3)
        ])
        chamber_volume = mp.sqrt(mp.det(gram))/6
        total_volume = 14400*chamber_volume
        s0 = (2*mp.pi**2/total_volume)**(mp.mpf(1)/3)

        inverse = gram**-1
        normals = (
            mp.matrix([1, 0, 0]), mp.matrix([0, 1, 0]),
            mp.matrix([0, 0, 1]), mp.matrix([-1, -1, -1]),
        )

        def inner(left, right):
            return (left.T*inverse*right)[0]

        spatial_angles = {}
        for left, right in PAIR4:
            omitted_a, omitted_b = [
                value for value in range(4) if value not in (left, right)
            ]
            cosine = -inner(normals[omitted_a], normals[omitted_b])/mp.sqrt(
                inner(normals[omitted_a], normals[omitted_a])
                * inner(normals[omitted_b], normals[omitted_b])
            )
            spatial_angles[left, right] = mp.acos(cosine)

        feasibility = json.loads(
            (HERE/"gravity_600cell_refined_canonical_map_feasibility.json").read_text()
        )
        populations = feasibility["levels"]["projected_barycentric"][
            "colour_pair_edge_populations"
        ]
        edge_counts = {
            pair: int(populations[f"{pair[0]}-{pair[1]}"])
            for pair in PAIR4
        }
        curvature_bar = (
            2*mp.pi*sum(
                edge_counts[pair]*mp.sqrt(unit_squares[pair]) for pair in PAIR4
            )
            - 14400*sum(
                mp.sqrt(unit_squares[pair])*spatial_angles[pair]
                for pair in PAIR4
            )
        )
        curvature = s0*curvature_bar
        mass = curvature/(8*mp.pi)
        return {
            "unit_squares": unit_squares,
            "edge_counts": edge_counts,
            "spatial_angles": spatial_angles,
            "chamber_volume": +chamber_volume,
            "total_volume": +total_volume,
            "s0": +s0,
            "curvature": +curvature,
            "mass": +mass,
        }


def variable_keys():
    return (
        tuple(("old",)+pair for pair in PAIR4)
        + tuple(("new",)+pair for pair in PAIR4)
        + tuple(("cross",)+pair for pair in PAIR4)
        + tuple(("rho", rank) for rank in range(4))
    )


VARIABLES = variable_keys()
INTERNAL_VARIABLES = tuple(("cross",)+pair for pair in PAIR4) + tuple(
    ("rho", rank) for rank in range(4)
)


def variable_label(key):
    if key[0] == "rho":
        return f"rho_{key[1]}"
    return f"{key[0]}_{key[1]}{key[2]}"


def base_coordinates(geometry):
    tau_square = mp.mpf(TAU_TEXT)**2
    result = {}
    for pair in PAIR4:
        spatial = geometry["s0"]**2*geometry["unit_squares"][pair]
        result[("old",)+pair] = spatial
        result[("new",)+pair] = spatial
        result[("cross",)+pair] = spatial-tau_square
    for rank in range(4):
        result["rho", rank] = tau_square
    return result


def edge_coordinate(left, right):
    rank_left, layer_left = left
    rank_right, layer_right = right
    if layer_left == layer_right:
        if rank_left == rank_right:
            raise ValueError("zero boundary edge")
        pair = tuple(sorted((rank_left, rank_right)))
        return (("old" if layer_left == 0 else "new",)+pair, mp.mpf(1))
    if rank_left == rank_right:
        return (("rho", rank_left), mp.mpf(-1))
    return (("cross",)+tuple(sorted((rank_left, rank_right))), mp.mpf(1))


def simplex_squared(states, coordinates):
    squared = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        key, jacobian = edge_coordinate(states[left], states[right])
        value = jacobian*coordinates[key]
        squared[left][right] = squared[right][left] = value
    return squared


def triangle_area_and_derivatives(states, coordinates):
    edge_records = []
    for left, right in combinations(range(3), 2):
        key, jacobian = edge_coordinate(states[left], states[right])
        edge_records.append((key, jacobian, jacobian*coordinates[key]))
    x, y, z = (record[2] for record in edge_records)
    area_square = (2*(x*y+x*z+y*z)-x*x-y*y-z*z)/16
    partials = ((y+z-x)/8, (x+z-y)/8, (x+y-z)/8)
    area = mp.sqrt(mp.mpc(area_square))
    derivatives = [
        (key, partial*jacobian/(2*area))
        for (key, jacobian, _), partial in zip(edge_records, partials)
    ]
    return area, derivatives


def evaluate_schedule(combinatorics, geometry, coordinates):
    angle_lookup = {}
    maximum_identity = mp.mpf(0)
    minimum_argument = mp.inf
    for simplex_record in combinatorics["simplex_types"]:
        states = simplex_record["states"]
        squared = simplex_squared(states, coordinates)
        angles, identity, argument = angle_record(squared)
        maximum_identity = max(maximum_identity, identity)
        minimum_argument = min(minimum_argument, argument)
        for local_triangle, angle in angles.items():
            triangle_states = tuple(sorted(states[index] for index in local_triangle))
            angle_lookup[(states, triangle_states)] = angle

    gravitational_sum = mp.mpc(0)
    gradient = {key: mp.mpc(0) for key in VARIABLES}
    curvature_records = []
    for triangle_record in combinatorics["triangle_types"]:
        states = triangle_record["states"]
        boundary = len({layer for _, layer in states}) == 1
        curvature = mp.pi if boundary else 2*mp.pi
        for contribution in triangle_record["contributions"]:
            curvature += contribution["multiplicity"]*angle_lookup[
                (contribution["simplex"], contribution["triangle"])
            ]
        area, derivatives = triangle_area_and_derivatives(states, coordinates)
        multiplicity = triangle_record["count"]
        gravitational_sum += multiplicity*area*curvature
        for key, area_derivative in derivatives:
            gradient[key] += (
                -mp.j*multiplicity*curvature*area_derivative*coordinates[key]
            )
        curvature_records.append((states, curvature))

    gravitational = -mp.j*gravitational_sum
    dust = mp.mpf(0)
    rank_mass = geometry["mass"]/4
    for rank in range(4):
        rho = coordinates["rho", rank]
        dust -= 8*mp.pi*rank_mass*mp.sqrt(rho)
        gradient["rho", rank] -= 4*mp.pi*rank_mass*mp.sqrt(rho)
    action = gravitational+dust
    maximum_imaginary_curvature = max(
        abs(mp.im(value)) for _, value in curvature_records
    )
    return {
        "action": action,
        "gravitational": gravitational,
        "dust": dust,
        "gradient": gradient,
        "maximum_angle_identity_residual": maximum_identity,
        "minimum_angle_argument": minimum_argument,
        "maximum_imaginary_curvature": maximum_imaginary_curvature,
    }


def finite_difference_control(combinatorics, geometry, coordinates):
    analytic = evaluate_schedule(combinatorics, geometry, coordinates)
    relative_errors = []
    records = {}
    for key in VARIABLES:
        derivatives = []
        for step_text in FD_STEP_TEXTS:
            step = mp.mpf(step_text)
            plus = dict(coordinates)
            minus = dict(coordinates)
            plus[key] *= mp.exp(step)
            minus[key] *= mp.exp(-step)
            value_plus = evaluate_schedule(combinatorics, geometry, plus)["action"]
            value_minus = evaluate_schedule(combinatorics, geometry, minus)["action"]
            derivatives.append((value_plus-value_minus)/(2*step))
        richardson = (4*derivatives[1]-derivatives[0])/3
        expected = analytic["gradient"][key]
        relative = abs(richardson-expected)/max(mp.mpf(1), abs(expected))
        relative_errors.append(relative)
        records[variable_label(key)] = {
            "analytic": complex_record(expected, 55),
            "richardson": complex_record(richardson, 55),
            "relative_error": mp_text(relative, 30),
        }
    return max(relative_errors), records


def internal_per_edge(result, geometry):
    per_edge = {}
    total = {}
    for key in INTERNAL_VARIABLES:
        value = result["gradient"][key]
        if key[0] == "cross":
            count = geometry["edge_counts"][key[1], key[2]]
        else:
            count = (120, 720, 1200, 600)[key[1]]
        total[key] = value
        per_edge[key] = value/count
    return total, per_edge


def induced_lapse_derivative(result, coordinates):
    value = sum(result["gradient"]["rho", rank] for rank in range(4))
    rho = mp.mpf(TAU_TEXT)**2
    for pair in PAIR4:
        value -= rho/coordinates[("cross",)+pair]*result["gradient"][("cross",)+pair]
    return value


print("="*78)
print("REFINED H4 STATIONARY-FILL CENSUS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "all frozen action, carrier and dust inputs have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "5518fa7" and PROTOCOL_COMMIT == "3d36c54",
    str(actual_hashes),
)

feasibility = json.loads(
    (HERE/"gravity_600cell_refined_canonical_map_feasibility.json").read_text()
)
acceleration = json.loads(
    (HERE/"gravity_600cell_projected_rank_edgewise_acceleration_blind.json").read_text()
)
local_dust = json.loads(
    (HERE/"gravity_600cell_projected_rank_edgewise_local_dust.json").read_text()
)
balanced = json.loads(
    (HERE/"gravity_600cell_projected_rank_edgewise_balanced_slab.json").read_text()
)
upstream_ok = check(
    "the frozen artifacts carry the required accepted but scoped outcomes",
    feasibility["outcome"] == "REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED"
    and feasibility["tests"] == {"passed": 8, "total": 8}
    and acceleration["outcome"]
        == "CANONICAL_CARRIER_ACCELERATION_COEFFICIENTS_DERIVED"
    and acceleration["passed"] == acceleration["tests"] == 8
    and local_dust["outcome"] == "P1_LOCAL_DUST_WEIGHTS_DERIVED_CONDITIONALLY"
    and local_dust["passed"] == local_dust["tests"] == 11
    and balanced["selection"]["ordered_slab_alternatives"] == 24
    and balanced["selection"]["existence_passes"],
)

vertices, adjacency, _ = build_600cell()
coarse_top = tetrahedra_from_adjacency(adjacency)
vertex_cells, top, colours = barycentric_chambers(coarse_top)
spatial_cells = all_simplices(tuple(map(tuple, top)))
topology_ok = check(
    "K0 is reconstructed with its exact flag counts and proper rank colouring",
    tuple(len(layer) for layer in spatial_cells) == EXPECTED_F
    and len(vertex_cells) == 2640 and len(top) == 14400
    and all({int(colours[vertex]) for vertex in tetrahedron} == set(range(4))
            for tetrahedron in top),
    f"f={tuple(len(layer) for layer in spatial_cells)}",
)

geometry100 = exact_geometry(PRIMARY_DPS)
geometry140 = exact_geometry(SECONDARY_DPS)
with mp.workdps(SECONDARY_DPS):
    frozen_volume = mp.mpf(str(
        acceleration["levels"]["projected_barycentric"]["volume_bar"]
    ))
    frozen_mass = mp.mpf(str(
        acceleration["levels"]["projected_barycentric"][
            "selected_total_dust_mass"
        ]
    ))
    volume_error = abs(geometry140["total_volume"]-frozen_volume)
    mass_error = abs(geometry140["mass"]-frozen_mass)
geometry_ok = check(
    "the exact rank geometry reproduces the frozen volume and selected mass",
    volume_error < mp.mpf("5e-13") and mass_error < mp.mpf("5e-12")
    and max(abs(geometry100["unit_squares"][pair]
                - geometry140["unit_squares"][pair]) for pair in PAIR4)
        < mp.mpf("1e-90"),
    f"V={mp_text(geometry140['total_volume'], 30)}, "
    f"M={mp_text(geometry140['mass'], 30)}, "
    f"errors=({mp_text(volume_error, 5)},{mp_text(mass_error, 5)})",
)

dust_orbits = {
    int(record["size"]): mp.mpf(str(record["weight_mean"]))
    for record in local_dust["levels"]["projected_barycentric"][
        "symmetry_orbits"
    ]
}
rank_sizes = (120, 720, 1200, 600)
with mp.workdps(SECONDARY_DPS):
    p1_weight_errors = []
    for size in rank_sizes:
        expected_weight = geometry140["total_volume"]/(4*size)
        p1_weight_errors.append(abs(expected_weight-dust_orbits[size]))
p1_ok = check(
    "the conditional P1 weights put exactly one quarter of mass in each rank",
    max(p1_weight_errors) < mp.mpf("2e-12"),
    f"max weight error={mp_text(max(p1_weight_errors), 6)}",
)

orders = tuple(permutations(range(4)))
combinatorics = []
combinatorics_ok = True
print("[INFO] constructing exact triangle signatures for all 24 schedules", flush=True)
for index, order in enumerate(orders):
    record = schedule_combinatorics(top, colours, order)
    combinatorics.append(record)
    ok = (
        record["pentachora"] == record["distinct_pentachora"] == 57600
        and record["boundary_triangles"] == 57600
        and record["mixed_triangle_types"] == 0
        and len(record["simplex_types"]) == 4
        and all(item["count"] == 14400 for item in record["simplex_types"])
        and all(len(item["incidence_values"]) == 1
                for item in record["triangle_types"])
    )
    combinatorics_ok &= ok
    if index in (0, 5, 11, 17, 23):
        print(f"[INFO] schedules completed: {index+1}/24", flush=True)
combinatorics_gate = check(
    "all 24 slabs have unmixed exact H4 triangle-incidence signatures",
    combinatorics_ok,
    str({str(record["order"]): {
        "triangles": record["triangles"],
        "triangle_types": len(record["triangle_types"]),
    } for record in combinatorics}),
)


def precision_census(dps, geometry):
    with mp.workdps(dps):
        coordinates = base_coordinates(geometry)
        evaluations = [
            evaluate_schedule(record, geometry, coordinates)
            for record in combinatorics
        ]
        internal = [internal_per_edge(value, geometry) for value in evaluations]
        lapse = [induced_lapse_derivative(value, coordinates)
                 for value in evaluations]
        return coordinates, evaluations, internal, lapse


print("[INFO] evaluating all schedules at 100 and 140 decimal digits", flush=True)
coordinates100, evaluations100, internal100, lapse100 = precision_census(
    PRIMARY_DPS, geometry100
)
coordinates140, evaluations140, internal140, lapse140 = precision_census(
    SECONDARY_DPS, geometry140
)

with mp.workdps(SECONDARY_DPS):
    max_identity = max(
        value["maximum_angle_identity_residual"] for value in evaluations140
    )
    minimum_argument = min(value["minimum_angle_argument"]
                           for value in evaluations140)
    maximum_imaginary_action = max(abs(mp.im(value["action"]))
                                   for value in evaluations140)
branch_ok = check(
    "every corrected Lorentzian angle and action stays on the frozen branch",
    max_identity < mp.mpf("1e-100")
    and minimum_argument > mp.mpf("1e-20")
    and maximum_imaginary_action < mp.mpf("1e-90"),
    f"identity={mp_text(max_identity, 6)}, min|arg|={mp_text(minimum_argument, 8)}, "
    f"max Im(S)={mp_text(maximum_imaginary_action, 6)}",
)

print("[INFO] finite-difference controls on first and reverse schedules", flush=True)
with mp.workdps(PRIMARY_DPS):
    fd_first, fd_first_records = finite_difference_control(
        combinatorics[0], geometry100, coordinates100
    )
    fd_reverse, fd_reverse_records = finite_difference_control(
        combinatorics[-1], geometry100, coordinates100
    )
fd_max = max(fd_first, fd_reverse)
fd_ok = check(
    "independent Richardson differences reproduce all 44 analytic derivatives",
    fd_max < mp.mpf(FD_GATE_TEXT),
    f"max relative error={mp_text(fd_max, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    precision_differences = []
    all_per_edge140 = []
    for (_, per100), (_, per140) in zip(internal100, internal140):
        for key in INTERNAL_VARIABLES:
            precision_differences.append(abs(per100[key]-per140[key]))
            all_per_edge140.append(abs(per140[key]))
    scale = max([mp.mpf(1)]+all_per_edge140)
    maximum_precision_difference = max(precision_differences)
    epsilon = 100*maximum_precision_difference+mp.mpf("1e-60")*scale

    action_precision_difference = max(
        abs(left["action"]-right["action"])
        for left, right in zip(evaluations100, evaluations140)
    )
    action_scale = max([mp.mpf(1)]+[abs(item["action"])
                                    for item in evaluations140])
    action_epsilon = (
        100*action_precision_difference+mp.mpf("1e-60")*action_scale
    )
    action_spread = max(
        abs(item["action"]-evaluations140[0]["action"])
        for item in evaluations140
    )
action_ok = check(
    "all 24 total actions agree at the common induced flat fill",
    action_spread <= action_epsilon,
    f"spread={mp_text(action_spread, 8)}, epsilon={mp_text(action_epsilon, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    reverse_differences = []
    order_index = {order: index for index, order in enumerate(orders)}
    for index, order in enumerate(orders):
        reverse_index = order_index[tuple(reversed(order))]
        left = internal140[index][1]
        right = internal140[reverse_index][1]
        reverse_differences.extend(
            abs(left[key]-right[key]) for key in INTERNAL_VARIABLES
        )
    maximum_reverse_difference = max(reverse_differences)
time_reversal_ok = check(
    "time reversal pairs all equal-boundary internal residual vectors",
    maximum_reverse_difference <= epsilon,
    f"max difference={mp_text(maximum_reverse_difference, 8)}, "
    f"epsilon={mp_text(epsilon, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    lapse_precision = max(abs(left-right) for left, right in zip(lapse100, lapse140))
    lapse_scale = max([mp.mpf(1)]+[abs(value) for value in lapse140])
    lapse_epsilon = 100*lapse_precision+mp.mpf("1e-60")*lapse_scale
    maximum_lapse = max(abs(value) for value in lapse140)
lapse_ok = check(
    "the induced common-lapse direction reproduces the frozen static equation",
    maximum_lapse <= lapse_epsilon,
    f"max={mp_text(maximum_lapse, 8)}, epsilon={mp_text(lapse_epsilon, 8)}",
)

with mp.workdps(SECONDARY_DPS):
    nonzero = []
    compatible = []
    vector_classes = []
    schedule_records = []
    for index, order in enumerate(orders):
        total, per_edge = internal140[index]
        vector = tuple(per_edge[key] for key in INTERNAL_VARIABLES)
        class_index = None
        for candidate_index, representative in enumerate(vector_classes):
            if max(abs(left-right) for left, right in zip(vector, representative)) <= epsilon:
                class_index = candidate_index
                break
        if class_index is None:
            class_index = len(vector_classes)
            vector_classes.append(vector)
        entries = {}
        for key in INTERNAL_VARIABLES:
            label = variable_label(key)
            certified = abs(per_edge[key]) > epsilon
            target = nonzero if certified else compatible
            target.append((index, key))
            entries[label] = {
                "total_log_residual": complex_record(total[key], 75),
                "per_edge_log_residual": complex_record(per_edge[key], 75),
                "certified_nonzero": certified,
            }
        schedule_records.append({
            "order": list(order),
            "residual_class": class_index,
            "action": complex_record(evaluations140[index]["action"], 75),
            "induced_lapse_derivative": complex_record(lapse140[index], 75),
            "internal": entries,
        })

    cross_nonzero = sum(key[0] == "cross" for _, key in nonzero)
    vertical_nonzero = sum(key[0] == "rho" for _, key in nonzero)
    minimum_nonzero = min(
        (abs(internal140[index][1][key]) for index, key in nonzero),
        default=mp.mpf(0),
    )
    maximum_residual = max(all_per_edge140)

controls_ok = all((
    provenance_ok, upstream_ok, topology_ok, geometry_ok, p1_ok,
    combinatorics_gate, branch_ok, fd_ok, action_ok, time_reversal_ok,
    lapse_ok,
))
if not controls_ok:
    outcome = "REFINED_H4_STATIONARY_FILL_CONTROL_FAILED"
elif nonzero:
    outcome = "REFINED_H4_INDUCED_FILL_OFF_SHELL"
else:
    outcome = "REFINED_H4_INDUCED_FILL_STATIONARY_CANDIDATE"
outcome_ok = check(
    "the frozen hierarchy assigns the H4 stationary-fill outcome",
    outcome in {
        "REFINED_H4_STATIONARY_FILL_CONTROL_FAILED",
        "REFINED_H4_INDUCED_FILL_OFF_SHELL",
        "REFINED_H4_INDUCED_FILL_STATIONARY_CANDIDATE",
    },
    outcome,
)

combinatorial_artifact = []
for record in combinatorics:
    combinatorial_artifact.append({
        "order": list(record["order"]),
        "pentachora": record["pentachora"],
        "distinct_pentachora": record["distinct_pentachora"],
        "triangles": record["triangles"],
        "boundary_triangles": record["boundary_triangles"],
        "mixed_triangle_types": record["mixed_triangle_types"],
        "simplex_types": [
            {"states": state_label(item["states"]), "count": item["count"]}
            for item in record["simplex_types"]
        ],
        "triangle_types": [
            {
                "states": state_label(item["states"]),
                "count": item["count"],
                "incidence_values": item["incidence_values"],
                "signature": [
                    {
                        "simplex": state_label(value["simplex"]),
                        "triangle": state_label(value["triangle"]),
                        "multiplicity": value["multiplicity"],
                    }
                    for value in item["contributions"]
                ],
            }
            for item in record["triangle_types"]
        ],
    })

artifact = {
    "title": "All-schedule H4 stationary-fill census on P(sd K_600)",
    "date": "2026-08-20",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "definitions": {
        "carrier": "K0=P(sd K_600)",
        "tau0_supplied": TAU_TEXT,
        "variables": [variable_label(key) for key in VARIABLES],
        "internal_variables": [variable_label(key) for key in INTERNAL_VARIABLES],
        "primary_decimal_digits": PRIMARY_DPS,
        "secondary_decimal_digits": SECONDARY_DPS,
        "finite_difference_steps": ["1e-15", "5e-16"],
        "root_hessian_spectrum_or_target_loaded": False,
    },
    "exact_geometry": {
        "unit_edge_squares": {
            f"{left}-{right}": mp_text(geometry140["unit_squares"][left, right], 75)
            for left, right in PAIR4
        },
        "edge_type_counts": {
            f"{left}-{right}": geometry140["edge_counts"][left, right]
            for left, right in PAIR4
        },
        "chamber_volume": mp_text(geometry140["chamber_volume"], 75),
        "total_chordal_volume": mp_text(geometry140["total_volume"], 75),
        "scale_s0": mp_text(geometry140["s0"], 75),
        "spatial_regge_curvature": mp_text(geometry140["curvature"], 75),
        "selected_total_mass": mp_text(geometry140["mass"], 75),
        "rank_mass": mp_text(geometry140["mass"]/4, 75),
        "frozen_volume_absolute_error": mp_text(volume_error, 30),
        "frozen_mass_absolute_error": mp_text(mass_error, 30),
    },
    "combinatorics": combinatorial_artifact,
    "precision": {
        "maximum_per_edge_100_vs_140_difference": mp_text(
            maximum_precision_difference, 40
        ),
        "per_edge_zero_envelope": mp_text(epsilon, 40),
        "maximum_angle_identity_residual": mp_text(max_identity, 40),
        "minimum_angle_argument": mp_text(minimum_argument, 40),
        "maximum_action_imaginary": mp_text(maximum_imaginary_action, 40),
        "maximum_finite_difference_relative_error": mp_text(fd_max, 40),
        "action_spread": mp_text(action_spread, 40),
        "action_envelope": mp_text(action_epsilon, 40),
        "maximum_time_reversal_residual_difference": mp_text(
            maximum_reverse_difference, 40
        ),
        "maximum_induced_lapse_residual": mp_text(maximum_lapse, 40),
        "induced_lapse_envelope": mp_text(lapse_epsilon, 40),
    },
    "finite_difference_controls": {
        "first_order": fd_first_records,
        "reverse_order": fd_reverse_records,
    },
    "census": {
        "schedule_count": len(orders),
        "internal_entries": len(orders)*len(INTERNAL_VARIABLES),
        "certified_nonzero_entries": len(nonzero),
        "zero_compatible_entries": len(compatible),
        "certified_nonzero_cross_entries": cross_nonzero,
        "certified_nonzero_vertical_entries": vertical_nonzero,
        "distinct_residual_vectors": len(vector_classes),
        "minimum_certified_nonzero_per_edge_absolute": mp_text(minimum_nonzero, 40),
        "maximum_per_edge_absolute": mp_text(maximum_residual, 40),
        "schedules": schedule_records,
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("-"*78)
print(f"NONZERO: {len(nonzero)}/240 "
      f"(cross={cross_nonzero}, vertical={vertical_nonzero})")
print(f"DISTINCT RESIDUAL VECTORS: {len(vector_classes)}")
print(f"OUTCOME: {outcome}")
print(f"RESULT: {passed}/{tests} checks passed")
sys.exit(0 if passed == tests and outcome_ok and controls_ok else 1)
