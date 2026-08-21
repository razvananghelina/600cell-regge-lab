#!/usr/bin/env python3
"""Arbitrary-precision direct-simplex adjudication of tick scale covariance.

Precision protocol commit: 093055e.
Pre-evaluation correction commit: 287dae9.
Recorded binary64 failure commit: 9951c7a.

The action is accumulated over all individual hinges after a literal loop over
all 2400 simplices. Neither prior action evaluator is called.
"""

import ast
from collections import Counter
import contextlib
from itertools import combinations
import io
import json
import re
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
COMBINATORIAL_SOURCE = HERE / "verify_gravity_global_boundary_legendre.py"
PRIMARY = HERE / "gravity_600cell_tick_scale_covariance.json"
BINARY64 = HERE / "gravity_600cell_tick_scale_covariance_adversarial.json"
OUTPUT = HERE / "gravity_600cell_tick_scale_covariance_precision.json"
PRECISION_PROTOCOL_COMMIT = "093055e"
PROTOCOL_CORRECTION_COMMIT = "287dae9"
BINARY64_FAILURE_COMMIT = "9951c7a"
DPS = 80
mp.mp.dps = DPS
ALPHAS = (mp.mpf(1), mp.mpf(3) / 5, mp.mpf(7) / 4)
SCALE_TOLERANCE = mp.mpf("1e-55")
PRIMARY_TOLERANCE = mp.mpf("1e-45")
IMAGINARY_TOLERANCE = mp.mpf("1e-60")
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


def text(value, digits=35):
    return mp.nstr(value, digits)


def relative_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


def load_combinatorics():
    """Load only the certified carrier and mappings, not its main audit."""
    tree = ast.parse(
        COMBINATORIAL_SOURCE.read_text(), filename=str(COMBINATORIAL_SOURCE)
    )
    cut = None
    for index, node in enumerate(tree.body):
        if (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "CONTROLS"
                for target in node.targets
            )
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("carrier-only cutoff not found")
    namespace = {
        "__file__": str(COMBINATORIAL_SOURCE),
        "__name__": "tick_scale_covariance_precision_carrier",
    }
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(COMBINATORIAL_SOURCE), "exec"), namespace)
    return namespace


print("Precision direct-2400-simplex tick covariance adjudication", flush=True)
carrier = load_combinatorics()
models = carrier["models"]
provenance_ok = bool(
    PRECISION_PROTOCOL_COMMIT == "093055e"
    and PROTOCOL_CORRECTION_COMMIT == "287dae9"
    and BINARY64_FAILURE_COMMIT == "9951c7a"
    and carrier["tests"] == carrier["passed"] == 7
    and set(models) == {"even", "odd"}
)
check("the corrected protocol and certified carrier provenance are intact", provenance_ok)


# Reconstruct constants without reading either covariance artifact.
M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON_3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON_3 * L0
TAU = mp.mpf("0.0102")
RHO = TAU**2
SLANT_SQUARE = L0_SQUARE - RHO
I = mp.mpc(0, 1)


def perturbed(values, modulus, center):
    return tuple(
        value * mp.exp(mp.mpf("1e-6") * ((index % modulus) - center))
        for index, value in enumerate(values)
    )


base_old = perturbed(tuple(L0_SQUARE for _ in range(30)), 7, 3)
base_variables = perturbed(
    tuple([SLANT_SQUARE] * 30 + [RHO] * 5 + [L0_SQUARE] * 30),
    5,
    2,
)
states = {
    alpha: {
        "old": tuple(alpha**2 * value for value in base_old),
        "variables": tuple(alpha**2 * value for value in base_variables),
        "mass": alpha * MASS,
    }
    for alpha in ALPHAS
}


def edge_square(model, edge, old_values, variables):
    edge = tuple(sorted(edge))
    if edge in model["all_edge_to_variable"]:
        index = model["all_edge_to_variable"][edge]
        return mp.mpf(model["all_edge_jacobian"][edge]) * variables[index]
    if edge in model["old_to_orbit"]:
        return old_values[model["old_to_orbit"][edge]]
    raise ValueError(f"edge absent from direct carrier: {edge}")


def simplex_squared(model, simplex, old_values, variables):
    squared = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    for left, right in combinations(range(5), 2):
        value = edge_square(
            model, (simplex[left], simplex[right]), old_values, variables
        )
        squared[left][right] = squared[right][left] = value
    return squared


def signed_volume_square(squared, vertices):
    vertices = list(vertices)
    dimension = len(vertices) - 1
    if dimension == 0:
        return mp.mpf(1)
    base = vertices[0]
    others = vertices[1:]
    gram = mp.matrix([
        [
            (
                squared[base][left]
                + squared[base][right]
                - squared[left][right]
            ) / 2
            for right in others
        ]
        for left in others
    ])
    return mp.det(gram) / mp.factorial(dimension) ** 2


def log_minus(value):
    scale = max(mp.mpf(1), abs(value))
    if abs(mp.im(value)) < mp.mpf("1e-65") * scale:
        real = mp.re(value)
        if real < 0:
            return mp.log(-real) - I * mp.pi
        return mp.log(real)
    return mp.log(value)


def angle_record(squared):
    gram = mp.matrix([
        [
            (
                squared[0][left]
                + squared[0][right]
                - squared[left][right]
            ) / 2
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
    leading_minors = []
    for size in range(1, 5):
        principal = mp.matrix([
            [gram[left, right] for right in range(size)]
            for left in range(size)
        ])
        leading_minors.append(mp.det(principal))
    signs = [1] + [1 if value > 0 else -1 if value < 0 else 0 for value in leading_minors]
    negative_directions = None if 0 in signs else sum(
        left != right for left, right in zip(signs, signs[1:])
    )
    angles = {}
    minimum_argument = mp.inf
    for omitted_a, omitted_b in combinations(range(5), 2):
        hinge = tuple(
            vertex for vertex in range(5)
            if vertex not in (omitted_a, omitted_b)
        )
        hinge_volume_square = signed_volume_square(squared, hinge)
        gram_derivative = mp.matrix(4, 4)
        opposite = {omitted_a, omitted_b}
        for left in range(1, 5):
            for right in range(1, 5):
                gram_derivative[left - 1, right - 1] = (
                    int({0, left} == opposite)
                    + int({0, right} == opposite)
                    - int(left != right and {left, right} == opposite)
                ) / 2
        product = inverse * gram_derivative
        volume_derivative = simplex_volume_square * sum(
            product[index, index] for index in range(4)
        )
        denominator = (
            mp.sqrt(mp.mpc(facet_volume_squares[omitted_a]))
            * mp.sqrt(mp.mpc(facet_volume_squares[omitted_b]))
        )
        cosine = 16 * volume_derivative / denominator
        sine = -mp.mpf(4) / 3 * (
            mp.sqrt(mp.mpc(hinge_volume_square))
            * mp.sqrt(mp.mpc(simplex_volume_square))
        ) / denominator
        argument = cosine + I * sine
        minimum_argument = min(minimum_argument, abs(argument))
        angles[hinge] = -I * log_minus(argument)
    return {
        "angles": angles,
        "negative_directions": negative_directions,
        "minimum_leading_minor": min(abs(value) for value in leading_minors),
        "minimum_argument": minimum_argument,
    }


def triangle_area_square(values):
    x, y, z = values
    return (2 * (x * y + x * z + y * z) - x**2 - y**2 - z**2) / 16


def direct_action(model, old_values, variables, mass):
    triangles = tuple(
        triangle for orbit in model["triangle_orbits"] for triangle in orbit
    )
    curvature = {
        triangle: mp.pi if triangle in model["boundary_triangles"] else 2 * mp.pi
        for triangle in triangles
    }
    negative_counts = Counter()
    minimum_leading_minor = mp.inf
    minimum_argument = mp.inf
    for simplex in model["slab"]:
        record = angle_record(simplex_squared(model, simplex, old_values, variables))
        negative_counts[record["negative_directions"]] += 1
        minimum_leading_minor = min(
            minimum_leading_minor, record["minimum_leading_minor"]
        )
        minimum_argument = min(minimum_argument, record["minimum_argument"])
        for local_hinge, angle in record["angles"].items():
            triangle = tuple(sorted(simplex[position] for position in local_hinge))
            curvature[triangle] += angle
    action_sum = mp.mpc(0)
    for triangle, triangle_curvature in curvature.items():
        edges = tuple(tuple(sorted(edge)) for edge in combinations(triangle, 2))
        values = tuple(
            edge_square(model, edge, old_values, variables) for edge in edges
        )
        area = mp.sqrt(mp.mpc(triangle_area_square(values)))
        action_sum += area * triangle_curvature
    gravitational = -I * action_sum
    dust = -(8 * mp.pi * mass / 5) * sum(
        mp.sqrt(value) for value in variables[30:35]
    )
    total = gravitational + dust
    return {
        "action": total,
        "negative_counts": negative_counts,
        "minimum_leading_minor": minimum_leading_minor,
        "minimum_argument": minimum_argument,
        "maximum_imaginary": abs(mp.im(total)),
        "triangle_count": len(triangles),
    }


# Complete all six independent direct actions before opening stored results.
direct = {parity: {} for parity in models}
for parity, model in models.items():
    for alpha in ALPHAS:
        print(f"  evaluating {parity}, alpha={text(alpha, 8)}", flush=True)
        state = states[alpha]
        direct[parity][alpha] = direct_action(
            model, state["old"], state["variables"], state["mass"]
        )


branch_ok = True
scale_ok = True
records = {}
for parity, evaluations in direct.items():
    base = evaluations[mp.mpf(1)]
    parity_records = {"states": {}}
    for alpha, record in evaluations.items():
        state_branch_ok = bool(
            record["negative_counts"] == Counter({1: 2400})
            and record["minimum_leading_minor"] > mp.mpf("1e-20")
            and record["minimum_argument"] > mp.mpf("1e-6")
            and record["maximum_imaginary"] < IMAGINARY_TOLERANCE
        )
        branch_ok &= state_branch_ok
        scale_error = relative_error(record["action"], alpha**2 * base["action"])
        if alpha != 1:
            scale_ok &= bool(scale_error < SCALE_TOLERANCE)
        parity_records["states"][text(alpha, 20)] = {
            "action": text(record["action"], 65),
            "negative_counts": {str(k): int(v) for k, v in record["negative_counts"].items()},
            "minimum_leading_minor": text(record["minimum_leading_minor"], 20),
            "minimum_argument": text(record["minimum_argument"], 20),
            "maximum_imaginary": text(record["maximum_imaginary"], 12),
            "triangle_count": record["triangle_count"],
            "branch_ok": state_branch_ok,
            "scale_error": text(scale_error, 15),
        }
        if alpha != 1:
            check(
                f"{parity}, alpha={text(alpha, 8)}: arbitrary-precision direct action has degree two",
                state_branch_ok and scale_error < SCALE_TOLERANCE,
                f"relative error={text(scale_error, 8)}",
            )
    records[parity] = parity_records

check(
    "all six direct states contain 2400 Lorentzian simplices on the certified branch",
    branch_ok,
)


def parse_mpc(value):
    match = re.fullmatch(r"\((\S+)\s+([+-])\s+(\S+)j\)", value.strip())
    if match is None:
        raise ValueError(f"cannot parse mpc text: {value}")
    real = mp.mpf(match.group(1))
    imaginary = mp.mpf(match.group(3))
    if match.group(2) == "-":
        imaginary = -imaginary
    return mp.mpc(real, imaginary)


# Stored results enter only after direct construction and within-parity scaling.
primary = json.loads(PRIMARY.read_text())
binary64 = json.loads(BINARY64.read_text())
artifact_controls_ok = bool(
    primary.get("outcome") == "TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED"
    and primary.get("tests") == primary.get("passed") == 12
    and binary64.get("outcome") == "TICK_SCALE_COVARIANCE_IMPLEMENTATIONS_DISAGREE"
    and binary64.get("tests") == 12
    and binary64.get("passed") == 10
)

comparison_ok = True
for parity in models:
    high = direct[parity][mp.mpf(1)]["action"]
    orbit = parse_mpc(primary["parities"][parity]["base"]["action"])
    binary_base = mp.mpc(
        binary64["parities"][parity]["base"]["action_real"],
        binary64["parities"][parity]["base"]["action_imaginary"],
    )
    orbit_error = relative_error(high, orbit)
    binary_error = relative_error(high, binary_base)
    recorded_error = mp.mpf(str(binary64["primary_action_errors"][parity]))
    reproduction_ratio = binary_error / recorded_error
    improvement = mp.inf if orbit_error == 0 else binary_error / orbit_error
    one_ok = bool(
        orbit_error < PRIMARY_TOLERANCE
        and mp.mpf("0.5") <= reproduction_ratio <= 2
        and improvement > mp.mpf("1e30")
    )
    comparison_ok &= one_ok
    records[parity]["artifact_comparison"] = {
        "high_precision_vs_primary_orbit_error": text(orbit_error, 15),
        "high_precision_vs_binary64_direct_error": text(binary_error, 15),
        "recorded_binary64_disagreement": text(recorded_error, 15),
        "binary64_error_reproduction_ratio": text(reproduction_ratio, 15),
        "precision_improvement_factor": text(improvement, 15),
        "pass": one_ok,
    }
    check(
        f"{parity}: direct arbitrary precision selects the primary orbit action and reproduces binary64 loss",
        one_ok,
        f"orbit={text(orbit_error, 6)}, binary64={text(binary_error, 6)}, improvement={text(improvement, 6)}",
    )

check("both stored outcomes have the required frozen provenance", artifact_controls_ok)

if not (provenance_ok and branch_ok and artifact_controls_ok):
    outcome = "TICK_SCALE_COVARIANCE_PRECISION_CONTROL_FAILED"
elif not (scale_ok and comparison_ok):
    outcome = "TICK_SCALE_COVARIANCE_HIGH_PRECISION_DISAGREES"
else:
    outcome = "ABSOLUTE_CLASSICAL_TICK_NO_GO_PRECISION_ADJUDICATED"

check(
    "precision adjudication corroborates the conditional absolute classical tick no-go",
    outcome == "ABSOLUTE_CLASSICAL_TICK_NO_GO_PRECISION_ADJUDICATED",
    outcome,
)

artifact = {
    "title": "Precision adjudication of direct 600-cell tick scale covariance",
    "precision_protocol_commit": PRECISION_PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "binary64_failure_commit": BINARY64_FAILURE_COMMIT,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
    "precision_digits": DPS,
    "mechanical_independence": {
        "individual_simplices": 2400,
        "individual_hinges_no_orbit_multiplicity": True,
        "primary_orbit_action_called": False,
        "binary64_full_evaluation_called": False,
        "artifacts_read_after_all_six_direct_actions": True,
        "shared": "certified carrier, edge lookup, action convention and dust normalization",
    },
    "thresholds": {
        "within_parity_scale": text(SCALE_TOLERANCE),
        "stored_primary_action": text(PRIMARY_TOLERANCE),
        "imaginary_action": text(IMAGINARY_TOLERANCE),
    },
    "parities": records,
    "interpretation": {
        "derived_exact_adversarially_corroborated": (
            "the stated scale-free classical Regge-dust action is globally scale "
            "covariant when geometrized masses scale with geometry"
        ),
        "derived_negative": (
            "it cannot select an absolute nonzero tick under those hypotheses"
        ),
        "not_excluded": [
            "tau/L",
            "tau_next/tau0",
            "relational dust time",
            "a scale fixed by independently supplied dimensionful physics",
        ],
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print(f"RESULT: {passed}/{tests}", flush=True)
print(f"OUTCOME: {outcome}", flush=True)
if passed != tests:
    raise SystemExit(1)

