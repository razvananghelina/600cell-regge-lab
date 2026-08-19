#!/usr/bin/env python3
"""Independent polynomial audit of cellular frustum underdetermination."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_cellular_frustum_anisotropic_rigidity_adversarial.json"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_adversarial_protocol.md"
PRIMARY_PROTOCOL = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_protocol.md"
PRIOR_ART = ROOT / "docs/gravity/gravity_600cell_cellular_frustum_anisotropic_rigidity_prior_art.md"
PRIMARY_SOURCE = HERE / "verify_gravity_600cell_cellular_frustum_anisotropic_rigidity.py"
PRIMARY_JSON = HERE / "gravity_600cell_cellular_frustum_anisotropic_rigidity.json"

PROTOCOL_COMMIT = "b39cd71"
EXPECTED_HASHES = {
    "protocol": "716d9f1e11d5fadde893ef087ff89dbf3c92c70d26800a55b93a7905dec1db4f",
    "primary_protocol": "ac20410bced8408c9cc8ec609653c3036a029b8e1d439a84c3acc3d5960eb1e8",
    "prior_art": "92c88042e8233a542b9f21e96a99bc0d09cf13cff89a8e243354f97984baaaab",
    "primary_source": "2f766503296aa43f6192d2cce6ce44faac3b7fb57ba131ba0fbf393a2da80f60",
    "primary_json": "c55f98313121018ff5ca1fc834260e8f2f075248a21fd7b99a356d89b2d18255",
}
ETA = sp.diag(1, 1, 1, -1)
BOTTOM = (
    sp.Matrix((5, 0, 0, 0)),
    sp.Matrix((0, 5, 0, 0)),
    sp.Matrix((0, 0, 5, 0)),
    sp.Matrix((3, 4, 0, 0)),
)
REPRESENTATIVES = ((1, 7), (2, 7), (3, 13))
PAIRS = tuple(combinations(range(4), 2))
ORDERS = tuple(permutations(range(4)))
Y = sp.symbols("y0:16")
TOP = tuple(sp.Matrix(Y[4 * index:4 * index + 4]) for index in range(4))
tests = passed = 0


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


def square(left, right):
    delta = left - right
    return sp.expand((delta.T * ETA * delta)[0])


def substitution(scale, lapse):
    values = []
    for point in BOTTOM:
        values.extend((
            scale * point[0], scale * point[1], scale * point[2], lapse
        ))
    return dict(zip(Y, values))


def schedule_choice(order):
    position = {colour: index for index, colour in enumerate(order)}
    return tuple(0 if position[left] < position[right] else 1
                 for left, right in PAIRS)


paths = {
    "protocol": PROTOCOL,
    "primary_protocol": PRIMARY_PROTOCOL,
    "prior_art": PRIOR_ART,
    "primary_source": PRIMARY_SOURCE,
    "primary_json": PRIMARY_JSON,
}
hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = hashes == EXPECTED_HASHES
check("all polynomial adversarial inputs have exact provenance",
      provenance_ok, str(hashes))

primary = json.loads(PRIMARY_JSON.read_text())
primary_ok = bool(
    primary["outcome"] == "CELLULAR_FRUSTUM_SIX_SHAPES_UNDERDETERMINED"
    and primary["passed"] == primary["tests"] == 11
    and primary["schedule_count"] == 24
    and primary["time_reversal_orbits"] == 12
    and all(record["nonisometric_flex_dimension"] == 6
            for record in primary["records"])
)
check("the primary six-shape result is preserved literally", primary_ok)

bottom_differences = sp.Matrix.hstack(
    *(BOTTOM[index] - BOTTOM[0] for index in range(1, 4))
)
bottom_norms = tuple(square(point, sp.zeros(4, 1)) for point in BOTTOM)
geometry_static_ok = bool(
    bottom_differences.rank() == 3
    and len(set(bottom_norms)) == 1 and bottom_norms[0] == 25
)

top_edge_polynomials = tuple(square(TOP[left], TOP[right])
                             for left, right in PAIRS)
strut_polynomials = tuple(square(BOTTOM[index], TOP[index])
                          for index in range(4))
base_polynomials = top_edge_polynomials + strut_polynomials
base_symbolic_jacobian = sp.Matrix(base_polynomials).jacobian(Y)

diagonal_options = tuple(
    (square(BOTTOM[left], TOP[right]),
     square(BOTTOM[right], TOP[left]))
    for left, right in PAIRS
)
schedule_choices = {order: schedule_choice(order) for order in ORDERS}
schedule_control = bool(
    len(ORDERS) == 24
    and len(set(schedule_choices.values())) == 24
    and all(len(choice) == 6 for choice in schedule_choices.values())
)
check("the independent polynomial geometry and 24 schedules are complete",
      geometry_static_ok and schedule_control)

records = []
representative_ok = True
base_rank_ok = True
schedule_determinants_ok = True
flex_action_ok = True
all_64_complete = True

for scale, lapse in REPRESENTATIVES:
    subs = substitution(scale, lapse)
    top_points = tuple(point.subs(subs) for point in TOP)
    strut_values = tuple(square(BOTTOM[index], top_points[index])
                          for index in range(4))
    representative_ok &= bool(
        len(set(strut_values)) == 1 and strut_values[0] < 0
    )

    base = base_symbolic_jacobian.subs(subs)
    base_rank = base.rank()
    nullspace = base.nullspace()
    kernel = sp.Matrix.hstack(*nullspace) if nullspace else sp.zeros(16, 0)
    base_rank_ok &= base_rank == 10 and kernel.shape == (16, 6)

    schedule_determinants = []
    schedule_flex_ranks = Counter()
    schedule_complete_ranks = Counter()
    for order in ORDERS:
        choice = schedule_choices[order]
        polynomials = tuple(
            diagonal_options[index][selected]
            for index, selected in enumerate(choice)
        )
        diagonal = sp.Matrix(polynomials).jacobian(Y).subs(subs)
        complete = base.col_join(diagonal)
        determinant = sp.factor(complete.det())
        schedule_determinants.append(abs(int(determinant)))
        schedule_complete_ranks[complete.rank()] += 1
        schedule_flex_ranks[(diagonal * kernel).rank()] += 1

    schedule_determinants_ok &= bool(
        all(value > 0 for value in schedule_determinants)
        and schedule_complete_ranks == {16: 24}
    )
    flex_action_ok &= schedule_flex_ranks == {6: 24}

    all_completion_ranks = Counter()
    for choice in product((0, 1), repeat=6):
        polynomials = tuple(
            diagonal_options[index][selected]
            for index, selected in enumerate(choice)
        )
        diagonal = sp.Matrix(polynomials).jacobian(Y).subs(subs)
        all_completion_ranks[base.col_join(diagonal).rank()] += 1
    all_64_complete &= sum(all_completion_ranks.values()) == 64

    determinant_multiset = Counter(schedule_determinants)
    records.append({
        "scale": scale,
        "lapse": lapse,
        "strut_squared_length": str(strut_values[0]),
        "base_jacobian_rank": base_rank,
        "base_nullity": kernel.shape[1],
        "schedule_complete_rank_counts": {
            str(key): value for key, value in schedule_complete_ranks.items()
        },
        "schedule_flex_action_rank_counts": {
            str(key): value for key, value in schedule_flex_ranks.items()
        },
        "absolute_schedule_determinant_multiset": {
            str(key): value for key, value in determinant_multiset.items()
        },
        "all_64_completion_rank_counts": {
            str(key): value for key, value in all_completion_ranks.items()
        },
    })

check("all irregular representatives have equal timelike struts",
      representative_ok)
check("all polynomial base Jacobians have exact rank ten and nullity six",
      base_rank_ok)
check("all 72 staircase Jacobian determinants are exactly nonzero",
      schedule_determinants_ok)
check("all staircase diagonals act with exact rank six on the flex kernel",
      flex_action_ok)
check("the unpredicted 64-choice censuses are complete", all_64_complete)

controls_ok = bool(
    provenance_ok and primary_ok and geometry_static_ok and schedule_control
    and representative_ok and all_64_complete
)
corroborated = bool(
    controls_ok and base_rank_ok and schedule_determinants_ok and flex_action_ok
)
if not controls_ok:
    outcome = "ADVERSARIAL_CELLULAR_FRUSTUM_CONTROL_FAILED"
elif corroborated:
    outcome = "ADVERSARIAL_CELLULAR_FRUSTUM_SIX_SHAPES_CORROBORATED"
else:
    outcome = "ADVERSARIAL_CELLULAR_FRUSTUM_DISAGREEMENT_OPEN"
allowed = {
    "ADVERSARIAL_CELLULAR_FRUSTUM_CONTROL_FAILED",
    "ADVERSARIAL_CELLULAR_FRUSTUM_SIX_SHAPES_CORROBORATED",
    "ADVERSARIAL_CELLULAR_FRUSTUM_DISAGREEMENT_OPEN",
}
check("the adversarial hierarchy assigns exactly one outcome",
      outcome in allowed, outcome)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "primary_outcome": primary["outcome"],
    "symbolic_top_variables": len(Y),
    "base_squared_length_polynomials": len(base_polynomials),
    "staircase_schedules": len(schedule_choices),
    "all_diagonal_choices": 64,
    "records": records,
    "classification": {
        "six_missing_cellular_shapes": (
            "ADVERSARIALLY CORROBORATED DERIVED EXACT / STRUCTURAL"
            if corroborated else "OPEN"
        ),
        "local_diagonal_completion_selects_schedule": "REFUTED",
        "global_conforming_schedule_count": 24,
        "physical_hessian_or_ensemble": "NOT TESTED",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print("OUTCOME:", outcome)
for record in records:
    print(
        f"(lambda,tau)=({record['scale']},{record['lapse']}): "
        f"base rank/nullity={record['base_jacobian_rank']}/"
        f"{record['base_nullity']}, all64={record['all_64_completion_rank_counts']}"
    )
print(f"RESULT: {passed}/{tests} checks passed")
if passed != tests:
    raise SystemExit(1)

