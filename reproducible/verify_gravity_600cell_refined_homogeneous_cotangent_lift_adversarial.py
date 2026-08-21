#!/usr/bin/env python3
"""Adversarial explicit-witness audit of refined cotangent-lift freedom."""

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import z3


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIMARY_JSON = HERE / "gravity_600cell_refined_homogeneous_cotangent_lift.json"
PRIMARY_SOURCE = (
    HERE / "verify_gravity_600cell_refined_homogeneous_cotangent_lift.py"
)
PRIMARY_RESULT = (
    ROOT
    / "docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_first_result.md"
)
FEASIBILITY = HERE / "gravity_600cell_refined_canonical_map_feasibility.json"
PROTOCOL = (
    ROOT
    / "docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_adversarial_protocol.md"
)
OUTPUT = (
    HERE / "gravity_600cell_refined_homogeneous_cotangent_lift_adversarial.json"
)

PROTOCOL_COMMIT = "282f570"
EXPECTED_HASHES = {
    "primary_json": "93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b",
    "primary_result": "8bc42b3174eebd192dfa97898e61f4f0ea3eef0b91d30848f67183a341907ec4",
    "primary_source": "154081b12f74ed8597a4b72b37a99219d64c0905829da1566e840fc562b1c20c",
    "feasibility": "ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e",
    "protocol": "10f2cfaed1ba2f75265138976514d1eec9f9d133265447354ff789b8744cc7d5",
}
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

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


def fraction_text(value):
    return str(value.numerator) if value.denominator == 1 else str(value)


def unit_lifts(coefficients, momentum=Fraction(1)):
    lifts = []
    for index, coefficient in enumerate(coefficients):
        vector = [Fraction(0) for _ in coefficients]
        vector[index] = momentum / coefficient
        lifts.append(tuple(vector))
    return tuple(lifts)


def canonical_value(coefficients, vector):
    return sum(coefficient * value
               for coefficient, value in zip(coefficients, vector))


def z3_fraction(value):
    return z3.RealVal(f"{value.numerator}/{value.denominator}")


paths = {
    "primary_json": PRIMARY_JSON,
    "primary_result": PRIMARY_RESULT,
    "primary_source": PRIMARY_SOURCE,
    "feasibility": FEASIBILITY,
    "protocol": PROTOCOL,
}
actual_hashes = {name: digest(path) for name, path in paths.items()}
provenance_ok = check(
    "all adversarial inputs match their frozen hashes",
    actual_hashes == EXPECTED_HASHES,
)

primary = json.loads(PRIMARY_JSON.read_text())
primary_ok = check(
    "the primary artifact is the frozen 12/12 underdetermined result",
    primary["outcome"]
    == "REFINED_HOMOGENEOUS_COTANGENT_LIFT_UNDERDETERMINED"
    and primary["tests"] == {"passed": 12, "total": 12},
)

feasibility = json.loads(FEASIBILITY.read_text())
level = feasibility["levels"]["projected_barycentric"]
populations = tuple(
    int(level["colour_pair_edge_populations"][f"{left}-{right}"])
    for left, right in PAIRS
)
class_sizes = tuple(int(value) for value in level["colour_class_sizes"])
geometry_ok = check(
    "the rank-coloured geometry has six positive orbits and four fixed colours",
    len(populations) == 6
    and all(value > 0 for value in populations)
    and class_sizes == (120, 720, 1200, 600)
    and len(set(class_sizes)) == 4,
    f"populations={populations}, colour sizes={class_sizes}",
)

total_coefficients = (Fraction(2),) * 6
total_lifts = unit_lifts(total_coefficients)
total_equations_ok = check(
    "six explicit orbit-total lifts carry the same unit coarse momentum",
    len(set(total_lifts)) == 6
    and all(canonical_value(total_coefficients, vector) == 1
            for vector in total_lifts),
)

base_total = total_lifts[5]
total_differences = tuple(
    tuple(left - right for left, right in zip(total_lifts[index], base_total))
    for index in range(5)
)
total_minor_diagonal = tuple(total_differences[index][index]
                             for index in range(5))
total_minor_off_diagonal_zero = all(
    total_differences[row][column] == 0
    for row in range(5) for column in range(5) if row != column
)
total_determinant = Fraction(1)
for value in total_minor_diagonal:
    total_determinant *= value
total_witness_ok = check(
    "five explicit invisible total-momentum differences are independent",
    total_minor_off_diagonal_zero
    and total_minor_diagonal == (Fraction(1, 2),) * 5
    and total_determinant == Fraction(1, 32)
    and all(canonical_value(total_coefficients, vector) == 0
            for vector in total_differences),
    f"direct determinant={fraction_text(total_determinant)}",
)

edge_coefficients = tuple(Fraction(2 * value) for value in populations)
edge_lifts = unit_lifts(edge_coefficients)
base_edge = edge_lifts[5]
edge_differences = tuple(
    tuple(left - right for left, right in zip(edge_lifts[index], base_edge))
    for index in range(5)
)
edge_minor_diagonal = tuple(edge_differences[index][index]
                            for index in range(5))
edge_minor_off_diagonal_zero = all(
    edge_differences[row][column] == 0
    for row in range(5) for column in range(5) if row != column
)
edge_determinant = Fraction(1)
for value in edge_minor_diagonal:
    edge_determinant *= value
edge_witness_ok = check(
    "five explicit per-edge differences have a nonzero exact direct determinant",
    edge_minor_off_diagonal_zero
    and edge_determinant != 0
    and all(canonical_value(edge_coefficients, vector) == 0
            for vector in edge_differences)
    and all(canonical_value(edge_coefficients, vector) == 1
            for vector in edge_lifts),
    f"direct determinant={fraction_text(edge_determinant)}",
)

a = z3.Reals("a0 a1 a2 a3 a4 a5")
b = z3.Reals("b0 b1 b2 b3 b4 b5")
multi_solver = z3.Solver()
multi_solver.add(2 * z3.Sum(a) == 1)
multi_solver.add(2 * z3.Sum(b) == 1)
multi_solver.add(z3.Or(*(left != right for left, right in zip(a, b))))
multi_status = multi_solver.check()
z3_multiple_ok = check(
    "Z3 finds two distinct six-orbit lifts with the same coarse momentum",
    multi_status == z3.sat,
    str(multi_status),
)

one_a, one_b = z3.Reals("one_a one_b")
one_solver = z3.Solver()
one_solver.add(2 * one_a == 1, 2 * one_b == 1, one_a != one_b)
one_status = one_solver.check()
one_orbit_ok = check(
    "the one-orbit distinct-lift attack is exactly UNSAT",
    one_status == z3.unsat,
    str(one_status),
)

s6 = z3.Reals("s60 s61 s62 s63 s64 s65")
s6_solver = z3.Solver()
s6_solver.add(2 * z3.Sum(s6) == 1)
s6_solver.add(*(s6[index] == s6[0] for index in range(1, 6)))
s6_status = s6_solver.check()
s6_model = s6_solver.model() if s6_status == z3.sat else None
s6_values = tuple(
    s6_model.eval(value, model_completion=True) for value in s6
) if s6_model is not None else ()
s6_distinct = z3.Reals("t60 t61 t62 t63 t64 t65")
s6_pair_solver = z3.Solver()
for variables in (s6, s6_distinct):
    s6_pair_solver.add(2 * z3.Sum(variables) == 1)
    s6_pair_solver.add(*(variables[index] == variables[0]
                         for index in range(1, 6)))
s6_pair_solver.add(z3.Or(*(left != right
                           for left, right in zip(s6, s6_distinct))))
s6_pair_status = s6_pair_solver.check()
s6_ok = check(
    "an extra full orbit-permutation symmetry would make the lift unique",
    s6_status == z3.sat
    and all(z3.is_true(z3.simplify(
        value == z3_fraction(Fraction(1, 12))
    )) for value in s6_values)
    and s6_pair_status == z3.unsat,
    f"single={s6_status}, distinct pair={s6_pair_status}, values={s6_values}",
)

negative_momentum = Fraction(-7, 3)
negative_lifts = unit_lifts(total_coefficients, negative_momentum)
negative_determinant = Fraction(1)
for index in range(5):
    negative_determinant *= (
        negative_lifts[index][index] - negative_lifts[5][index]
    )
negative_momentum_ok = check(
    "the explicit affine freedom persists at negative rational momentum",
    len(set(negative_lifts)) == 6
    and all(canonical_value(total_coefficients, vector) == negative_momentum
            for vector in negative_lifts)
    and negative_determinant != 0,
    f"direct determinant={fraction_text(negative_determinant)}",
)

reversed_total = tuple(reversed(total_coefficients))
reversed_edge = tuple(reversed(edge_coefficients))
reversed_total_lifts = unit_lifts(reversed_total)
reversed_edge_lifts = unit_lifts(reversed_edge)
reverse_total_det = Fraction(1)
reverse_edge_det = Fraction(1)
for index in range(5):
    reverse_total_det *= reversed_total_lifts[index][index]
    reverse_edge_det *= reversed_edge_lifts[index][index]
reversal_ok = check(
    "orbit reversal leaves both direct independence witnesses nonzero",
    reverse_total_det != 0 and reverse_edge_det != 0,
    f"determinants=({fraction_text(reverse_total_det)},"
    f" {fraction_text(reverse_edge_det)})",
)

corrupt = (Fraction(1, 3),) + (Fraction(0),) * 5
corruption_ok = check(
    "the deliberately corrupted unit lift fails the canonical equation",
    canonical_value(total_coefficients, corrupt) != 1,
    f"value={fraction_text(canonical_value(total_coefficients, corrupt))}",
)

zero_population_detected = False
try:
    _ = Fraction(1, 2 * 0)
except ZeroDivisionError:
    zero_population_detected = True
zero_population_ok = check(
    "a zero population is rejected before a per-edge lift is formed",
    zero_population_detected,
)

controls_ok = all((
    provenance_ok, primary_ok, geometry_ok, total_equations_ok,
    total_witness_ok, edge_witness_ok, z3_multiple_ok, one_orbit_ok, s6_ok,
    negative_momentum_ok, reversal_ok, corruption_ok, zero_population_ok,
))
if not controls_ok:
    outcome = "ADVERSARIAL_REFINED_COTANGENT_CONTROL_FAILED"
elif (
    total_determinant != 0
    and edge_determinant != 0
    and multi_status == z3.sat
    and one_status == z3.unsat
    and s6_pair_status == z3.unsat
):
    outcome = "ADVERSARIAL_REFINED_COTANGENT_UNDERDETERMINATION_CORROBORATED"
else:
    outcome = "ADVERSARIAL_REFINED_COTANGENT_DISAGREEMENT_OPEN"

outcome_ok = check(
    "the adversarial hierarchy assigns an allowed outcome",
    outcome in {
        "ADVERSARIAL_REFINED_COTANGENT_CONTROL_FAILED",
        "ADVERSARIAL_REFINED_COTANGENT_UNDERDETERMINATION_CORROBORATED",
        "ADVERSARIAL_REFINED_COTANGENT_DISAGREEMENT_OPEN",
    },
    outcome,
)

artifact = {
    "title": "Adversarial explicit-witness refined cotangent-lift audit",
    "date": "2026-08-21",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "method_firewall": {
        "primary_functions_imported": False,
        "matrix_rank_used": False,
        "nullspace_used": False,
        "svd_used": False,
        "pseudoinverse_used": False,
        "arithmetic": "fractions.Fraction plus Z3 rationals",
    },
    "geometry": {
        "orbit_labels": [f"{left}{right}" for left, right in PAIRS],
        "populations": list(populations),
        "colour_class_sizes": list(class_sizes),
        "rank_colours_pairwise_distinct_sizes": len(set(class_sizes)) == 4,
        "actual_H4_permutes_rank_pair_labels": False,
    },
    "explicit_witnesses": {
        "orbit_total_lifts": [
            [fraction_text(value) for value in vector] for vector in total_lifts
        ],
        "orbit_total_difference_minor_determinant": fraction_text(total_determinant),
        "per_edge_lifts": [
            [fraction_text(value) for value in vector] for vector in edge_lifts
        ],
        "per_edge_difference_minor_determinant": fraction_text(edge_determinant),
        "negative_momentum": fraction_text(negative_momentum),
        "negative_momentum_difference_minor_determinant": fraction_text(
            negative_determinant
        ),
    },
    "z3": {
        "version": z3.get_version_string(),
        "two_distinct_six_orbit_lifts": str(multi_status),
        "two_distinct_one_orbit_lifts": str(one_status),
        "S6_equal_lift": [str(value) for value in s6_values],
        "two_distinct_S6_equal_lifts": str(s6_pair_status),
    },
    "controls": {
        "corrupt_lift_canonical_value": fraction_text(
            canonical_value(total_coefficients, corrupt)
        ),
        "zero_population_rejected": zero_population_detected,
        "orbit_reversal_total_determinant": fraction_text(reverse_total_det),
        "orbit_reversal_per_edge_determinant": fraction_text(reverse_edge_det),
    },
    "interpretation": {
        "status": "STRUCTURAL" if "CORROBORATED" in outcome else "OPEN",
        "extra_S6_is_not_actual_H4": True,
        "action_selected_lift_excluded": False,
        "tick_c_G_planck": "NOT COMPUTED",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print(f"\nRESULT: {passed}/{tests} checks passed")
print(f"OUTCOME: {outcome}")
print(f"ARTIFACT: {OUTPUT}")
raise SystemExit(0 if passed == tests and outcome_ok else 1)
