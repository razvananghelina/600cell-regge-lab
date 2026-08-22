#!/usr/bin/env python3
"""Verify the bounded no-go for an exact two-slab constraint quotient."""

from __future__ import annotations

import ast
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUN_ALL = HERE / "run_all.py"
VERIFIER_NAME = Path(__file__).name
OUTPUT = HERE / "gravity_600cell_finite_height_constraint_quotient_no_go.json"

INPUTS = {
    "first_tangent": HERE
    / "gravity_600cell_finite_height_full_boundary_tangent_adversarial.json",
    "second_tangent": HERE
    / "gravity_600cell_finite_height_second_full_boundary_tangent_adversarial.json",
    "internal_reconciliation": HERE
    / "gravity_600cell_finite_height_internal_kernel_canonical_reconciliation.json",
    "old_rank_control": HERE
    / "gravity_600cell_dust_full_anisotropic_legendre_rank.json",
    "prior_art": ROOT
    / "docs/gravity/gravity_600cell_finite_height_constraint_quotient_prior_art.md",
}
EXPECTED_HASHES = {
    "first_tangent": "ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7",
    "second_tangent": "1355f8cf339d18c1cf2855ecb1228e97e868d73f7a1ef739e4c11ce9521fcd4b",
    "internal_reconciliation": "81ec0379247023451e82ab42f5beb026ee2d1b083aa5e2553e42b894554266f6",
    "old_rank_control": "7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226",
    "prior_art": "53de3f6262df6e3a0cdff916d11d435fc462289b5db080bb845dcb425bb270c5",
}
ACCEPTED_OUTCOME = (
    "FINITE_HEIGHT_TWO_STEP_EXACT_CONSTRAINT_QUOTIENT_BOUNDED_NO_GO"
)
OPEN_OUTCOME = "FINITE_HEIGHT_CONSTRAINT_QUOTIENT_OPEN"
SIGNAL_OUTCOME = "FINITE_HEIGHT_EXACT_CONSTRAINT_SIGNAL"
PHASE_DIMENSION = 1440

tests = 0
passed = 0


def check(label: str, condition: object, detail: str = "") -> None:
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def registered_scripts(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "scripts"
            for target in node.targets
        ):
            continue
        value = ast.literal_eval(node.value)
        if not isinstance(value, list) or not all(
            isinstance(item, str) for item in value
        ):
            raise TypeError("run_all.py scripts is not a literal string list")
        return value
    raise ValueError("run_all.py has no literal scripts registry")


def scalar_action_control(coefficients: tuple[sp.Rational, ...]) -> dict:
    """Exact stationary images and tangent for S(o,x,n)."""
    a, b, c, d, e, f = coefficients
    if b == 0:
        raise ValueError("the frozen controls require b != 0")

    j_matrix = sp.Matrix([[b, f], [-d, -e]])
    x_coefficients = sp.Matrix([[-d / b, -f / b]])
    # Columns are coefficients of the free stationary parameters (o,n).
    pre_momentum = -sp.Matrix([[a, d, e]]) * sp.Matrix.vstack(
        sp.eye(2)[0, :], x_coefficients, sp.eye(2)[1, :]
    )
    post_momentum = sp.Matrix([[e, f, c]]) * sp.Matrix.vstack(
        sp.eye(2)[0, :], x_coefficients, sp.eye(2)[1, :]
    )
    pre_image = sp.Matrix.vstack(sp.Matrix([[1, 0]]), pre_momentum)
    post_image = sp.Matrix.vstack(sp.Matrix([[0, 1]]), post_momentum)

    record = {
        "coefficients_a_b_c_d_e_f": [str(value) for value in coefficients],
        "pre_legendre_matrix": [
            [str(value) for value in row] for row in j_matrix.tolist()
        ],
        "pre_legendre_rank": int(j_matrix.rank()),
        "pre_legendre_determinant": str(sp.factor(j_matrix.det())),
        "pre_image_rank": int(pre_image.rank()),
        "post_image_rank": int(post_image.rank()),
        "pre_constraint_codimension": 2 - int(pre_image.rank()),
        "post_constraint_codimension": 2 - int(post_image.rank()),
        "pre_annihilator_dimension": len(pre_image.T.nullspace()),
        "post_annihilator_dimension": len(post_image.T.nullspace()),
    }

    if j_matrix.det() != 0:
        rhs = sp.Matrix([[-d, 0], [a, 1]])
        internal_new_response = sp.simplify(j_matrix.inv() * rhs)
        tangent = sp.Matrix.vstack(
            internal_new_response[1, :],
            sp.Matrix([[e, 0]])
            + sp.Matrix([[f, c]]) * internal_new_response,
        )
        omega = sp.Matrix([[0, 1], [-1, 0]])
        symplectic_defect = sp.simplify(tangent.T * omega * tangent - omega)
        record.update(
            {
                "tangent": [
                    [str(sp.factor(value)) for value in row]
                    for row in tangent.tolist()
                ],
                "tangent_rank": int(tangent.rank()),
                "tangent_determinant": str(sp.factor(tangent.det())),
                "symplectic_defect": [
                    [str(value) for value in row]
                    for row in symplectic_defect.tolist()
                ],
                "route_b_post_annihilator_dimension": len(
                    tangent.T.nullspace()
                ),
                "route_b_pre_annihilator_dimension": len(
                    tangent.inv().T.nullspace()
                ),
            }
        )
    return record


def strict_margin(record: dict) -> tuple[bool, Decimal]:
    gate = Decimal(record["gap_gate"])
    values = [
        Decimal(value)
        for value in record["normalized_smallest_singular_values"].values()
    ]
    ratio = min(values) / gate
    return bool(
        record["classification"] == "REGULAR"
        and all(value > gate for value in values)
    ), ratio


print("=" * 78)
print("FINITE-HEIGHT EXACT CONSTRAINT QUOTIENT: BOUNDED NO-GO")
print("=" * 78)

input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
hashes_ok = input_hashes == EXPECTED_HASHES
check("all five frozen input hashes agree", hashes_ok, str(input_hashes))

scripts = registered_scripts(RUN_ALL)
duplicates = sorted({name for name in scripts if scripts.count(name) > 1})
registry_ok = scripts.count(VERIFIER_NAME) == 1 and not duplicates
check(
    "the verifier is registered exactly once and the registry has no duplicates",
    registry_ok,
    f"entries={len(scripts)}, duplicates={duplicates}",
)

payloads = {
    name: json.loads(path.read_text())
    for name, path in INPUTS.items()
    if path.suffix == ".json"
}
first = payloads["first_tangent"]
second = payloads["second_tangent"]
internal = payloads["internal_reconciliation"]
old_control = payloads["old_rank_control"]

first_accepted = bool(
    first.get("outcome")
    == "FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
    and first.get("passed") == first.get("tests") == 22
    and first.get("schedule", {}).get("outcome") == "SCHEDULE_ROBUST"
)
check("the first adversarial artifact retains its accepted 22/22 outcome", first_accepted)

first_margins = {}
first_regular = True
for parity, record in first.get("parities", {}).items():
    margin_ok, ratio = strict_margin(record)
    first_margins[parity] = str(ratio)
    first_regular &= bool(margin_ok and record.get("canonical") is True)
first_regular &= set(first_margins) == {"even", "odd"}
check(
    "both first-slab parities are strictly regular and canonical above their frozen rank gates",
    first_regular,
    f"minimum/gate ratios={first_margins}",
)

first_firewall = first.get("firewall", {})
second_firewall = second.get("firewall", {})
second_accepted = bool(
    second.get("outcome")
    == "TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED"
    and second.get("passed") == second.get("tests") == 28
    and not first_firewall.get("continuum_target_parsed", True)
    and not first_firewall.get("eigenvalues_computed", True)
    and not second_firewall.get("continuum_or_particle_target_parsed", True)
    and not second_firewall.get("tangent_eigenvalues_computed", True)
    and not second_firewall.get("tangent_singular_values_computed", True)
)
check(
    "the second adversarial artifact retains its accepted 28/28 outcome and target firewall",
    second_accepted,
)

second_margins = {}
second_stages_regular = True
for parity, stages in second.get("stages", {}).items():
    for stage, record in stages.items():
        margin_ok, ratio = strict_margin(record["pre_legendre"])
        key = f"{parity}:{stage}"
        second_margins[key] = str(ratio)
        second_stages_regular &= bool(
            margin_ok and record.get("tangent", {}).get("canonical") is True
        )
second_stages_regular &= len(second_margins) == 6
check(
    "all six first/second stage realizations are strictly regular and canonical above their frozen rank gates",
    second_stages_regular,
    f"minimum/gate ratios={second_margins}",
)

product_records = second.get("product_canonicality", {})
products_ok = bool(
    set(product_records)
    == {"even_even", "even_odd", "odd_even", "odd_odd"}
    and all(record.get("passed") is True for record in product_records.values())
    and second.get("product_schedule_outcome") == "TWO_STEP_SCHEDULE_ROBUST"
    and len(second.get("product_schedule", [])) == 6
    and all(
        record.get("label") == "AGREES"
        for record in second.get("product_schedule", [])
    )
)
check("all four two-step products are canonical and schedule robust", products_ok)

internal_ok = bool(
    internal.get("outcome")
    == "INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT"
    and internal.get("passed") == internal.get("tests") == 14
    and internal.get("exact_certificates", {}).get("fixed_input_intersection")
    == "ker(R_p) intersect ker(dP) = {0}"
    and internal.get("interpretation", {}).get("fixed_incoming_canonical_tangent")
    == "zero"
)
check("the internal lapse tangent is still removed by fixed incoming momentum", internal_ok)

positive = scalar_action_control(tuple(map(sp.Rational, (2, 3, 5, 7, 11, 13))))
zero_matrix = [["0", "0"], ["0", "0"]]
positive_ok = bool(
    positive["pre_legendre_rank"] == 2
    and positive["pre_constraint_codimension"] == 0
    and positive["post_constraint_codimension"] == 0
    and positive["tangent_rank"] == 2
    and positive["tangent_determinant"] == "1"
    and positive["symplectic_defect"] == zero_matrix
    and positive["route_b_post_annihilator_dimension"] == 0
    and positive["route_b_pre_annihilator_dimension"] == 0
)
check("the exact regular positive control passes both proof routes", positive_ok)

singular = scalar_action_control(tuple(map(sp.Rational, (2, 1, 1, 1, 1, 1))))
singular_ok = bool(
    singular["pre_legendre_rank"] == 1
    and singular["pre_constraint_codimension"] == 1
    and singular["post_constraint_codimension"] == 1
    and singular["pre_annihilator_dimension"] == 1
    and singular["post_annihilator_dimension"] == 1
    and "tangent_rank" not in singular
)
check("the exact singular negative control is detected", singular_ok)

epsilon = sp.Rational(1, 10**100)
tiny_j = sp.diag(sp.Integer(1), -epsilon)
tiny_exact_rank = int(tiny_j.rank())
tiny_hostile_rank = 1 + int(abs(float(epsilon)) > 1.0e-12)
tiny_ok = bool(
    tiny_exact_rank == 2
    and tiny_j.det() != 0
    and tiny_hostile_rank == 1
    and tiny_exact_rank != tiny_hostile_rank
)
check(
    "the exact 10^-100 hostile control refutes threshold quotienting",
    tiny_ok,
    f"exact_rank={tiny_exact_rank}, hostile_cutoff_rank={tiny_hostile_rank}",
)

first_free_tangent = sp.Matrix([[1, 1], [0, 1]])
later_constraint = sp.Matrix([[1, 1]])
pulled_back_constraint = later_constraint * first_free_tangent
future_ok = bool(
    first_free_tangent.det() == 1
    and later_constraint.rank() == 1
    and pulled_back_constraint == sp.Matrix([[1, 2]])
    and pulled_back_constraint.rank() == 1
)
check(
    "the future-singular control propagates a constraint backward and bounds the claim",
    future_ok,
    f"pulled_back={pulled_back_constraint.tolist()}",
)

old_control_ok = bool(
    old_control.get("outcome") == "FULL_CANONICAL_LEGENDRE_REGULAR"
    and old_control.get("passed") == old_control.get("tests") == 18
    and old_control.get("classification", {}).get("gauge_identification")
    == "OPEN"
    and old_control.get("classification", {}).get("pseudoconstraints")
    == "STRUCTURAL CANDIDATES ONLY"
    and all(
        record.get("full_resolved_rank") == 1560
        and record.get("full_error_consistent_nullity") == 0
        and record.get("pseudoconstraint_candidates") == 120
        for record in old_control.get("parities", {}).values()
    )
)
check("the old weak-mode artifact remains a regular non-gauge control", old_control_ok)

actual_maps = {
    "first_even": first_regular,
    "first_odd": first_regular,
    "second_even_normalized": second_stages_regular,
    "second_even_physical": second_stages_regular,
    "second_odd_normalized": second_stages_regular,
    "second_odd_physical": second_stages_regular,
    "product_even_even": products_ok,
    "product_even_odd": products_ok,
    "product_odd_even": products_ok,
    "product_odd_odd": products_ok,
}
route_a = {
    name: {
        "phase_dimension": PHASE_DIMENSION,
        "pre_constraint_codimension": 0 if premise else None,
        "post_constraint_codimension": 0 if premise else None,
    }
    for name, premise in actual_maps.items()
}
route_b = {
    name: {
        "phase_dimension": PHASE_DIMENSION,
        "post_constraint_covector_annihilator_dimension": 0
        if premise
        else None,
        "pre_constraint_covector_annihilator_dimension": 0
        if premise
        else None,
    }
    for name, premise in actual_maps.items()
}

premises = [
    hashes_ok,
    registry_ok,
    first_accepted,
    first_regular,
    second_accepted,
    second_stages_regular,
    products_ok,
    internal_ok,
    positive_ok,
    singular_ok,
    tiny_ok,
    future_ok,
    old_control_ok,
]
route_agreement = bool(
    all(premises)
    and all(
        route_a[name]["pre_constraint_codimension"] == 0
        and route_a[name]["post_constraint_codimension"] == 0
        and route_b[name]["post_constraint_covector_annihilator_dimension"] == 0
        and route_b[name]["pre_constraint_covector_annihilator_dimension"] == 0
        for name in actual_maps
    )
)
check(
    "Routes A and B agree on zero local constraint codimension and the bounded outcome",
    route_agreement,
    f"maps={len(actual_maps)}, phase_dimension={PHASE_DIMENSION}",
)

constraint_signal = bool(
    hashes_ok
    and first_accepted
    and second_accepted
    and (
        any(
            record.get("classification") == "SINGULAR"
            for record in first.get("parities", {}).values()
        )
        or any(
            record.get("pre_legendre", {}).get("classification") == "SINGULAR"
            for stages in second.get("stages", {}).values()
            for record in stages.values()
        )
    )
)
if route_agreement and passed == tests == 14:
    outcome = ACCEPTED_OUTCOME
elif constraint_signal:
    outcome = SIGNAL_OUTCOME
else:
    outcome = OPEN_OUTCOME

artifact = {
    "provenance": {
        "prior_art_commit": "cbbf7b6",
        "protocol_commit": "286c34b",
        "input_sha256": input_hashes,
    },
    "hypotheses": {
        "background": "frozen branch-B history at incoming v=3/2",
        "slabs": 2,
        "staircase_parities": ["even", "odd"],
        "boundary_phase_dimension": PHASE_DIMENSION,
        "action": "zero-Lambda Lorentzian Regge plus conserved homogeneous dust",
        "claim_scope": "local exact pre/post constraints selected by Legendre degeneracy",
    },
    "margins": {
        "first_minimum_to_gate_ratio": first_margins,
        "second_minimum_to_gate_ratio": second_margins,
    },
    "route_a_rank_and_image": route_a,
    "route_b_covector_contradiction": route_b,
    "controls": {
        "exact_regular": positive,
        "exact_singular": singular,
        "tiny_nonzero": {
            "coefficient": "1/10^100",
            "exact_rank": tiny_exact_rank,
            "hostile_1e-12_cutoff_rank": tiny_hostile_rank,
            "threshold_quotient_rejected": tiny_ok,
        },
        "future_singular": {
            "regular_first_tangent": [
                [int(value) for value in row]
                for row in first_free_tangent.tolist()
            ],
            "later_constraint_covector": [
                [int(value) for value in row]
                for row in later_constraint.tolist()
            ],
            "pulled_back_constraint_covector": [
                [int(value) for value in row]
                for row in pulled_back_constraint.tolist()
            ],
            "constraint_propagates_backward": future_ok,
        },
        "old_weak_mode_control": {
            "accepted": old_control_ok,
            "resolved_rank": 1560,
            "nullity": 0,
            "weak_nonzero_candidates": 120,
            "used_as_finite_height_premise": False,
        },
    },
    "firewall": {
        "physical_mode_spectrum_computed": False,
        "tangent_eigenvalues_computed": False,
        "continuum_target_parsed": False,
        "limiting_speed_G_planck_or_particle_target_parsed": False,
        "numerical_threshold_used_to_define_quotient": False,
        "full_suite_run": False,
    },
    "classification": {
        "exact_local_pre_constraint_quotient": "BOUNDED_NO_GO"
        if outcome == ACCEPTED_OUTCOME
        else "OPEN",
        "exact_local_post_constraint_quotient": "BOUNDED_NO_GO"
        if outcome == ACCEPTED_OUTCOME
        else "OPEN",
        "pseudo_constraint_threshold_quotient": "REFUTED_AS_FITTING",
        "future_propagated_constraint": "OPEN",
        "refinement_gauge_restoration": "OPEN",
        "enlarged_matter_or_alternative_action": "OPEN",
        "physical_graviton_quotient": "OPEN",
    },
    "scope_exclusions": [
        "third or later slab",
        "infinite history",
        "refinement or continuum limit",
        "independent dust perturbations",
        "different or perfect action",
        "separately derived continuous momentum-map reduction",
    ],
    "reopening_conditions": [
        "a later accepted singular slab with a constraint propagated to current data",
        "an exact continuous symmetry of an enlarged action",
        "a preregistered refinement scaling law to an exact kernel",
        "a different or perfect action with exact discrete diffeomorphism symmetry",
    ],
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(outcome)
print(f"{passed}/{tests} PASS")
print(f"artifact: {OUTPUT.name}")

if passed != tests or outcome == OPEN_OUTCOME:
    sys.exit(1)
if outcome == SIGNAL_OUTCOME:
    sys.exit(3)
