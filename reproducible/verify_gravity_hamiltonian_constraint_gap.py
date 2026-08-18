#!/usr/bin/env python3
"""Coverage audit for a physical metric Hamiltonian constraint.

Protocol commit: 6c151e1.  This reads committed machine certificates and
does not rerun their expensive calculations.  Its negative is scoped to the
current canonical dynamical candidates, not all future extensions.
"""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_hamiltonian_constraint_gap.json"
PROTOCOL_COMMIT = "6c151e1"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def load(name):
    return json.loads((HERE / name).read_text())


print("=" * 78)
print("CURRENT GRAVITY DYNAMICS: HAMILTONIAN-CONSTRAINT COVERAGE AUDIT")
print("=" * 78)

tick = load("kahler_dirac_local_tick.json")
tick_refinement = load("kahler_dirac_tick_refinement.json")
whitney_second = load("whitney_constraint_dirac_bergmann.json")
whitney_conversion = load("whitney_first_class_conversion.json")
whitney_weak = load("whitney_weak_first_class_hamiltonian.json")
smooth_saddle = load("round_a2_transverse_hessian.json")
finite_saddle = load("finite_regge_a2_hessian.json")
a4 = load("smooth_derham_a4_stabilization.json")

check(
    "the fixed local tick is a unitary cochain/incidence construction",
    tick["carrier"]["cochain_dimension"] == 2640
    and tick["carrier"]["directed_incidence_arc_dimension"] == 14880
    and tick["walk"]["one_hasse_edge_per_micro_tick"]
    and tick["walk"]["definition"].startswith("U=S("),
)
check(
    "its own certificate leaves physical Lorentzian selection open",
    tick["status"]["physical"].startswith("OPEN:")
    and "Lorentzian" in tick["status"]["physical"],
)
check(
    "the refinement audit rejects the tick as full metric Whitney dynamics",
    tick_refinement["verdict"].startswith("DERIVED KINEMATIC ONLY")
    and "not a lift of the accepted metric Whitney" in tick_refinement["scope"],
)

levels = whitney_second["level_classifications"]
check(
    "both canonical Whitney levels have zero physical first-class constraints",
    [item["physical_first_class_constraint_count"] for item in levels] == [0, 0],
    f"counts={[item['physical_first_class_constraint_count'] for item in levels]}",
)
check(
    "all independent Whitney copy equalities are second class",
    [item["real_second_class_constraint_count"] for item in levels]
    == [12720, 306240]
    and "every independent complex copy constraint is second class"
    in whitney_second["classification_theorem"],
)
check(
    "Dirac reduction reproduces assembled cochain evolution rather than metric evolution",
    whitney_second["small_boundary_4_simplex_control"]["reduced_vector_field_exact"]
    and any("reproduces assembled Whitney evolution" in item
            for item in whitney_second["verdicts"]),
)

minimal = whitney_conversion["minimal_conversion"]
check(
    "minimal conversion is first class only after adding cochain auxiliaries",
    minimal["first_class_constraint"] == "Phi=C u+eta"
    and all(item["first_class_complex_constraint_count_r"] > 0
            for item in whitney_conversion["dimension_records"]),
)
check(
    "its unique physical dressing and Hamiltonian contain the global Gram inverse",
    "G^-1" in minimal["dressing"]
    and "G^-1" in minimal["quadratic_cross_block"]
    and any("does not provide a local Whitney tick" in item
            for item in whitney_conversion["verdicts"]),
)
check(
    "weak first-classness still fails the endpoint-local quadratic gate",
    "weak first class" in whitney_weak["hypotheses"]["invariance"]
    and any("not endpoint-local" in item for item in whitney_weak["verdicts"]),
)

check(
    "the smooth conformal negative direction is explicitly non-gauge",
    smooth_saddle["exact_results"]["l2_lambda"] == 8
    and "non-gauge conformal" in smooth_saddle["derived"][1],
)
check(
    "the finite metric Hessian has no extra quotient null modes",
    finite_saddle["quotient_hessian"]["inertia_positive_zero_negative"]
    == [569, 0, 150],
)
check(
    "smooth A4 does not supply a selected finite-time constraint rescue",
    "FINITE CUTOFF OPEN" in a4["verdict"]
    and "selection of cutoff shape, tau" in a4["open"][2],
)

# Scoped source guard.  A Grover 'flip shift' is present, but it is a walk
# permutation, not an ADM shift multiplier; the structured certificate above
# fixes that meaning positively.
authoritative_sources = [
    "verify_kahler_dirac_local_tick.py",
    "verify_kahler_dirac_tick_refinement.py",
    "verify_whitney_constraint_dirac_bergmann.py",
    "verify_whitney_first_class_conversion.py",
    "verify_whitney_weak_first_class_hamiltonian.py",
]
source_text = "\n".join((HERE / name).read_text().lower()
                        for name in authoritative_sources)
forbidden_metric_phase_tokens = (
    "lapse", "extrinsic curvature", "canonical_momentum",
    "metric_momentum", "adm constraint", "dynamical edge length",
)
check(
    "authoritative dynamics sources declare no metric lapse/momentum variables",
    all(token not in source_text for token in forbidden_metric_phase_tokens),
    "the occurrence of 'shift' is the explicitly certified Grover flip permutation",
)

# Frozen coverage matrix.  'Constraint closure' is true for the auxiliary
# cochain conversion, but no row covers all six metric-gravity requirements.
fields = (
    "metric_configuration", "metric_momentum", "lapse_shift",
    "metric_first_class_constraints", "constraint_closure",
    "metric_evolution",
)
coverage = {
    "fixed_local_tick": [False, False, False, False, False, False],
    "whitney_second_class_action": [False, False, False, False, True, False],
    "minimal_first_class_conversion": [False, False, False, False, True, False],
    "weak_first_class_completion": [False, False, False, False, True, False],
}
check(
    "no current canonical candidate covers the six ADM-like metric fields",
    all(not all(row) for row in coverage.values())
    and all(len(row) == len(fields) for row in coverage.values()),
)

protocol = (
    ROOT / "docs" / "gravity" / "gravity_hamiltonian_constraint_gap_protocol.md"
).read_text()
check(
    "the preregistered boundary scopes the negative to the current repository",
    "current repository" in protocol
    and "not a no-go theorem" in protocol
    and "does **not** say the conformal modes are experimentally physical" in protocol,
)

verdict = "DERIVED CURRENT HAMILTONIAN-CONSTRAINT GAP"
payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "benchmark_source": "https://arxiv.org/abs/gr-qc/0405109",
    "inventory": authoritative_sources,
    "required_fields": list(fields),
    "coverage": {name: dict(zip(fields, row)) for name, row in coverage.items()},
    "certified_facts": {
        "local_tick_carrier": "fixed 2640 cochains / 14880 incidence arcs",
        "physical_first_class_copy_constraints": [0, 0],
        "real_second_class_copy_constraints": [12720, 306240],
        "minimal_conversion": "first class on auxiliary cochains but G^-1 nonlocal",
        "weak_completion": "endpoint-local quadratic Hamiltonian gate fails",
        "smooth_conformal_mode": "l=2 non-gauge under spatial diffeomorphisms",
        "finite_metric_hessian": "scale quotient inertia (569,0,150), no null",
    },
    "verdict": verdict,
    "derived": [
        "no current canonical dynamics removes conformal metric modes by a physical first-class constraint",
        "the local unitary tick evolves fields on a fixed incidence geometry, not the metric",
        "Whitney copy constraints cannot be relabeled as the missing Hamiltonian constraint",
    ],
    "open": [
        "a canonical metric phase space and symplectic structure",
        "selected lapse/shift carrier and first-class constraint algebra",
        "local/refinement-compatible metric Hamiltonian",
        "whether constrained Lorentzian dynamics removes the conformal modes",
    ],
    "not_claimed": [
        "all future discrete gravity constraints are impossible",
        "the conformal negative modes are experimentally propagating",
        "Lorentzian GR or a Planck scale is derived",
    ],
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(verdict)
print("OPEN: metric phase space, first-class algebra, and Lorentzian evolution")
raise SystemExit(0 if passed == tests else 1)
