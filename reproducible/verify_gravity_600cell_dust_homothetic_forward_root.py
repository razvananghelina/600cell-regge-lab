#!/usr/bin/env python3
"""Bracketed non-static homothetic forward root for the 600-cell dust slab.

Prior-art commit: 2ce58bd.
Protocol commit: 29653b9.

The bracket and 80-step bisection are frozen.  Acceptance requires all 35
internal equations and the independently committed 30-component canonical
junction; a zero of the summed lapse equation alone does not pass.
"""

import ast
import contextlib
import io
import json
from pathlib import Path

import mpmath as arb


HERE = Path(__file__).resolve().parent
RESPONSE_SOURCE = (
    HERE / "verify_gravity_600cell_dust_homothetic_mass_conservation.py"
)
RESPONSE_ARTIFACT = (
    HERE / "gravity_600cell_dust_homothetic_mass_conservation.json"
)
GLUING_ARTIFACT = HERE / "gravity_600cell_dust_two_slab_gluing.json"
OUTPUT = HERE / "gravity_600cell_dust_homothetic_forward_root.json"
PRIOR_ART_COMMIT = "2ce58bd"
PROTOCOL_COMMIT = "29653b9"
RESPONSE_RESULT_COMMIT = "086009a"
DPS = 100
arb.mp.dps = DPS
S_LEFT = -arb.mpf(1)/40000
S_RIGHT = -arb.mpf(1)/640000
BISECTIONS = 80
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


def load_response_prefix():
    """Load the committed evaluator through its fresh static controls only."""
    tree = ast.parse(RESPONSE_SOURCE.read_text(), filename=str(RESPONSE_SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "scale_points"
            for target in node.targets
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("response evaluator cutoff was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(RESPONSE_SOURCE),
        "__name__": "forward_root_imported_response_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(RESPONSE_SOURCE), "exec"), namespace)
    return namespace


def norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def infinity_norm(values):
    return max(abs(value) for value in values)


def spread(values):
    mean = sum(values, arb.mpf(0))/len(values)
    return max(abs(value-mean) for value in values)


def text(value, digits=50):
    return arb.nstr(value, digits)


def real_vector(raw):
    return tuple(arb.mpf(item["real"]) for item in raw)


print("Homothetic 600-cell forward-root audit", flush=True)
response_core = load_response_prefix()
evaluate_homothetic = response_core["evaluate_homothetic"]
models = response_core["models"]
core = response_core["core"]

check(
    "the imported evaluator retains its six exact/static controls",
    response_core["tests"] == response_core["passed"] == 6,
)
check(
    "the root mission provenance and two parity carriers are fixed",
    PRIOR_ART_COMMIT == "2ce58bd"
    and PROTOCOL_COMMIT == "29653b9"
    and RESPONSE_RESULT_COMMIT == "086009a"
    and set(models) == {"even", "odd"},
)

response_artifact = json.loads(RESPONSE_ARTIFACT.read_text())
gluing_artifact = json.loads(GLUING_ARTIFACT.read_text())
upstream_ok = bool(
    response_artifact["outcome"]
        == "LOCAL_HOMOTHETIC_STATIC_ONLY_GLOBAL_AND_LOCAL"
    and response_artifact["passed"] == response_artifact["tests"] == 10
    and response_artifact["exact_geometry"]["passed"] is True
    and gluing_artifact["outcome"] == "TWO_SLAB_GLUING_CONTROL_PASSED"
    and gluing_artifact["passed"] == gluing_artifact["tests"]
)
check("the response and two-slab artifacts authorize the root test", upstream_ok)


targets = {}
target_uncertainties = {}
target_maps_ok = True
for parity in models:
    parity_data = gluing_artifact["parities"][parity]
    old_to_final = tuple(parity_data["geometry"]["old_to_final_orbit_map"])
    target_maps_ok &= sorted(old_to_final) == list(range(30))
    static_post = real_vector(parity_data["momenta"]["post"])
    targets[parity] = tuple(static_post[index] for index in old_to_final)
    target_uncertainties[parity] = arb.mpf(
        parity_data["momenta"]["cusp_uncertainty_norm"]
    )
check(
    "the committed old-to-final maps are permutations and define 30 targets",
    target_maps_ok and all(len(target) == 30 for target in targets.values()),
)


def sign(value):
    return -1 if value < 0 else 1 if value > 0 else 0


root_records = {}
all_branch_ok = True
brackets_ok = True
bisections_ok = True

for parity, model in models.items():
    print(f"  {parity}: evaluating frozen bracket", flush=True)
    left = S_LEFT
    right = S_RIGHT
    left_record = evaluate_homothetic(model, left)
    right_record = evaluate_homothetic(model, right)
    all_branch_ok &= left_record["branch_pass"] and right_record["branch_pass"]
    bracket_ok = bool(
        left < right < 0
        and sign(left_record["lapse"])*sign(right_record["lapse"]) == -1
    )
    brackets_ok &= bracket_ok
    print(
        "    E(left)={} E(right)={} bracket={}".format(
            text(left_record["lapse"], 12),
            text(right_record["lapse"], 12),
            "PASS" if bracket_ok else "FAIL",
        ),
        flush=True,
    )

    evaluated = 2
    early_exact = False
    if bracket_ok:
        left_sign = sign(left_record["lapse"])
        for iteration in range(BISECTIONS):
            midpoint = (left+right)/2
            midpoint_record = evaluate_homothetic(model, midpoint)
            evaluated += 1
            all_branch_ok &= midpoint_record["branch_pass"]
            midpoint_sign = sign(midpoint_record["lapse"])
            if midpoint_sign == 0:
                left = right = midpoint
                left_record = right_record = midpoint_record
                early_exact = True
                break
            if midpoint_sign == left_sign:
                left = midpoint
                left_record = midpoint_record
                left_sign = midpoint_sign
            else:
                right = midpoint
                right_record = midpoint_record
            if (iteration+1) % 20 == 0:
                print(
                    f"    bisection {iteration+1:2d}: width={float(right-left):.3e}",
                    flush=True,
                )

        root = (left+right)/2
        root_record = evaluate_homothetic(model, root)
        evaluated += 1
        all_branch_ok &= root_record["branch_pass"]
        width = right-left
        root_gate = bool(
            width < arb.mpf("2e-28")
            and abs(root_record["lapse"]) < arb.mpf("1e-25")
            and abs(root) > arb.mpf("1e-7")
            and root < 0
        )
        bisections_ok &= root_gate
    else:
        root = None
        root_record = None
        width = None
        root_gate = False
        bisections_ok = False

    if root_record is not None:
        diagonal = root_record["local"][:30]
        poles = root_record["local"][30:]
        internal_gate = bool(
            infinity_norm(diagonal) < arb.mpf("1e-60")
            and infinity_norm(poles) < arb.mpf("1e-25")
            and infinity_norm(root_record["local"]) < arb.mpf("1e-25")
            and spread(diagonal) < arb.mpf("1e-60")
            and spread(poles) < arb.mpf("1e-60")
        )
        momentum_residual = tuple(
            value-target
            for value, target in zip(root_record["pre"], targets[parity])
        )
        momentum_norm = norm(momentum_residual)
        momentum_spread = spread(momentum_residual)
        momentum_bound = 10*target_uncertainties[parity]
        canonical_gate = bool(
            momentum_norm <= momentum_bound
            and momentum_spread <= momentum_bound
        )
        expected = core["ARB_EPSILON_3"]*core["ARB_L0"]*core["ARB_TAU"]/4
        exact_formula_error = infinity_norm(tuple(
            value-expected for value in root_record["pre"]
        ))
        print(
            "    root={} width={} max35={} ||junction||={}".format(
                text(root, 18), text(width, 6),
                text(infinity_norm(root_record["local"]), 6),
                text(momentum_norm, 6),
            ),
            flush=True,
        )
    else:
        diagonal = poles = momentum_residual = ()
        internal_gate = canonical_gate = False
        momentum_norm = momentum_spread = momentum_bound = None
        exact_formula_error = None

    root_records[parity] = {
        "left": left,
        "right": right,
        "width": width,
        "root": root,
        "record": root_record,
        "bracket_pass": bracket_ok,
        "root_gate": root_gate,
        "internal_gate": internal_gate,
        "canonical_gate": canonical_gate,
        "momentum_residual": momentum_residual,
        "momentum_norm": momentum_norm,
        "momentum_spread": momentum_spread,
        "momentum_bound": momentum_bound,
        "exact_formula_error": exact_formula_error,
        "evaluations": evaluated,
        "early_exact": early_exact,
    }

check(
    "every evaluated root state remains on the certified Lorentzian branch",
    all_branch_ok,
)


if all(root_records[parity]["root"] is not None for parity in models):
    even = root_records["even"]
    odd = root_records["odd"]
    overlap = bool(
        even["left"]-arb.mpf("1e-28")
            <= odd["right"]+arb.mpf("1e-28")
        and odd["left"]-arb.mpf("1e-28")
            <= even["right"]+arb.mpf("1e-28")
    )
    root_difference = abs(even["root"]-odd["root"])
    scale_difference = abs(arb.exp(even["root"])-arb.exp(odd["root"]))
    pre_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["record"]["pre"], odd["record"]["pre"]
        )
    ))
    post_difference = infinity_norm(tuple(
        left-right for left, right in zip(
            even["record"]["post"], odd["record"]["post"]
        )
    ))
    parity_gate = bool(
        overlap
        and root_difference < arb.mpf("1e-27")
        and scale_difference < arb.mpf("1e-27")
        and pre_difference < arb.mpf("1e-24")
        and post_difference < arb.mpf("1e-24")
    )
else:
    overlap = parity_gate = False
    root_difference = scale_difference = pre_difference = post_difference = None


controls_ok = bool(upstream_ok and target_maps_ok and all_branch_ok)
internal_ok = all(record["internal_gate"] for record in root_records.values())
canonical_ok = all(record["canonical_gate"] for record in root_records.values())

if not controls_ok:
    outcome = "HOMOTHETIC_FORWARD_CONTROL_FAILED"
elif not brackets_ok:
    outcome = "HOMOTHETIC_FORWARD_BRACKET_REFUTED"
elif not bisections_ok:
    outcome = "HOMOTHETIC_FORWARD_ROOT_NUMERICALLY_OPEN"
elif not internal_ok:
    outcome = "HOMOTHETIC_MINISUPERSPACE_ONLY"
elif not canonical_ok:
    outcome = "HOMOTHETIC_STATIONARY_NOT_CANONICAL"
elif not parity_gate:
    outcome = "HOMOTHETIC_FORWARD_SCHEDULE_DEPENDENT"
else:
    outcome = "HOMOTHETIC_FORWARD_ROOT_ACCEPTED"

check(
    "the frozen hierarchy assigns exactly one forward-root outcome",
    outcome in {
        "HOMOTHETIC_FORWARD_CONTROL_FAILED",
        "HOMOTHETIC_FORWARD_BRACKET_REFUTED",
        "HOMOTHETIC_FORWARD_ROOT_NUMERICALLY_OPEN",
        "HOMOTHETIC_MINISUPERSPACE_ONLY",
        "HOMOTHETIC_STATIONARY_NOT_CANONICAL",
        "HOMOTHETIC_FORWARD_SCHEDULE_DEPENDENT",
        "HOMOTHETIC_FORWARD_ROOT_ACCEPTED",
    },
    outcome,
)


s_cont = -core["ARB_ZETA"]**2*core["ARB_RHO"]/(2*core["ARB_L0_SQUARE"])


def serialize_branch(record):
    return {
        "negative_direction_counts": {
            "NONE" if key is None else str(key): value
            for key, value in record["branch"]["negative_counts"].items()
        },
        "minimum_leading_minor": text(
            record["branch"]["minimum_leading_minor"], 40
        ),
        "minimum_angle_argument": text(
            record["branch"]["minimum_argument"], 40
        ),
        "maximum_imaginary": text(record["maximum_imaginary"], 40),
        "passed": record["branch_pass"],
    }


serialized_roots = {}
for parity, result in root_records.items():
    record = result["record"]
    serialized_roots[parity] = {
        "initial_bracket": [text(S_LEFT, 30), text(S_RIGHT, 30)],
        "final_bracket": (
            None if result["root"] is None
            else [text(result["left"], 60), text(result["right"], 60)]
        ),
        "width": None if result["width"] is None else text(result["width"], 40),
        "root": None if result["root"] is None else text(result["root"], 60),
        "scale_ratio": (
            None if result["root"] is None else text(arb.exp(result["root"]), 60)
        ),
        "root_over_continuum_control": (
            None if result["root"] is None else text(result["root"]/s_cont, 40)
        ),
        "evaluations": result["evaluations"],
        "early_exact": result["early_exact"],
        "bracket_pass": result["bracket_pass"],
        "root_gate": result["root_gate"],
        "internal_gate": result["internal_gate"],
        "canonical_gate": result["canonical_gate"],
        "lapse_equation": (
            None if record is None else text(record["lapse"], 50)
        ),
        "diagonal_residuals": (
            None if record is None
            else [text(value, 40) for value in record["local"][:30]]
        ),
        "pole_residuals": (
            None if record is None
            else [text(value, 40) for value in record["local"][30:]]
        ),
        "maximum_internal_residual": (
            None if record is None
            else text(infinity_norm(record["local"]), 40)
        ),
        "pre_momentum": (
            None if record is None
            else [text(value, 50) for value in record["pre"]]
        ),
        "post_momentum": (
            None if record is None
            else [text(value, 50) for value in record["post"]]
        ),
        "target_momentum": [text(value, 50) for value in targets[parity]],
        "junction_residual": [
            text(value, 40) for value in result["momentum_residual"]
        ],
        "junction_norm": (
            None if result["momentum_norm"] is None
            else text(result["momentum_norm"], 40)
        ),
        "junction_spread": (
            None if result["momentum_spread"] is None
            else text(result["momentum_spread"], 40)
        ),
        "junction_bound": (
            None if result["momentum_bound"] is None
            else text(result["momentum_bound"], 40)
        ),
        "exact_static_formula_error": (
            None if result["exact_formula_error"] is None
            else text(result["exact_formula_error"], 40)
        ),
        "branch": None if record is None else serialize_branch(record),
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "response_result_commit": RESPONSE_RESULT_COMMIT,
    "frozen_parameters": {
        "dps": DPS,
        "bisections": BISECTIONS,
        "tau": text(core["ARB_TAU"], 30),
        "rho": text(core["ARB_RHO"], 40),
        "mass": text(core["ARB_MASS"], 50),
        "continuum_scale_log_control": text(s_cont, 50),
    },
    "roots": serialized_roots,
    "parity_gate": {
        "overlap": overlap,
        "root_difference": None if root_difference is None else text(root_difference, 40),
        "scale_ratio_difference": (
            None if scale_difference is None else text(scale_difference, 40)
        ),
        "pre_momentum_infinity_difference": (
            None if pre_difference is None else text(pre_difference, 40)
        ),
        "post_momentum_infinity_difference": (
            None if post_difference is None else text(post_difference, 40)
        ),
        "passed": parity_gate,
    },
    "classification": {
        "all_35_internal_equations_pass": internal_ok,
        "canonical_junction_pass": canonical_ok,
        "parity_pass": parity_gate,
        "lapse_selected": False,
        "chosen_lapse": True,
        "multi_tick_evolution": False,
        "physical_clock": "OPEN",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)
if passed != tests:
    raise SystemExit(1)
