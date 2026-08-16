#!/usr/bin/env python3
"""Local conserved-mass homothetic audit of the 600-cell dust slab.

Prior-art commit: 8865346.
Protocol commits: 428330e, ff8d352.

No root, mass, internal length or derivative coefficient is fitted.  The
only non-static evaluations are the six scale displacements frozen in the
protocol and their lapse-chain validation points.
"""

import ast
from collections import Counter
import contextlib
import io
import json
from pathlib import Path

import mpmath as arb
import sympy as sp


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_gravity_600cell_dust_canonical_continuation.py"
UPSTREAM = HERE / "gravity_600cell_dust_regular_lapse_identity.json"
OUTPUT = HERE / "gravity_600cell_dust_homothetic_mass_conservation.json"
PRIOR_ART_COMMIT = "8865346"
PROTOCOL_COMMIT = "428330e"
PROTOCOL_CLARIFICATION_COMMIT = "ff8d352"
DPS = 100
arb.mp.dps = DPS
H = arb.mpf("1e-4")
K = arb.mpf("1e-6")
ARITHMETIC_FLOOR = arb.mpf("1e-70")
NONZERO_FACTOR = arb.mpf(100)
SCHEDULE_FACTOR = arb.mpf(10)
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


def load_continuation_core():
    """Execute only the definition/control prefix of the committed verifier."""
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    cut = None
    for index, node in enumerate(tree.body):
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "print"
        ):
            cut = index
            break
    if cut is None:
        raise RuntimeError("continuation main marker was not found")
    prefix = ast.Module(body=tree.body[:cut], type_ignores=[])
    namespace = {
        "__file__": str(SOURCE),
        "__name__": "homothetic_imported_continuation_core",
    }
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(prefix, str(SOURCE), "exec"), namespace)
    return namespace


def vector_norm(values):
    return arb.sqrt(sum(abs(value)**2 for value in values))


def maximum_imaginary(action, gradient):
    return max(abs(arb.im(action)), *(abs(arb.im(value)) for value in gradient))


def relative_error(left, right):
    return abs(left-right)/max(arb.mpf(1), abs(left), abs(right))


def string(value, digits=50):
    return arb.nstr(value, digits)


print("Conserved-mass homothetic 600-cell slab audit", flush=True)
core = load_continuation_core()
models = core["models"]
action_and_gradient = core["action_and_gradient"]
branch_pass = core["branch_pass"]

check(
    "the imported canonical-action prefix retains its four controls",
    core["tests"] == core["passed"] == 4,
)
check(
    "the prior-art, protocol and parity carrier provenance is frozen",
    PRIOR_ART_COMMIT == "8865346"
    and PROTOCOL_COMMIT == "428330e"
    and PROTOCOL_CLARIFICATION_COMMIT == "ff8d352"
    and set(models) == {"even", "odd"},
)

upstream = json.loads(UPSTREAM.read_text())
upstream_ok = bool(
    upstream["outcome"] == "REGULAR_LAPSE_IDENTITY_PROVED"
    and upstream["passed"] == upstream["tests"] == 13
    and upstream["protocol_commit"] == "cf492b9"
)
check("the exact regular-lapse theorem is an upstream control", upstream_ok)


# Exact geometry, independent of the action implementation.
phi = (sp.Integer(1)+sp.sqrt(5))/2
L_minus, L_plus, rho_symbol = sp.symbols(
    "L_minus L_plus rho", positive=True
)
adjacent_dot = phi/2
radius_minus = phi*L_minus
radius_plus = phi*L_plus
edge_square = sp.factor(2*radius_minus**2*(1-adjacent_dot))
outer_time_square = sp.factor(
    rho_symbol+(radius_plus-radius_minus)**2
)
cross_square = sp.factor(
    radius_minus**2+radius_plus**2
    - 2*radius_minus*radius_plus*adjacent_dot
    - outer_time_square
)
geometry_ok = bool(
    sp.simplify(2*(1-adjacent_dot)-phi**-2) == 0
    and sp.simplify(edge_square-L_minus**2) == 0
    and sp.simplify(cross_square-(L_minus*L_plus-rho_symbol)) == 0
)
check(
    "the homothetic diagonal L_minus*L_plus-rho follows exactly",
    geometry_ok,
    f"edge^2={edge_square}; cross^2={cross_square}",
)


ARB_L0 = core["ARB_L0"]
ARB_L0_SQUARE = core["ARB_L0_SQUARE"]
ARB_RHO = core["ARB_RHO"]
ARB_EPSILON_3 = core["ARB_EPSILON_3"]
ARB_TAU = core["ARB_TAU"]
ARB_BASE_OLD = core["ARB_BASE_OLD"]


def evaluate_homothetic(model, scale_log, rho_log=None):
    """Evaluate the complete action on the fixed geometric homothetic state."""
    if rho_log is None:
        rho_value = ARB_RHO
    else:
        rho_value = arb.exp(rho_log)
    scale = arb.exp(scale_log)
    diagonal = scale*ARB_L0_SQUARE-rho_value
    final_square = scale**2*ARB_L0_SQUARE
    if diagonal <= 0 or rho_value <= 0 or final_square <= 0:
        raise ValueError("homothetic magnitude left the positive domain")
    internal = tuple([diagonal]*30+[rho_value]*5)
    final = tuple([final_square]*30)
    action, gradient, branch = action_and_gradient(
        model, ARB_BASE_OLD, internal, final
    )
    local = tuple(arb.re(value) for value in gradient[30:65])
    pre = tuple(-arb.re(value) for value in gradient[:30])
    post = tuple(arb.re(value) for value in gradient[65:95])
    lapse = (
        sum(local[:30], arb.mpf(0))*(-rho_value/diagonal)
        + sum(local[30:], arb.mpf(0))
    )
    imaginary = maximum_imaginary(action, gradient)
    return {
        "scale_log": scale_log,
        "rho": rho_value,
        "diagonal": diagonal,
        "final_square": final_square,
        "action": action,
        "gradient": gradient,
        "local": local,
        "pre": pre,
        "post": post,
        "lapse": lapse,
        "branch": branch,
        "maximum_imaginary": imaginary,
        "branch_pass": branch_pass(branch, imaginary),
    }


def static_control(parity, model):
    record = evaluate_homothetic(model, arb.mpf(0))
    expected_pre = -ARB_EPSILON_3*ARB_L0*ARB_TAU/4
    expected_post = -expected_pre
    local_error = max(abs(value) for value in record["local"])
    pre_error = max(relative_error(value, expected_pre) for value in record["pre"])
    post_error = max(
        relative_error(value, expected_post) for value in record["post"]
    )
    ok = bool(
        record["branch_pass"]
        and local_error < arb.mpf("1e-60")
        and pre_error < arb.mpf("1e-60")
        and post_error < arb.mpf("1e-60")
    )
    check(
        f"{parity}: a fresh static evaluation reproduces all 35 equations and momenta",
        ok,
        "max local={} pre rel={} post rel={}".format(
            string(local_error, 8), string(pre_error, 8), string(post_error, 8)
        ),
    )
    return record, ok, local_error, pre_error, post_error


static_records = {}
static_ok = True
for parity, model in models.items():
    record, ok, *_ = static_control(parity, model)
    static_records[parity] = record
    static_ok &= ok


scale_points = (-H, -H/2, -H/4, H/4, H/2, H)
records = {}
chain_records = {}
all_branch_ok = True
all_chain_ok = True
base_rho_log = arb.log(ARB_RHO)

for parity, model in models.items():
    print(f"  evaluating non-static {parity} schedule", flush=True)
    records[parity] = {}
    chain_records[parity] = {}
    for scale_log in scale_points:
        key = string(scale_log, 20)
        central = evaluate_homothetic(model, scale_log)
        records[parity][key] = central
        all_branch_ok &= central["branch_pass"]

        actions = {}
        validation_branch = True
        for lapse_step in (K, K/2):
            for sign in (-1, 1):
                displaced = evaluate_homothetic(
                    model, scale_log, base_rho_log+sign*lapse_step
                )
                actions[(string(lapse_step, 20), sign)] = arb.re(
                    displaced["action"]
                )
                validation_branch &= displaced["branch_pass"]
                all_branch_ok &= displaced["branch_pass"]

        q_k = (
            actions[(string(K, 20), 1)]-actions[(string(K, 20), -1)]
        )/(2*K)
        q_half = (
            actions[(string(K/2, 20), 1)]
            - actions[(string(K/2, 20), -1)]
        )/K
        chain_difference = abs(q_half-24*central["lapse"])
        chain_proxy = 10*abs(q_half-q_k)+arb.mpf("1e-60")
        chain_ok = bool(validation_branch and chain_difference <= chain_proxy)
        all_chain_ok &= chain_ok
        chain_records[parity][key] = {
            "q_k": q_k,
            "q_half": q_half,
            "difference": chain_difference,
            "bound": chain_proxy,
            "passed": chain_ok,
        }
        print(
            f"    s={float(scale_log):+.2e} "
            f"local={float(vector_norm(central['local'])):.3e} "
            f"E_lapse={float(central['lapse']):+.3e} "
            f"chain={'PASS' if chain_ok else 'FAIL'}",
            flush=True,
        )

check(
    "all frozen displaced evaluations retain the Lorentzian branch",
    all_branch_ok,
)
check(
    "the independent restricted-action derivatives reproduce 24*E_lapse",
    all_chain_ok,
)


def central_derivative(positive, negative, step):
    if isinstance(positive, tuple):
        return tuple(
            (left-right)/(2*step) for left, right in zip(positive, negative)
        )
    return (positive-negative)/(2*step)


def richardson_record(parity, field):
    parity_records = records[parity]

    def value(scale):
        return parity_records[string(scale, 20)][field]

    d_h = central_derivative(value(H), value(-H), H)
    d_half = central_derivative(value(H/2), value(-H/2), H/2)
    d_quarter = central_derivative(value(H/4), value(-H/4), H/4)

    if isinstance(d_h, tuple):
        r4_h = tuple((4*middle-large)/3 for middle, large in zip(d_half, d_h))
        operational = tuple(
            (4*small-middle)/3 for small, middle in zip(d_quarter, d_half)
        )
        proxy = tuple(
            abs(small-large)+ARITHMETIC_FLOOR
            for small, large in zip(operational, r4_h)
        )
        norm = vector_norm(operational)
        proxy_norm = vector_norm(proxy)
        component_resolved = tuple(
            abs(value) > NONZERO_FACTOR*error
            for value, error in zip(operational, proxy)
        )
        resolved = bool(
            norm > NONZERO_FACTOR*proxy_norm and any(component_resolved)
        )
        return {
            "operational": operational,
            "proxy": proxy,
            "norm": norm,
            "proxy_norm": proxy_norm,
            "component_resolved": component_resolved,
            "resolved_nonzero": resolved,
        }

    r4_h = (4*d_half-d_h)/3
    operational = (4*d_quarter-d_half)/3
    proxy = abs(operational-r4_h)+ARITHMETIC_FLOOR
    resolved = bool(abs(operational) > NONZERO_FACTOR*proxy)
    return {
        "operational": operational,
        "proxy": proxy,
        "resolved_nonzero": resolved,
    }


derivatives = {
    parity: {
        field: richardson_record(parity, field)
        for field in ("local", "lapse", "pre", "post")
    }
    for parity in models
}


def vector_schedule_gate(left, right):
    difference = vector_norm(tuple(
        a-b for a, b in zip(left["operational"], right["operational"])
    ))
    bound = (
        SCHEDULE_FACTOR*(left["proxy_norm"]+right["proxy_norm"])
        + ARITHMETIC_FLOOR
    )
    return bool(difference <= bound), difference, bound


def scalar_schedule_gate(left, right):
    difference = abs(left["operational"]-right["operational"])
    bound = SCHEDULE_FACTOR*(left["proxy"]+right["proxy"])+ARITHMETIC_FLOOR
    return bool(difference <= bound), difference, bound


local_schedule_ok, local_schedule_difference, local_schedule_bound = (
    vector_schedule_gate(derivatives["even"]["local"], derivatives["odd"]["local"])
)
lapse_schedule_ok, lapse_schedule_difference, lapse_schedule_bound = (
    scalar_schedule_gate(derivatives["even"]["lapse"], derivatives["odd"]["lapse"])
)
schedule_ok = bool(local_schedule_ok and lapse_schedule_ok)
check(
    "the even/odd local and global decision derivatives agree within proxies",
    schedule_ok,
    "local diff/bound={}/{}; lapse diff/bound={}/{}".format(
        string(local_schedule_difference, 8), string(local_schedule_bound, 8),
        string(lapse_schedule_difference, 8), string(lapse_schedule_bound, 8),
    ),
)


local_resolved = all(
    derivatives[parity]["local"]["resolved_nonzero"] for parity in models
)
lapse_resolved = all(
    derivatives[parity]["lapse"]["resolved_nonzero"] for parity in models
)
controls_ok = bool(
    upstream_ok and geometry_ok and static_ok and all_branch_ok and all_chain_ok
)

if not controls_ok:
    outcome = "HOMOTHETIC_CONTROL_FAILED"
elif not schedule_ok:
    outcome = "HOMOTHETIC_SCHEDULE_DEPENDENT"
elif local_resolved and lapse_resolved:
    outcome = "LOCAL_HOMOTHETIC_STATIC_ONLY_GLOBAL_AND_LOCAL"
elif local_resolved and not lapse_resolved:
    outcome = "LOCAL_HOMOTHETIC_STATIC_ONLY_LOCAL"
elif lapse_resolved and not local_resolved:
    outcome = "LOCAL_HOMOTHETIC_GLOBAL_ONLY_LOCAL_TANGENT"
elif not local_resolved and not lapse_resolved:
    outcome = "LOCAL_HOMOTHETIC_TANGENT_SURVIVES"
else:
    outcome = "HOMOTHETIC_NUMERICALLY_OPEN"

check(
    "the frozen hierarchy assigns one mechanical outcome",
    outcome in {
        "HOMOTHETIC_CONTROL_FAILED",
        "HOMOTHETIC_SCHEDULE_DEPENDENT",
        "LOCAL_HOMOTHETIC_STATIC_ONLY_GLOBAL_AND_LOCAL",
        "LOCAL_HOMOTHETIC_STATIC_ONLY_LOCAL",
        "LOCAL_HOMOTHETIC_GLOBAL_ONLY_LOCAL_TANGENT",
        "LOCAL_HOMOTHETIC_TANGENT_SURVIVES",
        "HOMOTHETIC_NUMERICALLY_OPEN",
    },
    outcome,
)


def serialize_branch(record):
    return {
        "negative_direction_counts": {
            "NONE" if key is None else str(key): value
            for key, value in record["branch"]["negative_counts"].items()
        },
        "minimum_leading_minor": string(
            record["branch"]["minimum_leading_minor"], 40
        ),
        "minimum_angle_argument": string(
            record["branch"]["minimum_argument"], 40
        ),
        "maximum_imaginary": string(record["maximum_imaginary"], 40),
        "passed": record["branch_pass"],
    }


def serialize_derivative(record):
    if isinstance(record["operational"], tuple):
        return {
            "operational": [string(value, 50) for value in record["operational"]],
            "proxy": [string(value, 30) for value in record["proxy"]],
            "norm": string(record["norm"], 50),
            "proxy_norm": string(record["proxy_norm"], 30),
            "component_resolved": list(record["component_resolved"]),
            "resolved_nonzero": record["resolved_nonzero"],
        }
    return {
        "operational": string(record["operational"], 50),
        "proxy": string(record["proxy"], 30),
        "resolved_nonzero": record["resolved_nonzero"],
    }


artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_clarification_commit": PROTOCOL_CLARIFICATION_COMMIT,
    "upstream_static_outcome": upstream["outcome"],
    "exact_geometry": {
        "adjacent_dot": str(adjacent_dot),
        "unit_circumradius_edge_square": str(sp.factor(edge_square/L_minus**2)),
        "circumradius": "phi*L",
        "outer_time_square": str(outer_time_square),
        "cross_diagonal_square": str(cross_square),
        "passed": geometry_ok,
    },
    "frozen_parameters": {
        "dps": DPS,
        "scale_log_h": string(H, 20),
        "lapse_log_k": string(K, 20),
        "L0": string(ARB_L0, 50),
        "rho0": string(ARB_RHO, 50),
        "mass": string(core["ARB_MASS"], 50),
    },
    "static_controls": {
        parity: {
            "local_norm": string(vector_norm(record["local"]), 40),
            "pre_mean": string(sum(record["pre"])/30, 40),
            "post_mean": string(sum(record["post"])/30, 40),
            "branch": serialize_branch(record),
        }
        for parity, record in static_records.items()
    },
    "displaced_states": {
        parity: {
            key: {
                "rho": string(record["rho"], 40),
                "diagonal": string(record["diagonal"], 40),
                "final_square": string(record["final_square"], 40),
                "local": [string(value, 40) for value in record["local"]],
                "local_norm": string(vector_norm(record["local"]), 40),
                "lapse_equation": string(record["lapse"], 40),
                "pre": [string(value, 40) for value in record["pre"]],
                "post": [string(value, 40) for value in record["post"]],
                "branch": serialize_branch(record),
                "chain_validation": {
                    item: (
                        value if isinstance(value, bool) else string(value, 40)
                    )
                    for item, value in chain_records[parity][key].items()
                },
            }
            for key, record in records[parity].items()
        }
        for parity in models
    },
    "derivatives": {
        parity: {
            field: serialize_derivative(record)
            for field, record in parity_records.items()
        }
        for parity, parity_records in derivatives.items()
    },
    "schedule_gate": {
        "local_passed": local_schedule_ok,
        "local_difference": string(local_schedule_difference, 40),
        "local_bound": string(local_schedule_bound, 40),
        "lapse_passed": lapse_schedule_ok,
        "lapse_difference": string(lapse_schedule_difference, 40),
        "lapse_bound": string(lapse_schedule_bound, 40),
        "passed": schedule_ok,
    },
    "classification": {
        "local_derivative_resolved_nonzero": local_resolved,
        "global_lapse_derivative_resolved_nonzero": lapse_resolved,
        "physical_tick_accepted": False,
        "scope": "local rigid homothetic ansatz at the published sandwich",
    },
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("\nDerivative summary", flush=True)
for parity in models:
    local = derivatives[parity]["local"]
    lapse = derivatives[parity]["lapse"]
    print(
        f"  {parity}: ||dR/ds||={float(local['norm']):.12e} "
        f"proxy={float(local['proxy_norm']):.3e} "
        f"dE/ds={float(lapse['operational']):+.12e} "
        f"proxy={float(lapse['proxy']):.3e}",
        flush=True,
    )
print(f"\nOUTCOME: {outcome}", flush=True)
print(f"RESULT: {passed}/{tests}", flush=True)

if passed != tests:
    raise SystemExit(1)
