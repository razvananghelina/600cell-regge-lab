#!/usr/bin/env python3
"""Test transport of the homogeneous curvature-kernel line across two slabs.

Prior-art commit: ed51853.
Preregistered protocol commit: 8fe4176.
The target equation T_1 K_1 = K_2 is disclosed.
"""

import ast
from collections import Counter, defaultdict
import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
import json
import math
from pathlib import Path
import sys

from flint import acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la


HERE = Path(__file__).resolve().parent
TWO_INPUT = HERE / "gravity_600cell_dust_homogeneous_two_by_two.json"
TWO_SOURCE = HERE / "verify_gravity_600cell_dust_homogeneous_two_by_two.py"
FIRST_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
SECOND_INPUT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
SECOND_SOURCE = HERE / "verify_gravity_600cell_dust_second_tick_local_correction.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
CURVATURE_SOURCE = HERE / "verify_gravity_600cell_dust_internal_curvature_response.py"
OUTPUT = HERE / "gravity_600cell_dust_curvature_kernel_transport.json"

PRIOR_ART_COMMIT = "ed51853"
PROTOCOL_COMMIT = "8fe4176"
EXPECTED_HASHES = {
    "two": "d0017d4cfdf3a8833cf19bfcd287b21ac91a7f631c803d5d67114fdf64b77622",
    "two_source": "b97793c99ad2a24d5fd744f6a2e029b8fb51a40c632598c3860aeea602f6c816",
    "first": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "second_source": "cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "curvature_source": "276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66",
}

DPS = 100
BALL_DPS = 80
mp.mp.dps = DPS
ctx.dps = BALL_DPS
DERIVATIVE_STEPS = {
    "operational_primary": mp.mpf("1e-20"),
    "operational_shadow": mp.mpf("1e-15"),
    "validation_primary": mp.mpf("3e-20"),
    "validation_shadow": mp.mpf("3e-15"),
}
VARIANTS = tuple(DERIVATIVE_STEPS)
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_HINGES = tuple(combinations(range(5), 3))
LOCAL_HINGE_INDEX = {hinge: index for index, hinge in enumerate(LOCAL_HINGES)}
I = mp.mpc(0, 1)
ARITHMETIC_FLOOR = mp.mpf("1e-70")
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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_named_functions(source, wanted):
    tree = ast.parse(source.read_text(), filename=str(source))
    body = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in wanted]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch in {source.name}: {wanted-found}")
    exec(compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"), globals())


load_named_functions(RANK_SOURCE, {
    "orbit_sort_key", "augment_boundary_orbits", "log_minus",
    "signed_volume_square", "angle_record", "area_data",
    "extended_edge_image", "group_and_index_data", "prepare_geometry",
})
load_named_functions(TANGENT_SOURCE, {
    "mp_frobenius", "mp_submatrix", "cluster_sorted",
    "high_precision_sector_bases", "high_precision_pattern_cache",
    "assemble_full_representative_kernels", "project_full_kernel",
    "mp_to_acb", "mp_matrix_to_acb", "expanded_types", "build_tangent_ball",
})
load_named_functions(CURVATURE_SOURCE, {
    "triangle_area_square", "extended_triangle_image", "group_inverses",
    "triangle_response_data", "project_curvature_kernel", "singular_record",
})
load_named_functions(TWO_SOURCE, {
    "mp_fro", "mp_spectral_two_columns", "mp_singular_two_columns",
    "acb_number_mid_radius", "acb_matrix_mid_radius", "mp_to_numpy_high",
    "acb_uniform_basis", "response_ball", "boundary_mapping",
    "variant_uncertainty", "zero_label", "line_distance", "line_comparison",
    "sf", "sc",
})


def slab_index_data(model, a_old, a_new, r_value):
    data = dict(group_and_index_data(model, (mp.mpf(0), r_value)))
    rho = RHO0 * mp.exp(r_value)
    values = {
        "old": mp.exp(2 * a_old) * L0_SQUARE,
        "internal": mp.exp(a_old + a_new) * L0_SQUARE - rho,
        "pole": -rho,
        "new": mp.exp(2 * a_new) * L0_SQUARE,
    }
    if min(values["old"], values["internal"], rho, values["new"]) <= 0:
        raise RuntimeError("slab left the positive magnitude domain")
    data["rho"] = rho
    data["signed_base"] = tuple(values[kind] for kind in data["edge_kind"])
    return data, values


def reconstruct_slab(parity, label, state, need_tangent):
    a_old, a_new, r_value = state
    model = models[parity]
    index_data, kind_values = slab_index_data(model, a_old, a_new, r_value)
    geometry = prepare_geometry(model, index_data)
    mapping = boundary_mapping(index_data)
    sectors, sector_control = high_precision_sector_bases(index_data)
    trivial = [index for index, sector in enumerate(sectors)
               if sector["constant_overlap"] > mp.mpf("0.5")]
    sector_index = trivial[0] if len(trivial) == 1 else -1
    sector = sectors[sector_index]
    carrier_ok = bool(
        len(trivial) == 1 and sector["dimension"] == 1
        and sorted(mapping) == list(range(30))
        and len(geometry["triangle_records"]) == 6240
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
    )

    pattern_cache, branch = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    curvature_data = triangle_response_data(
        model, index_data, geometry, pattern_cache
    )
    kernels, hessian = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    blocks = {name: project_full_kernel(kernel, sector)
              for name, kernel in kernels.items()}
    curvature_blocks = {
        name: project_curvature_kernel(kernel, sector)
        for name, kernel in curvature_data["kernels"].items()
    }

    variant = {}
    determinants_ok = True
    for name in VARIANTS:
        det_response, y_ball = response_ball(blocks[name])
        determinants_ok &= not det_response.contains(0)
        tangent_ball = None
        if need_tangent:
            _, det_tangent, tangent_ball, _ = build_tangent_ball(
                blocks[name], 1, mapping
            )
            determinants_ok &= not det_tangent.contains(0)

        z_ball = acb_mat(95, 60)
        for index in range(30):
            z_ball[index, index] = 1
        for row in range(65):
            for column in range(60):
                z_ball[30 + row, column] = y_ball[row, column]
        d_ball = mp_matrix_to_acb(curvature_blocks[name])
        f_ball = d_ball * z_ball
        b_ball = f_ball * u_ball
        f_mid, f_radius = acb_matrix_mid_radius(f_ball)
        b_mid, b_radius = acb_matrix_mid_radius(b_ball)
        y_mid, y_radius = acb_matrix_mid_radius(y_ball)
        record = {
            "F_mp": f_mid,
            "F_radius": f_radius,
            "F": mp_to_numpy_high(f_mid),
            "D_norm": float(la.svdvals(mp_to_numpy_high(curvature_blocks[name]))[0]),
            "Z_radius": float(mp_fro(y_radius)),
            "B": b_mid,
            "B_radius": b_radius,
        }
        if need_tangent:
            record["T"], record["T_radius"] = acb_matrix_mid_radius(tangent_ball)
        variant[name] = record

    reconstruction_ok = bool(
        carrier_ok and branch["entry_pass"]
        and branch["base_negative_counts"] == Counter({1: 2400})
        and hessian["maximum_imaginary"] < ARITHMETIC_FLOOR
        and all(value < ARITHMETIC_FLOOR for value in
                curvature_data["maximum_derivative_imaginary"].values())
        and all(value < ARITHMETIC_FLOOR for value in
                curvature_data["maximum_equivariance_residual"].values())
        and determinants_ok
    )
    full_rank = singular_record(variant)
    rank_ok = bool(
        full_rank["resolved_rank"] == 59
        and full_rank["zero_count"] == 1
        and full_rank["open_count"] == 0
        and full_rank["columns"] == 60
    )
    matrices_b = {name: data["B"] for name, data in variant.items()}
    radii_b = {name: data["B_radius"] for name, data in variant.items()}
    epsilon_b = variant_uncertainty(matrices_b, radii_b)
    singulars = {}
    kernel = {}
    for name, data in variant.items():
        singulars[name], vectors = mp_singular_two_columns(data["B"])
        kernel[name] = vectors[:, 0] / mp_fro(vectors[:, 0])
    op_singular = singulars["operational_primary"]
    rank_one = bool(
        op_singular[0] <= 10 * epsilon_b
        and op_singular[1] > 100 * epsilon_b
    )
    gap = op_singular[1] - op_singular[0]
    kernel_error = epsilon_b / gap if gap > 100 * epsilon_b else mp.inf
    variant_distance = max(
        line_distance(kernel["operational_primary"], kernel[name])
        for name in VARIANTS
    )
    kernel_error += variant_distance + mp.mpf("1e-70")
    controls_ok = bool(reconstruction_ok and rank_ok and rank_one)
    return {
        "label": label,
        "state": state,
        "mapping": mapping,
        "old_orbits": index_data["orbit_edges"][:30],
        "controls_ok": controls_ok,
        "reconstruction_ok": reconstruction_ok,
        "rank_ok": rank_ok,
        "rank_one": rank_one,
        "full_rank": full_rank,
        "B_singular": op_singular,
        "epsilon_B": epsilon_b,
        "kernel": kernel,
        "kernel_error": kernel_error,
        "kernel_variant_distance": variant_distance,
        "variant": variant,
        "sector_index": sector_index,
    }


def transport_record(first_slab, second_slab):
    distances = {}
    curvature = {}
    transported = {}
    transported_errors = []
    for name in VARIANTS:
        first = first_slab["variant"][name]
        second = second_slab["variant"][name]
        k1 = first_slab["kernel"][name]
        k2 = second_slab["kernel"][name]
        v1 = u_mid * k1
        v2 = u_mid * k2
        w = first["T"] * v1
        transported[name] = w / mp_fro(w)
        distances[name] = line_distance(w, v2)
        curvature[name] = mp_fro(second["F_mp"] * w) / mp_fro(w)
        t_u = first["T"] * u_mid
        t_u_norm = mp_spectral_two_columns(t_u)
        t_radius = mp_fro(first["T_radius"])
        transported_errors.append(
            (t_u_norm * first_slab["kernel_error"] + t_radius)
            / max(mp.mpf("1e-80"), mp_fro(w))
        )

    op_distance = distances["operational_primary"]
    epsilon_step = (
        abs(op_distance - distances["operational_shadow"])
        + abs(distances["validation_primary"] - distances["validation_shadow"])
        + abs(op_distance - distances["validation_primary"])
    )
    transported_error = max(transported_errors)
    epsilon_line = (
        epsilon_step + 2 * transported_error + second_slab["kernel_error"]
        + mp.mpf("1e-70")
    )
    if epsilon_line >= mp.mpf("1e-2"):
        line_label = "NUMERICALLY_OPEN"
    elif op_distance <= 10 * epsilon_line:
        line_label = "TRANSPORT_IDENTIFIED"
    elif op_distance > 100 * epsilon_line:
        line_label = "TRANSPORT_SEPARATED"
    else:
        line_label = "NUMERICALLY_OPEN"

    op_curve = curvature["operational_primary"]
    epsilon_curve_step = (
        abs(op_curve - curvature["operational_shadow"])
        + abs(curvature["validation_primary"] - curvature["validation_shadow"])
        + abs(op_curve - curvature["validation_primary"])
    )
    f2_radius = max(mp_fro(data["F_radius"])
                    for data in second_slab["variant"].values())
    f2_norm = max(mp.mpf(la.svdvals(data["F"])[0])
                  for data in second_slab["variant"].values())
    epsilon_curve = (
        epsilon_curve_step + f2_radius
        + 2 * f2_norm * transported_error + mp.mpf("1e-70")
    )
    curve_label = zero_label(op_curve, epsilon_curve)
    consistent = not (
        (line_label == "TRANSPORT_IDENTIFIED" and curve_label == "NONZERO")
        or (line_label == "TRANSPORT_SEPARATED" and curve_label == "ZERO")
    )
    return {
        "distance": op_distance,
        "epsilon_distance": epsilon_line,
        "line_label": line_label,
        "curvature_response": op_curve,
        "epsilon_curvature_response": epsilon_curve,
        "curvature_label": curve_label,
        "consistent": consistent,
        "transported": transported,
        "transported_error": transported_error,
        "variant_distances": distances,
        "variant_curvature": curvature,
    }


def schedule_line_record(left, right, epsilon):
    distance, label = line_comparison(left, right, epsilon)
    return {"distance": distance, "epsilon": epsilon, "label": label}


hashes = {
    "two": sha256(TWO_INPUT),
    "two_source": sha256(TWO_SOURCE),
    "first": sha256(FIRST_INPUT),
    "second": sha256(SECOND_INPUT),
    "second_source": sha256(SECOND_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "tangent_source": sha256(TANGENT_SOURCE),
    "curvature_source": sha256(CURVATURE_SOURCE),
}
two = json.loads(TWO_INPUT.read_text())
first_tick = json.loads(FIRST_INPUT.read_text())
second_tick = json.loads(SECOND_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and two["outcome"] == "HOMOGENEOUS_2X2_KERNEL_NOT_EIGENLINE"
    and two["passed"] == two["tests"] == 17
    and first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["fixed_mass"] is True
    and second_tick["mass_recomputed_from_later_scale"] is False
)
check("all target-disclosed inputs have exact frozen provenance",
      provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_curvature_kernel_transport", GEOMETRY_SOURCE
)
gro = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gro
try:
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(gro)
except SystemExit as upstream_exit:
    if upstream_exit.code not in (None, 0):
        raise
check("the direct one-slab geometry import retains all 43 certificates",
      gro.tests == gro.passed == 43)

M_STAR = mp.mpf(10)
ZETA = (mp.pi**2 * mp.sqrt(2) / 50) ** (mp.mpf(1) / 3)
R0 = 4 * M_STAR / (3 * mp.pi)
L0 = ZETA * R0
L0_SQUARE = L0**2
EPSILON3 = 2 * mp.pi - 5 * mp.acos(mp.mpf(1) / 3)
MASS = (90 / mp.pi) * EPSILON3 * L0
RHO0 = mp.mpf("0.0102") ** 2
models = {parity: augment_boundary_orbits(model)
          for parity, model in gro.models.items()}
u_ball = acb_uniform_basis()
u_mid, u_radius = acb_matrix_mid_radius(u_ball)

seam_residuals = {}
seam_ok = True
for parity in ("even", "odd"):
    temporary_data, _ = slab_index_data(
        models[parity], mp.mpf(0), mp.mpf(first_tick["solutions"][parity]["state"][0]),
        mp.mpf(first_tick["solutions"][parity]["state"][1])
    )
    mapping = boundary_mapping(temporary_data)
    post = [mp.mpf(value) for value in first_tick["solutions"][parity]["post_momentum"]]
    pre = [mp.mpf(value) for value in second_tick["solutions"][parity]["pre_momentum"]]
    residual = max(abs(pre[index] - post[mapping[index]]) for index in range(30))
    bound = max(
        mp.mpf(first_tick["solutions"][parity]["junction_bound"]),
        mp.mpf(second_tick["solutions"][parity]["junction_bound"]),
    )
    seam_residuals[parity] = (residual, bound)
    seam_ok &= residual <= bound
check("the complete committed 30-component first-to-second seam passes",
      seam_ok,
      ", ".join(f"{p}: {sf(r,8)}/{sf(b,8)}" for p,(r,b) in seam_residuals.items()))

print("=" * 78)
print("TRANSPORT OF THE HOMOGENEOUS CURVATURE-KERNEL LINE")
print("=" * 78)

slabs = {}
transports = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43 and seam_ok
for parity in ("even", "odd"):
    a1, r1 = [mp.mpf(value) for value in first_tick["solutions"][parity]["state"]]
    a2, r2 = [mp.mpf(value) for value in
              second_tick["solutions"][parity]["state_absolute"]]
    first_state = (mp.mpf(0), a1, r1)
    second_state = (a1, a2, r2)
    print(f"[{parity}] reconstructing slab 1", flush=True)
    slab1 = reconstruct_slab(parity, "slab1", first_state, True)
    check(f"{parity}: slab 1 independently reproduces a unique homogeneous kernel",
          slab1["controls_ok"],
          f"rank={slab1['full_rank']['resolved_rank']}, "
          f"B=({sf(slab1['B_singular'][0],7)},{sf(slab1['B_singular'][1],7)})")
    print(f"[{parity}] reconstructing slab 2", flush=True)
    slab2 = reconstruct_slab(parity, "slab2", second_state, False)
    check(f"{parity}: slab 2 has a resolved unique homogeneous kernel",
          slab2["controls_ok"],
          f"rank={slab2['full_rank']['resolved_rank']}, "
          f"B=({sf(slab2['B_singular'][0],7)},{sf(slab2['B_singular'][1],7)})")
    transport = transport_record(slab1, slab2)
    check(f"{parity}: line transport and fresh curvature response agree",
          transport["consistent"],
          f"d={sf(transport['distance'],9)} ({transport['line_label']}), "
          f"F2w={sf(transport['curvature_response'],9)} "
          f"({transport['curvature_label']})")
    slabs[parity] = {"slab1": slab1, "slab2": slab2}
    transports[parity] = transport
    global_controls &= bool(slab1["controls_ok"] and slab2["controls_ok"]
                            and transport["consistent"])

even_orbits = slabs["even"]["slab1"]["old_orbits"]
odd_orbits = slabs["odd"]["slab1"]["old_orbits"]
even_to_odd = []
literal_map_ok = True
for even_orbit in even_orbits:
    matches = [index for index, odd_orbit in enumerate(odd_orbits)
               if frozenset(odd_orbit) == frozenset(even_orbit)]
    literal_map_ok &= len(matches) == 1
    if len(matches) == 1:
        even_to_odd.append(matches[0])
literal_map_ok &= sorted(even_to_odd) == list(range(30))
check("the target-independent even-to-odd orbit map reconstructs",
      literal_map_ok, str(even_to_odd))

schedule_records = {}
for name, getter, error_getter in (
    ("K1", lambda p: slabs[p]["slab1"]["kernel"]["operational_primary"],
     lambda p: slabs[p]["slab1"]["kernel_error"]),
    ("K2", lambda p: slabs[p]["slab2"]["kernel"]["operational_primary"],
     lambda p: slabs[p]["slab2"]["kernel_error"]),
    ("T1K1", lambda p: transports[p]["transported"]["operational_primary"],
     lambda p: transports[p]["transported_error"]),
):
    epsilon = error_getter("even") + error_getter("odd") + mp.mpf("1e-70")
    schedule_records[name] = schedule_line_record(
        getter("even"), getter("odd"), epsilon
    )
check("K1, K2 and T1K1 receive complete calibrated schedule labels",
      all(record["label"] in {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"}
          for record in schedule_records.values()),
      str({name: record["label"] for name, record in schedule_records.items()}))

transport_labels = {record["line_label"] for record in transports.values()}
curve_labels = {record["curvature_label"] for record in transports.values()}
schedule_separated = any(record["label"] == "SEPARATED"
                         for record in schedule_records.values())
if not global_controls or not literal_map_ok:
    outcome = "CURVATURE_KERNEL_TRANSPORT_CONTROL_FAILED"
elif schedule_separated or len(transport_labels) > 1 or len(curve_labels) > 1:
    outcome = "CURVATURE_KERNEL_TRANSPORT_SCHEDULE_DEPENDENT"
elif transport_labels == {"TRANSPORT_IDENTIFIED"} and curve_labels == {"ZERO"}:
    outcome = "CURVATURE_KERNEL_TRANSPORT_IDENTIFIED"
elif transport_labels == {"TRANSPORT_SEPARATED"} and curve_labels == {"NONZERO"}:
    outcome = "CURVATURE_KERNEL_TRANSPORT_REFUTED"
else:
    outcome = "CURVATURE_KERNEL_TRANSPORT_NUMERICALLY_OPEN"

allowed = {
    "CURVATURE_KERNEL_TRANSPORT_CONTROL_FAILED",
    "CURVATURE_KERNEL_TRANSPORT_SCHEDULE_DEPENDENT",
    "CURVATURE_KERNEL_TRANSPORT_IDENTIFIED",
    "CURVATURE_KERNEL_TRANSPORT_REFUTED",
    "CURVATURE_KERNEL_TRANSPORT_NUMERICALLY_OPEN",
}
check("the preregistered hierarchy assigns the transport outcome",
      outcome in allowed, f"outcome={outcome}")

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": "T_1 K_1 = K_2",
    "fixed_mass": True,
    "outcome": outcome,
    "literal_even_to_odd_orbit_map": even_to_odd,
    "seam": {
        parity: {"residual": sf(value[0]), "bound": sf(value[1])}
        for parity, value in seam_residuals.items()
    },
    "schedule_comparisons": {
        name: {"distance": sf(record["distance"]),
               "epsilon": sf(record["epsilon"]), "label": record["label"]}
        for name, record in schedule_records.items()
    },
    "parities": {
        parity: {
            "slab1_state": [sf(value, 60) for value in slabs[parity]["slab1"]["state"]],
            "slab2_state": [sf(value, 60) for value in slabs[parity]["slab2"]["state"]],
            "slab1_full_rank": slabs[parity]["slab1"]["full_rank"]["resolved_rank"],
            "slab2_full_rank": slabs[parity]["slab2"]["full_rank"]["resolved_rank"],
            "slab1_B_singular": [sf(value) for value in slabs[parity]["slab1"]["B_singular"]],
            "slab2_B_singular": [sf(value) for value in slabs[parity]["slab2"]["B_singular"]],
            "slab1_kernel_qp": [sc(slabs[parity]["slab1"]["kernel"]["operational_primary"][i]) for i in range(2)],
            "slab2_kernel_qp": [sc(slabs[parity]["slab2"]["kernel"]["operational_primary"][i]) for i in range(2)],
            "transport_distance": sf(transports[parity]["distance"]),
            "epsilon_transport_distance": sf(transports[parity]["epsilon_distance"]),
            "transport_label": transports[parity]["line_label"],
            "fresh_F2_response": sf(transports[parity]["curvature_response"]),
            "epsilon_fresh_F2_response": sf(transports[parity]["epsilon_curvature_response"]),
            "fresh_F2_label": transports[parity]["curvature_label"],
            "transported_line_error": sf(transports[parity]["transported_error"]),
        }
        for parity in ("even", "odd")
    },
    "classification": {
        "line_transport": "DERIVED COMPUTATIONAL",
        "constraint_or_gauge": "OPEN",
        "later_ticks": "OPEN",
        "refinement": "OPEN",
        "external_novelty": "OPEN",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
