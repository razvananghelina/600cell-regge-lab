#!/usr/bin/env python3
"""Blind full second-slab tangent and two-step cocycle census.

Prior-art commit: d51fbbf.
Preregistered protocol commit: 4ec04b1.
No continuum, particle, speed or stability target is parsed.
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
import zipfile

from flint import arb, acb, acb_mat, ctx
import mpmath as mp
import numpy as np
import scipy.linalg as la
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
FIRST_FULL_INPUT = HERE / "gravity_600cell_dust_full_boundary_tangent.json"
FIRST_FULL_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
FIRST_TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
SECOND_TICK_INPUT = HERE / "gravity_600cell_dust_second_tick_local_correction.json"
SECOND_TICK_SOURCE = HERE / "verify_gravity_600cell_dust_second_tick_local_correction.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
OUTPUT = HERE / "gravity_600cell_dust_two_step_full_tangent.json"
NUMERIC_OUTPUT = HERE / "gravity_600cell_dust_two_step_full_tangent.npz"

PRIOR_ART_COMMIT = "d51fbbf"
PROTOCOL_COMMIT = "4ec04b1"
EXPECTED_HASHES = {
    "first_full": "4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5",
    "first_full_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "first_tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
    "second_tick": "936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70",
    "second_tick_source": "cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
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
    "norm2", "mp_frobenius", "mp_submatrix", "cluster_sorted",
    "high_precision_sector_bases", "high_precision_pattern_cache",
    "assemble_full_representative_kernels", "project_full_kernel",
    "mp_to_numpy", "mp_to_acb", "mp_matrix_to_acb",
    "acb_midpoint_and_radii", "expanded_types", "build_tangent_ball",
    "optimal_spectral_distance", "tangent_analysis", "serialize_float",
    "serialize_mp", "serialize_complex", "deterministic_npz",
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
        raise RuntimeError("second slab left the positive magnitude domain")
    data["rho"] = rho
    data["signed_base"] = tuple(values[kind] for kind in data["edge_kind"])
    return data, values


def boundary_mapping(index_data):
    mapping = []
    for old_type in range(30):
        shifted = tuple(tuple(vertex + 120 for vertex in edge)
                        for edge in index_data["orbit_edges"][old_type])
        matches = [final for final in range(30)
                   if shifted == index_data["orbit_edges"][65 + final]]
        if len(matches) != 1:
            raise RuntimeError(f"boundary mapping is not unique: {old_type}, {matches}")
        mapping.append(matches[0])
    return tuple(mapping)


def radius_ball(midpoint, radius):
    center = arb(format(float(midpoint), ".17g"))
    if radius <= 0:
        return center
    return center + arb(0, format(float(radius), ".17g"))


def reenclose_binary_matrix(midpoint, radii):
    rows, columns = midpoint.shape
    matrix = acb_mat(rows, columns)
    for row in range(rows):
        for column in range(columns):
            value = complex(midpoint[row, column])
            stored = float(radii[row, column])
            real_radius = stored + 0.5 * abs(float(np.spacing(value.real)))
            imag_radius = stored + 0.5 * abs(float(np.spacing(value.imag)))
            matrix[row, column] = acb(
                radius_ball(value.real, real_radius),
                radius_ball(value.imag, imag_radius),
            )
    return matrix


def omega_ball(size):
    matrix = acb_mat(2 * size, 2 * size)
    for index in range(size):
        matrix[index, size + index] = 1
        matrix[size + index, index] = -1
    return matrix


def map_branch_counts(analysis, multiplicity):
    values = analysis["eigen_arrays"]["operational_primary"]
    epsilon = analysis["epsilon_eigenvalue"]
    distances = np.abs(np.abs(values) - 1)
    if math.isfinite(epsilon):
        unit = distances < 10 * epsilon
        resolved = distances > 100 * epsilon
    else:
        unit = np.zeros(len(values), dtype=bool)
        resolved = np.zeros(len(values), dtype=bool)
    open_flags = ~(unit | resolved)
    contracting = resolved & (np.abs(values) < 1)
    expanding = resolved & (np.abs(values) > 1)
    resolved_moduli = np.abs(values[resolved])
    return {
        "minimal_dimension": len(values),
        "unit_minimal": int(np.sum(unit)),
        "open_minimal": int(np.sum(open_flags)),
        "resolved_minimal": int(np.sum(resolved)),
        "contracting_minimal": int(np.sum(contracting)),
        "expanding_minimal": int(np.sum(expanding)),
        "unit_weighted": int(multiplicity * np.sum(unit)),
        "open_weighted": int(multiplicity * np.sum(open_flags)),
        "resolved_weighted": int(multiplicity * np.sum(resolved)),
        "contracting_weighted": int(multiplicity * np.sum(contracting)),
        "expanding_weighted": int(multiplicity * np.sum(expanding)),
        "minimum_resolved_modulus": (
            float(np.min(resolved_moduli)) if len(resolved_moduli) else None
        ),
        "maximum_resolved_modulus": (
            float(np.max(resolved_moduli)) if len(resolved_moduli) else None
        ),
    }


def schedule_label(distance, epsilon):
    if not math.isfinite(epsilon):
        return "SCHEDULE_OPEN"
    if distance <= 10 * epsilon:
        return "SCHEDULE_ROBUST"
    if distance > 100 * epsilon:
        return "SCHEDULE_DEPENDENT"
    return "SCHEDULE_OPEN"


def schedule_record(left, right):
    singular_distance = float(np.max(np.abs(
        left["singular_arrays"]["operational_primary"]
        - right["singular_arrays"]["operational_primary"]
    )))
    singular_epsilon = float(
        left["epsilon_singular"] + right["epsilon_singular"]
        + 10 * np.finfo(float).eps
        * max(1.0, left["singular_arrays"]["operational_primary"][0],
              right["singular_arrays"]["operational_primary"][0])
    )
    eigen_distance = optimal_spectral_distance(
        left["eigen_arrays"]["operational_primary"],
        right["eigen_arrays"]["operational_primary"],
    )
    eigen_epsilon = float(
        left["epsilon_eigenvalue"] + right["epsilon_eigenvalue"]
    )
    return {
        "singular_distance": singular_distance,
        "singular_epsilon": singular_epsilon,
        "singular_label": schedule_label(singular_distance, singular_epsilon),
        "eigen_distance": eigen_distance,
        "eigen_epsilon": eigen_epsilon,
        "eigen_label": schedule_label(eigen_distance, eigen_epsilon),
    }


def public_analysis(analysis, counts):
    return {
        "epsilon_tangent": serialize_float(analysis["epsilon_t"]),
        "maximum_tangent_ball_radius": serialize_float(
            analysis["maximum_tangent_ball_radius"]
        ),
        "maximum_defect_ball_radius": serialize_float(
            analysis["maximum_defect_ball_radius"]
        ),
        "symplectic_norm": serialize_float(analysis["symplectic_norm"]),
        "epsilon_symplectic": serialize_float(analysis["epsilon_sym"]),
        "symplectic_ok": analysis["symplectic_ok"],
        "reciprocal_singular_norm": serialize_float(analysis["reciprocal_norm"]),
        "epsilon_reciprocal_singular": serialize_float(analysis["epsilon_reciprocal"]),
        "reciprocal_svd_floor": serialize_float(analysis["reciprocal_svd_floor"]),
        "reciprocal_singular_ok": analysis["reciprocal_ok"],
        "determinant_log_modulus": serialize_float(
            analysis["determinant_log_moduli"]["operational_primary"]
        ),
        "epsilon_determinant_log_modulus": serialize_float(
            analysis["epsilon_determinant_log_modulus"]
        ),
        "determinant_ok": analysis["determinant_ok"],
        "tangent_condition": serialize_float(analysis["condition_tangent"]),
        "eigenvector_condition": serialize_float(analysis["eigenvector_condition"]),
        "epsilon_eigenvalue": serialize_float(analysis["epsilon_eigenvalue"]),
        "reciprocal_conjugate_eigenvalue_distance": serialize_float(
            analysis["reciprocal_eigenvalue_distance"]
        ),
        "singular_values": [serialize_float(value) for value in
                            analysis["singular_arrays"]["operational_primary"]],
        "eigenvalues": [serialize_complex(value) for value in
                         analysis["eigen_arrays"]["operational_primary"]],
        "counts": counts,
    }


hashes = {
    "first_full": sha256(FIRST_FULL_INPUT),
    "first_full_numeric": sha256(FIRST_FULL_NUMERIC),
    "tangent_source": sha256(TANGENT_SOURCE),
    "first_tick": sha256(FIRST_TICK_INPUT),
    "second_tick": sha256(SECOND_TICK_INPUT),
    "second_tick_source": sha256(SECOND_TICK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
}
first_full = json.loads(FIRST_FULL_INPUT.read_text())
first_numeric = np.load(FIRST_FULL_NUMERIC)
first_tick = json.loads(FIRST_TICK_INPUT.read_text())
second_tick = json.loads(SECOND_TICK_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and first_full["outcome"] == "FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED"
    and first_full["numeric_archive_arrays"] == len(first_numeric.files) == 224
    and first_tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
    and second_tick["outcome"] == "SECOND_HOMOTHETIC_TICK_ACCEPTED"
    and second_tick["fixed_mass"] is True
    and second_tick["mass_recomputed_from_later_scale"] is False
)
check("all blind inputs have exact frozen provenance", provenance_ok, str(hashes))

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_two_step_full_tangent", GEOMETRY_SOURCE
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

seam_ok = True
seam_records = {}
for parity in ("even", "odd"):
    a1, r1 = [mp.mpf(value) for value in first_tick["solutions"][parity]["state"]]
    temporary, _ = slab_index_data(models[parity], 0, a1, r1)
    mapping = boundary_mapping(temporary)
    post = [mp.mpf(value) for value in first_tick["solutions"][parity]["post_momentum"]]
    pre = [mp.mpf(value) for value in second_tick["solutions"][parity]["pre_momentum"]]
    residual = max(abs(pre[index] - post[mapping[index]]) for index in range(30))
    bound = max(mp.mpf(first_tick["solutions"][parity]["junction_bound"]),
                mp.mpf(second_tick["solutions"][parity]["junction_bound"]))
    seam_records[parity] = {"residual": residual, "bound": bound}
    seam_ok &= residual <= bound
check("both complete 30-component canonical seams pass", seam_ok,
      str({p: (mp.nstr(v['residual'], 5), mp.nstr(v['bound'], 5))
           for p, v in seam_records.items()}))

print("=" * 78)
print("BLIND FULL SECOND-SLAB TANGENT AND TWO-STEP COCYCLE")
print("=" * 78)

records = {}
numeric_arrays = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43 and seam_ok
common_sector_signature = None

for parity in ("even", "odd"):
    a1, _ = [mp.mpf(value) for value in first_tick["solutions"][parity]["state"]]
    a2, r2 = [mp.mpf(value) for value in
              second_tick["solutions"][parity]["state_absolute"]]
    index_data, kind_values = slab_index_data(models[parity], a1, a2, r2)
    geometry = prepare_geometry(models[parity], index_data)
    mapping = boundary_mapping(index_data)
    carrier_ok = bool(
        len(models[parity]["old_edges"]) == 720
        and len(models[parity]["internal_edges"]) == 840
        and len(models[parity]["new_edges"]) == 720
        and len(geometry["triangle_records"]) == 6240
        and sorted(mapping) == list(range(30))
    )
    sectors, sector_control = high_precision_sector_bases(index_data)
    signature = tuple((sector["dimension"],
                       mp.nstr(sector["old_central_eigenvalue"], 70),
                       sector["splitter"]) for sector in sectors)
    if common_sector_signature is None:
        common_sector_signature = signature
    basis_ok = bool(
        sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
        and sum(60 * sector["dimension"] ** 2 for sector in sectors) == 1440
        and signature == common_sector_signature
        and all(sector["dimension"] ==
                first_full["parities"][parity]["sectors"][index]["dimension"]
                for index, sector in enumerate(sectors))
    )
    check(f"{parity}: carrier and seven deterministic sectors exhaust 1440 dimensions",
          carrier_ok and basis_ok,
          f"dims={[sector['dimension'] for sector in sectors]}, map={mapping}")

    print(f"[{parity}] differentiating second-slab local patterns", flush=True)
    pattern_cache, branch = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    kernel_ok = bool(
        branch["entry_pass"]
        and branch["base_negative_counts"] == Counter({1: 2400})
        and branch["minimum_leading_minor"] > 0
        and branch["minimum_argument"] > mp.mpf("1e-6")
        and kernel_control["maximum_imaginary"] < ARITHMETIC_FLOOR
        and len(set(kernel_control["nonzero_entries"].values())) == 1
    )
    check(f"{parity}: second-slab branch, derivative and reality controls pass",
          kernel_ok, f"entries={kernel_control['nonzero_entries']}")

    sector_records = []
    all_determinants = True
    all_t2_canonical = True
    all_product_canonical = True
    for sector_index, sector in enumerate(sectors):
        dimension = sector["dimension"]
        print(f"[{parity}] sector {sector_index+1}/7 d={dimension}", flush=True)
        blocks = {name: project_full_kernel(kernel, sector)
                  for name, kernel in kernels.items()}
        t2_ball_records = {}
        product_ball_records = {}
        determinant_flags = {}
        for name, block in blocks.items():
            _, determinant, t2_ball, t2_defect_ball = build_tangent_ball(
                block, dimension, mapping
            )
            determinant_flags[name] = not determinant.contains(0)
            prefix1 = f"{parity}_sector{sector_index}_{name}"
            t1_mid = first_numeric[f"{prefix1}_tangent_midpoint"]
            t1_rad = first_numeric[f"{prefix1}_tangent_radii"]
            t1_ball = reenclose_binary_matrix(t1_mid, t1_rad)
            product_ball = t2_ball * t1_ball
            omega = omega_ball(30 * dimension)
            product_defect_ball = (
                product_ball.transpose().conjugate() * omega * product_ball - omega
            )

            map_balls = {
                "t2": (t2_ball, t2_defect_ball),
                "product": (product_ball, product_defect_ball),
            }
            for map_name, (map_ball, defect_ball) in map_balls.items():
                midpoint, radii = acb_midpoint_and_radii(map_ball)
                defect_midpoint, defect_radii = acb_midpoint_and_radii(defect_ball)
                target = t2_ball_records if map_name == "t2" else product_ball_records
                target[name] = {
                    "midpoint": midpoint,
                    "radii": radii,
                    "defect_midpoint": defect_midpoint,
                    "defect_radii": defect_radii,
                    "det_j": determinant,
                }
                prefix = f"{parity}_sector{sector_index}_{map_name}_{name}"
                numeric_arrays[f"{prefix}_midpoint"] = midpoint
                numeric_arrays[f"{prefix}_radii"] = radii
                numeric_arrays[f"{prefix}_defect_midpoint"] = defect_midpoint
                numeric_arrays[f"{prefix}_defect_radii"] = defect_radii

        determinants_ok = all(determinant_flags.values())
        t2_analysis = tangent_analysis(t2_ball_records)
        product_analysis = tangent_analysis(product_ball_records)
        for analysis, ball_records in (
            (t2_analysis, t2_ball_records),
            (product_analysis, product_ball_records),
        ):
            analysis["maximum_tangent_ball_radius"] = max(
                la.norm(item["radii"], "fro") for item in ball_records.values()
            )
            analysis["maximum_defect_ball_radius"] = max(
                la.norm(item["defect_radii"], "fro")
                for item in ball_records.values()
            )
        all_determinants &= determinants_ok
        all_t2_canonical &= t2_analysis["canonicality_ok"]
        all_product_canonical &= product_analysis["canonicality_ok"]
        sector_records.append({
            "sector_index": sector_index,
            "dimension": dimension,
            "center": sector["old_central_eigenvalue"],
            "splitter": sector["splitter"],
            "constant_overlap": sector["constant_overlap"],
            "determinant_flags": determinant_flags,
            "t2": t2_analysis,
            "product": product_analysis,
            "t2_counts": map_branch_counts(t2_analysis, dimension),
            "product_counts": map_branch_counts(product_analysis, dimension),
        })

    check(f"{parity}: all 28 second-slab pre-Legendre determinant balls exclude zero",
          all_determinants)
    check(f"{parity}: all seven second-slab tangent blocks are canonical",
          all_t2_canonical)
    check(f"{parity}: all seven rigorous two-step product blocks are canonical",
          all_product_canonical)
    controls_ok = bool(carrier_ok and basis_ok and kernel_ok and all_determinants
                       and all_t2_canonical and all_product_canonical)
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "mapping": mapping,
        "state": (a1, a2, r2),
        "sectors": sector_records,
    }

archive_ok = len(numeric_arrays) == 448
deterministic_npz(NUMERIC_OUTPUT, numeric_arrays)
check("the deterministic two-map archive contains exactly 448 arrays",
      archive_ok, f"arrays={len(numeric_arrays)}, sha={sha256(NUMERIC_OUTPUT)}")

schedule = []
for sector_index in range(7):
    left = records["even"]["sectors"][sector_index]
    right = records["odd"]["sectors"][sector_index]
    for map_name in ("t2", "product"):
        comparison = schedule_record(left[map_name], right[map_name])
        schedule.append({
            "sector_index": sector_index,
            "dimension": left["dimension"],
            "map": map_name,
            **comparison,
        })
primary_counts = Counter(item["singular_label"] for item in schedule)
eigen_counts = Counter(item["eigen_label"] for item in schedule)
complete_schedule = len(schedule) == 14 and sum(primary_counts.values()) == 14
check("all fourteen blind schedule comparisons receive calibrated labels",
      complete_schedule,
      f"singular={dict(primary_counts)}, eigen={dict(eigen_counts)}")

resolved_dependent = any(
    item["singular_label"] == "SCHEDULE_DEPENDENT"
    or item["eigen_label"] == "SCHEDULE_DEPENDENT"
    for item in schedule
)
primary_open = any(item["singular_label"] == "SCHEDULE_OPEN" for item in schedule)
if not global_controls or not archive_ok or not complete_schedule:
    outcome = "TWO_STEP_FULL_TANGENT_CONTROL_FAILED"
elif resolved_dependent:
    outcome = "TWO_STEP_FULL_TANGENT_SCHEDULE_DEPENDENT"
elif primary_open:
    outcome = "TWO_STEP_FULL_TANGENT_SCHEDULE_OPEN"
else:
    outcome = "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED"

allowed = {
    "TWO_STEP_FULL_TANGENT_CONTROL_FAILED",
    "TWO_STEP_FULL_TANGENT_SCHEDULE_DEPENDENT",
    "TWO_STEP_FULL_TANGENT_SCHEDULE_OPEN",
    "TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED",
}
check("the frozen hierarchy assigns the blind cocycle outcome",
      outcome in allowed, f"outcome={outcome}")

def aggregate_counts(parity, map_name):
    fields = ("unit_weighted", "open_weighted", "resolved_weighted",
              "contracting_weighted", "expanding_weighted")
    result = {field: sum(sector[f"{map_name}_counts"][field]
                         for sector in records[parity]["sectors"])
              for field in fields}
    moduli = [sector[f"{map_name}_counts"][key]
              for sector in records[parity]["sectors"]
              for key in ("minimum_resolved_modulus", "maximum_resolved_modulus")
              if sector[f"{map_name}_counts"][key] is not None]
    result["minimum_resolved_modulus"] = min(moduli) if moduli else None
    result["maximum_resolved_modulus"] = max(moduli) if moduli else None
    result["total"] = result["unit_weighted"] + result["open_weighted"] + result["resolved_weighted"]
    return result

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "blind": True,
    "target_comparisons_performed": False,
    "outcome": outcome,
    "full_phase_dimension": 1440,
    "numeric_archive": NUMERIC_OUTPUT.name,
    "numeric_archive_arrays": len(numeric_arrays),
    "numeric_archive_sha256": sha256(NUMERIC_OUTPUT),
    "seam": {parity: {key: serialize_mp(value)
                       for key, value in record.items()}
             for parity, record in seam_records.items()},
    "schedule_label_counts": {
        "singular": dict(primary_counts), "eigenvalue": dict(eigen_counts)
    },
    "schedule_comparisons": [
        {key: (serialize_float(value) if isinstance(value, float) else value)
         for key, value in item.items()}
        for item in schedule
    ],
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "state_a1_a2_r2": [serialize_mp(value) for value in record["state"]],
            "boundary_mapping": list(record["mapping"]),
            "t2_full_counts": aggregate_counts(parity, "t2"),
            "product_full_counts": aggregate_counts(parity, "product"),
            "sectors": [{
                "sector_index": sector["sector_index"],
                "dimension": sector["dimension"],
                "central_eigenvalue": serialize_complex(complex(
                    float(mp.re(sector["center"])), float(mp.im(sector["center"]))
                )),
                "splitter_group_index": sector["splitter"],
                "constant_overlap": serialize_mp(sector["constant_overlap"]),
                "pre_legendre_determinants_exclude_zero": sector["determinant_flags"],
                "t2": public_analysis(sector["t2"], sector["t2_counts"]),
                "product": public_analysis(sector["product"], sector["product_counts"]),
            } for sector in record["sectors"]],
        }
        for parity, record in records.items()
    },
    "classification": {
        "finite_canonical_cocycle": "DERIVED COMPUTATIONAL",
        "euclidean_amplification_norm": "STRUCTURAL",
        "physical_mode_interpretation": "OPEN",
        "continuum_or_refinement": "OPEN",
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
print(f"Numeric archive: {NUMERIC_OUTPUT}")
if passed != tests:
    raise SystemExit(1)
