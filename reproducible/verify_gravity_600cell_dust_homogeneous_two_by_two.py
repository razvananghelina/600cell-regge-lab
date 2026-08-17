#!/usr/bin/env python3
"""Calibrate the homogeneous 2x2 tangent and curvature-kernel line.

Prior-art commit: 24eed99.
Preregistered protocol commit: f139f60.
The uniform basis and the previously observed near-minus target are disclosed.
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
KERNEL_INPUT = HERE / "gravity_600cell_dust_homogeneous_curvature_kernel.json"
KERNEL_SOURCE = HERE / "verify_gravity_600cell_dust_homogeneous_curvature_kernel.py"
CURVATURE_INPUT = HERE / "gravity_600cell_dust_internal_curvature_response.json"
CURVATURE_SOURCE = HERE / "verify_gravity_600cell_dust_internal_curvature_response.py"
TANGENT_NUMERIC = HERE / "gravity_600cell_dust_full_boundary_tangent.npz"
TANGENT_SOURCE = HERE / "verify_gravity_600cell_dust_full_boundary_tangent.py"
ALIGNMENT_SOURCE = HERE / "verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py"
RANK_SOURCE = HERE / "verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py"
GEOMETRY_SOURCE = HERE / "verify_gravity_global_regge_orbits.py"
TICK_INPUT = HERE / "gravity_600cell_dust_homothetic_canonical_lapse.json"
OUTPUT = HERE / "gravity_600cell_dust_homogeneous_two_by_two.json"

PRIOR_ART_COMMIT = "24eed99"
PROTOCOL_COMMIT = "f139f60"
EXPECTED_HASHES = {
    "kernel": "b55887ff3905afd94e86821852d58f0d60c227b52dfbd945044874bfe87540e9",
    "kernel_source": "43837b4d97fcf21cc6de9e4debea0c22bc827d5186d0c75ca07dfe5c799e1a15",
    "curvature": "95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc",
    "curvature_source": "276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66",
    "tangent_numeric": "816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b",
    "tangent_source": "c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571",
    "alignment_source": "e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047",
    "rank_source": "834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5",
    "geometry_source": "ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf",
    "tick": "4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9",
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
    body = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in wanted
    ]
    found = {node.name for node in body}
    if found != wanted:
        raise RuntimeError(f"audited function mismatch in {source.name}: {wanted-found}")
    exec(
        compile(ast.Module(body=body, type_ignores=[]), str(source), "exec"),
        globals(),
    )


def mp_fro(matrix):
    return mp.sqrt(mp.fsum(
        abs(matrix[row, column]) ** 2
        for row in range(matrix.rows) for column in range(matrix.cols)
    ))


def mp_spectral_two_columns(matrix):
    gram = matrix.H * matrix
    values, _ = mp.eighe(gram)
    largest = max(mp.mpf(0), mp.re(values[values.rows - 1]))
    return mp.sqrt(largest)


def mp_singular_two_columns(matrix):
    gram = matrix.H * matrix
    values, vectors = mp.eighe(gram)
    values = [max(mp.mpf(0), mp.re(values[index])) for index in range(2)]
    return (mp.sqrt(values[0]), mp.sqrt(values[1])), vectors


def acb_number_mid_radius(value):
    real = mp.mpf(value.real.mid().str(90, radius=False))
    imag = mp.mpf(value.imag.mid().str(90, radius=False))
    real_radius = mp.mpf(value.real.rad().upper().str(90, radius=False))
    imag_radius = mp.mpf(value.imag.rad().upper().str(90, radius=False))
    return mp.mpc(real, imag), mp.sqrt(real_radius**2 + imag_radius**2)


def acb_matrix_mid_radius(matrix):
    midpoint = mp.matrix(matrix.nrows(), matrix.ncols())
    radii = mp.matrix(matrix.nrows(), matrix.ncols())
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            midpoint[row, column], radii[row, column] = acb_number_mid_radius(
                matrix[row, column]
            )
    return midpoint, radii


def mp_to_numpy_high(matrix):
    return np.array([
        [complex(float(mp.re(matrix[row, column])), float(mp.im(matrix[row, column])))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ], dtype=np.complex128)


def acb_identity(size):
    matrix = acb_mat(size, size)
    for index in range(size):
        matrix[index, index] = 1
    return matrix


def acb_uniform_basis():
    coefficient = mp.mpf(1) / mp.sqrt(30)
    coefficient_ball = acb(mp.nstr(coefficient, 90))
    basis = acb_mat(60, 2)
    for index in range(30):
        basis[index, 0] = coefficient_ball
        basis[30 + index, 1] = coefficient_ball
    return basis


def response_ball(block):
    old = expanded_types(0, 30, 1)
    internal = expanded_types(30, 65, 1)
    new = expanded_types(65, 95, 1)
    k_xx = mp_submatrix(block, internal, internal)
    k_xn = mp_submatrix(block, internal, new)
    k_ox = mp_submatrix(block, old, internal)
    k_on = mp_submatrix(block, old, new)
    j_matrix = mp.matrix(65, 65)
    for row in range(35):
        for column in range(35):
            j_matrix[row, column] = k_xx[row, column]
        for column in range(30):
            j_matrix[row, 35 + column] = k_xn[row, column]
    for row in range(30):
        for column in range(35):
            j_matrix[35 + row, column] = -k_ox[row, column]
        for column in range(30):
            j_matrix[35 + row, 35 + column] = -k_on[row, column]
    rhs = mp.matrix(65, 60)
    k_xo = mp_submatrix(block, internal, old)
    k_oo = mp_submatrix(block, old, old)
    for row in range(35):
        for column in range(30):
            rhs[row, column] = -k_xo[row, column]
    for row in range(30):
        for column in range(30):
            rhs[35 + row, column] = k_oo[row, column]
        rhs[35 + row, 30 + row] = 1
    j_ball = mp_matrix_to_acb(j_matrix)
    return j_ball.det(), j_ball.solve(mp_matrix_to_acb(rhs))


def boundary_mapping(index_data):
    mapping = []
    for old_type in range(30):
        shifted = tuple(
            tuple(vertex + 120 for vertex in edge)
            for edge in index_data["orbit_edges"][old_type]
        )
        matches = [
            final_type for final_type in range(30)
            if shifted == index_data["orbit_edges"][65 + final_type]
        ]
        if len(matches) != 1:
            raise RuntimeError(f"boundary orbit map is not unique: {old_type}, {matches}")
        mapping.append(matches[0])
    return tuple(mapping)


def variant_uncertainty(matrices, radii):
    op = matrices["operational_primary"]
    return (
        mp_fro(op - matrices["operational_shadow"])
        + mp_fro(matrices["validation_primary"]
                 - matrices["validation_shadow"])
        + mp_fro(op - matrices["validation_primary"])
        + max(mp_fro(value) for value in radii.values())
        + mp.mpf("1e-70")
    )


def zero_label(value, epsilon, bounded=False):
    if bounded and epsilon >= mp.mpf("1e-2"):
        return "NUMERICALLY_OPEN"
    if value <= 10 * epsilon:
        return "ZERO"
    if value > 100 * epsilon:
        return "NONZERO"
    return "NUMERICALLY_OPEN"


def line_distance(left, right):
    left = left / mp_fro(left)
    right = right / mp_fro(right)
    overlap = abs((left.H * right)[0])
    return mp.sqrt(max(mp.mpf(0), 1 - min(mp.mpf(1), overlap) ** 2))


def line_comparison(left, right, epsilon):
    distance = line_distance(left, right)
    if epsilon >= mp.mpf("1e-2"):
        label = "NUMERICALLY_OPEN"
    elif distance <= 10 * epsilon:
        label = "IDENTIFIED"
    elif distance > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return distance, label


def matrix_comparison(left, right, epsilon):
    distance = mp_fro(left - right)
    if distance <= 10 * epsilon:
        label = "IDENTIFIED"
    elif distance > 100 * epsilon:
        label = "SEPARATED"
    else:
        label = "NUMERICALLY_OPEN"
    return distance, label


def sf(value, digits=40):
    return mp.nstr(value, digits)


def sc(value, digits=40):
    return {"real": sf(mp.re(value), digits), "imaginary": sf(mp.im(value), digits)}


def sm(matrix, digits=40):
    return [[sc(matrix[row, column], digits) for column in range(matrix.cols)]
            for row in range(matrix.rows)]


hashes = {
    "kernel": sha256(KERNEL_INPUT),
    "kernel_source": sha256(KERNEL_SOURCE),
    "curvature": sha256(CURVATURE_INPUT),
    "curvature_source": sha256(CURVATURE_SOURCE),
    "tangent_numeric": sha256(TANGENT_NUMERIC),
    "tangent_source": sha256(TANGENT_SOURCE),
    "alignment_source": sha256(ALIGNMENT_SOURCE),
    "rank_source": sha256(RANK_SOURCE),
    "geometry_source": sha256(GEOMETRY_SOURCE),
    "tick": sha256(TICK_INPUT),
}
kernel_input = json.loads(KERNEL_INPUT.read_text())
curvature_input = json.loads(CURVATURE_INPUT.read_text())
tick = json.loads(TICK_INPUT.read_text())
provenance_ok = bool(
    hashes == EXPECTED_HASHES
    and kernel_input["outcome"] == "HOMOGENEOUS_CURVATURE_KERNEL_SUBSPACE_LOCALIZED"
    and kernel_input["passed"] == kernel_input["tests"] == 14
    and curvature_input["outcome"] == "STRONG_TANGENT_CURVATURE_INJECTIVE"
    and tick["outcome"] == "HOMOTHETIC_CANONICAL_LAPSE_SELECTED"
)
check("all target-disclosed inputs have exact frozen provenance", provenance_ok, str(hashes))


load_named_functions(RANK_SOURCE, {
    "orbit_sort_key", "augment_boundary_orbits", "log_minus",
    "signed_volume_square", "angle_record", "area_data",
    "extended_edge_image", "group_and_index_data", "prepare_geometry",
})
load_named_functions(TANGENT_SOURCE, {
    "mp_submatrix", "cluster_sorted", "high_precision_sector_bases",
    "high_precision_pattern_cache", "assemble_full_representative_kernels",
    "project_full_kernel", "mp_to_acb", "mp_matrix_to_acb", "expanded_types",
    "build_tangent_ball",
})
load_named_functions(CURVATURE_SOURCE, {
    "triangle_area_square", "extended_triangle_image", "group_inverses",
    "triangle_response_data", "project_curvature_kernel", "singular_record",
})

spec = importlib.util.spec_from_file_location(
    "global_regge_orbits_homogeneous_two_by_two", GEOMETRY_SOURCE
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
u_orthogonality = mp_fro(u_mid.H * u_mid - mp.eye(2))
check("the disclosed uniform q/p basis is orthonormal at the arithmetic floor",
      u_orthogonality < ARITHMETIC_FLOOR and mp_fro(u_radius) < ARITHMETIC_FLOOR,
      f"orthogonality={sf(u_orthogonality, 8)}")

print("=" * 78)
print("HOMOGENEOUS 2x2 TANGENT AND INTERNAL-CURVATURE KERNEL")
print("=" * 78)

records = {}
old_boundary_orderings = {}
global_controls = provenance_ok and gro.tests == gro.passed == 43

for parity in ("even", "odd"):
    print(f"[{parity}] reconstructing high-precision trivial-sector operators", flush=True)
    model = models[parity]
    state = tick["solutions"][parity]["state"]
    index_data = group_and_index_data(model, state)
    geometry = prepare_geometry(model, index_data)
    mapping = boundary_mapping(index_data)
    sectors, sector_control = high_precision_sector_bases(index_data)
    trivial_indices = [index for index, sector in enumerate(sectors)
                       if sector["constant_overlap"] > mp.mpf("0.5")]
    trivial_ok = len(trivial_indices) == 1
    sector_index = trivial_indices[0]
    sector = sectors[sector_index]
    old_boundary_orderings[parity] = index_data["orbit_edges"][:30]
    carrier_ok = bool(
        trivial_ok and sector["dimension"] == 1
        and sorted(mapping) == list(range(30))
        and len(geometry["triangle_records"]) == 6240
        and sector_control["irrep_dimensions"] == [1, 1, 1, 2, 2, 2, 3]
    )
    check(f"{parity}: carrier, trivial sector and literal boundary map reconstruct",
          carrier_ok, f"sector={sector_index}, mapping={mapping}")

    s = mp.mpf(state[0])
    kind_values = {
        "old": L0_SQUARE,
        "internal": mp.exp(s) * L0_SQUARE - index_data["rho"],
        "pole": -index_data["rho"],
        "new": mp.exp(2 * s) * L0_SQUARE,
    }
    pattern_cache, branch_control = high_precision_pattern_cache(
        geometry["patterns"], kind_values
    )
    curvature_data = triangle_response_data(
        model, index_data, geometry, pattern_cache
    )
    hessian_kernels, hessian_control = assemble_full_representative_kernels(
        index_data, geometry, pattern_cache
    )
    blocks = {name: project_full_kernel(kernel, sector)
              for name, kernel in hessian_kernels.items()}
    curvature_blocks = {
        name: project_curvature_kernel(kernel, sector)
        for name, kernel in curvature_data["kernels"].items()
    }

    variant = {}
    determinants_ok = True
    for name in VARIANTS:
        det_response, y_ball = response_ball(blocks[name])
        _, det_tangent, tangent_ball, _ = build_tangent_ball(
            blocks[name], 1, mapping
        )
        determinants_ok &= bool(
            not det_response.contains(0) and not det_tangent.contains(0)
        )

        z_ball = acb_mat(95, 60)
        for index in range(30):
            z_ball[index, index] = 1
        for row in range(65):
            for column in range(60):
                z_ball[30 + row, column] = y_ball[row, column]
        d_ball = mp_matrix_to_acb(curvature_blocks[name])
        f_ball = d_ball * z_ball

        projector = acb_identity(60) - u_ball * u_ball.transpose().conjugate()
        a_ball = u_ball.transpose().conjugate() * tangent_ball * u_ball
        l_ball = projector * tangent_ball * u_ball
        b_ball = f_ball * u_ball

        tangent_mid, tangent_radius = acb_matrix_mid_radius(tangent_ball)
        f_mid, f_radius = acb_matrix_mid_radius(f_ball)
        a_mid, a_radius = acb_matrix_mid_radius(a_ball)
        l_mid, l_radius = acb_matrix_mid_radius(l_ball)
        b_mid, b_radius = acb_matrix_mid_radius(b_ball)
        y_mid, y_radius = acb_matrix_mid_radius(y_ball)
        d_numpy = mp_to_numpy_high(curvature_blocks[name])
        variant[name] = {
            "T": tangent_mid, "T_radius": tangent_radius,
            "F_mp": f_mid, "F_radius": f_radius,
            "F": mp_to_numpy_high(f_mid),
            "D_norm": float(la.svdvals(d_numpy)[0]),
            "Z_radius": float(mp_fro(y_radius)),
            "A": a_mid, "A_radius": a_radius,
            "L": l_mid, "L_radius": l_radius,
            "B": b_mid, "B_radius": b_radius,
        }

    reconstruction_ok = bool(
        branch_control["entry_pass"]
        and hessian_control["maximum_imaginary"] < ARITHMETIC_FLOOR
        and all(value < ARITHMETIC_FLOOR for value in
                curvature_data["maximum_derivative_imaginary"].values())
        and all(value < ARITHMETIC_FLOOR for value in
                curvature_data["maximum_equivariance_residual"].values())
        and determinants_ok
    )
    check(f"{parity}: branch, reality, equivariance and ball determinants pass",
          reconstruction_ok)

    full_rank = singular_record(variant)
    rank_ok = bool(
        full_rank["resolved_rank"] == 59
        and full_rank["zero_count"] == 1
        and full_rank["open_count"] == 0
        and full_rank["columns"] == 60
    )
    check(f"{parity}: direct response reproduces calibrated rank 59/nullity 1",
          rank_ok,
          f"eps={full_rank['epsilon_singular']:.3e}, "
          f"smin={full_rank['minimum_singular']:.3e}")

    matrices_a = {name: data["A"] for name, data in variant.items()}
    radii_a = {name: data["A_radius"] for name, data in variant.items()}
    matrices_l = {name: data["L"] for name, data in variant.items()}
    radii_l = {name: data["L_radius"] for name, data in variant.items()}
    matrices_b = {name: data["B"] for name, data in variant.items()}
    radii_b = {name: data["B_radius"] for name, data in variant.items()}
    epsilon_a = variant_uncertainty(matrices_a, radii_a)
    epsilon_l = variant_uncertainty(matrices_l, radii_l)
    epsilon_b = variant_uncertainty(matrices_b, radii_b)

    leakage = mp_spectral_two_columns(matrices_l["operational_primary"])
    leakage_label = zero_label(leakage, epsilon_l)
    singulars = {}
    kernels = {}
    for name, data in variant.items():
        singulars[name], vectors = mp_singular_two_columns(data["B"])
        kernels[name] = vectors[:, 0] / mp_fro(vectors[:, 0])
    op_singular = singulars["operational_primary"]
    b_rank_one = bool(
        op_singular[0] <= 10 * epsilon_b
        and op_singular[1] > 100 * epsilon_b
    )
    singular_gap = op_singular[1] - op_singular[0]
    kernel_error = (
        epsilon_b / singular_gap if singular_gap > 100 * epsilon_b else mp.inf
    )
    kernel_variant_distance = max(
        line_distance(kernels["operational_primary"], kernels[name])
        for name in VARIANTS
    )
    kernel_error += kernel_variant_distance + mp.mpf("1e-70")
    check(f"{parity}: plane leakage and 160x2 rank receive frozen labels",
          leakage_label in {"ZERO", "NONZERO", "NUMERICALLY_OPEN"}
          and (b_rank_one or op_singular[0] > 10 * epsilon_b
               or op_singular[1] <= 100 * epsilon_b),
          f"leak={sf(leakage, 8)} ({leakage_label}), "
          f"B singular=({sf(op_singular[0], 8)},{sf(op_singular[1], 8)}), "
          f"epsB={sf(epsilon_b, 8)}")

    residuals = {}
    multipliers = {}
    tu_norms = {}
    for name, data in variant.items():
        k = kernels[name]
        vector = u_mid * k
        image = data["T"] * vector
        mu = (vector.H * image)[0]
        residual = image - vector * mu
        residuals[name] = mp_fro(residual)
        multipliers[name] = mu
        tu_norms[name] = mp_spectral_two_columns(data["T"] * u_mid)

    op_residual = residuals["operational_primary"]
    epsilon_residual_step = (
        abs(op_residual - residuals["operational_shadow"])
        + abs(residuals["validation_primary"]
              - residuals["validation_shadow"])
        + abs(op_residual - residuals["validation_primary"])
    )
    epsilon_tu_ball = max(mp_fro(data["T_radius"])
                          for data in variant.values())
    epsilon_residual = (
        epsilon_residual_step + epsilon_tu_ball
        + 4 * max(tu_norms.values()) * kernel_error + mp.mpf("1e-70")
    )
    eigenline_label = zero_label(op_residual, epsilon_residual, bounded=True)

    op_mu = multipliers["operational_primary"]
    epsilon_mu_step = (
        abs(op_mu - multipliers["operational_shadow"])
        + abs(multipliers["validation_primary"]
              - multipliers["validation_shadow"])
        + abs(op_mu - multipliers["validation_primary"])
    )
    epsilon_mu = (
        epsilon_mu_step + epsilon_tu_ball
        + 2 * max(tu_norms.values()) * kernel_error + mp.mpf("1e-70")
    )
    minus_distance = abs(op_mu + 1)
    minus_zero = zero_label(minus_distance, epsilon_mu)
    minus_label = {
        "ZERO": "EXACT_MINUS_ONE_WITHIN_ERROR",
        "NONZERO": "RESOLVED_NOT_MINUS_ONE",
        "NUMERICALLY_OPEN": "NUMERICALLY_OPEN",
    }[minus_zero]

    a_op = matrices_a["operational_primary"]
    trace_a = a_op[0, 0] + a_op[1, 1]
    determinant_a = a_op[0, 0] * a_op[1, 1] - a_op[0, 1] * a_op[1, 0]
    discriminant = trace_a**2 - 4 * determinant_a
    root = mp.sqrt(discriminant)
    eigenvalues_a = ((trace_a - root) / 2, (trace_a + root) / 2)
    eigenvalues_a = tuple(sorted(eigenvalues_a, key=lambda value: abs(value)))
    symplectic = mp.matrix([[0, 1], [-1, 0]])
    symplectic_defect = mp_fro(a_op.H * symplectic * a_op - symplectic)
    diagnostics_ok = bool(
        eigenline_label in {"ZERO", "NONZERO", "NUMERICALLY_OPEN"}
        and minus_label in {
            "EXACT_MINUS_ONE_WITHIN_ERROR", "RESOLVED_NOT_MINUS_ONE",
            "NUMERICALLY_OPEN",
        }
        and all(mp.isfinite(value) for value in (
            leakage, epsilon_l, epsilon_b, op_residual, epsilon_residual,
            minus_distance, epsilon_mu, symplectic_defect,
        ))
    )
    check(f"{parity}: 2x2 eigenline, multiplier and symplectic diagnostics calibrate",
          diagnostics_ok,
          f"res={sf(op_residual, 8)} ({eigenline_label}), "
          f"mu={sf(op_mu, 14)}, mu+1={sf(minus_distance, 8)} ({minus_label})")

    controls_ok = bool(carrier_ok and reconstruction_ok and rank_ok and diagnostics_ok)
    global_controls &= controls_ok
    records[parity] = {
        "controls_ok": controls_ok,
        "sector_index": sector_index,
        "mapping": mapping,
        "variant": variant,
        "A": a_op,
        "epsilon_A": epsilon_a,
        "plane_leakage": leakage,
        "epsilon_plane_leakage": epsilon_l,
        "plane_label": leakage_label,
        "B_singular": op_singular,
        "epsilon_B": epsilon_b,
        "B_rank_one": b_rank_one,
        "kernel": kernels["operational_primary"],
        "kernel_error": kernel_error,
        "kernel_variant_distance": kernel_variant_distance,
        "residual": op_residual,
        "epsilon_residual": epsilon_residual,
        "eigenline_label": eigenline_label,
        "mu": op_mu,
        "epsilon_mu": epsilon_mu,
        "minus_distance": minus_distance,
        "minus_label": minus_label,
        "trace": trace_a,
        "determinant": determinant_a,
        "eigenvalues": eigenvalues_a,
        "symplectic_defect": symplectic_defect,
        "full_rank": full_rank,
    }


even_orbits = old_boundary_orderings["even"]
odd_orbits = old_boundary_orderings["odd"]
even_to_odd = []
literal_map_ok = True
for even_orbit in even_orbits:
    matches = [index for index, odd_orbit in enumerate(odd_orbits)
               if frozenset(odd_orbit) == frozenset(even_orbit)]
    literal_map_ok &= len(matches) == 1
    if len(matches) == 1:
        even_to_odd.append(matches[0])
literal_map_ok &= sorted(even_to_odd) == list(range(30))
uniform_fixed = bool(
    literal_map_ok
    and all(abs(u_mid[index, 0] - u_mid[even_to_odd[index], 0])
            < ARITHMETIC_FLOOR for index in range(30))
)
check("the literal even-to-odd orbit permutation fixes both uniform vectors",
      uniform_fixed, str(even_to_odd))

kernel_schedule_epsilon = (
    records["even"]["kernel_error"] + records["odd"]["kernel_error"]
    + mp.mpf("1e-70")
)
kernel_schedule_distance, kernel_schedule_label = line_comparison(
    records["even"]["kernel"], records["odd"]["kernel"],
    kernel_schedule_epsilon,
)
check("the two 2x2 curvature-kernel lines receive a calibrated schedule label",
      kernel_schedule_label in {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"},
      f"distance={sf(kernel_schedule_distance, 8)}, "
      f"epsilon={sf(kernel_schedule_epsilon, 8)}, label={kernel_schedule_label}")

a_schedule_epsilon = (
    records["even"]["epsilon_A"] + records["odd"]["epsilon_A"]
    + mp.mpf("1e-70")
)
a_schedule_distance, a_schedule_label = matrix_comparison(
    records["even"]["A"], records["odd"]["A"], a_schedule_epsilon
)
check("the two compressed tangent matrices receive a calibrated schedule label",
      a_schedule_label in {"IDENTIFIED", "SEPARATED", "NUMERICALLY_OPEN"},
      f"distance={sf(a_schedule_distance, 8)}, "
      f"epsilon={sf(a_schedule_epsilon, 8)}, label={a_schedule_label}")


if not global_controls:
    outcome = "HOMOGENEOUS_2X2_CONTROL_FAILED"
elif "SEPARATED" in {kernel_schedule_label, a_schedule_label}:
    outcome = "HOMOGENEOUS_2X2_SCHEDULE_DEPENDENT"
elif not all(record["B_rank_one"] for record in records.values()):
    outcome = "HOMOGENEOUS_2X2_CURVATURE_KERNEL_OPEN"
elif any(record["eigenline_label"] == "NONZERO" for record in records.values()):
    outcome = "HOMOGENEOUS_2X2_KERNEL_NOT_EIGENLINE"
elif all(record["eigenline_label"] == "ZERO" for record in records.values()):
    plane_labels = {record["plane_label"] for record in records.values()}
    if "NONZERO" in plane_labels:
        outcome = "HOMOGENEOUS_2X2_EIGENLINE_PLANE_LEAKS"
    elif plane_labels == {"ZERO"}:
        outcome = "HOMOGENEOUS_2X2_INVARIANT_EIGENLINE"
    else:
        outcome = "HOMOGENEOUS_2X2_EIGENLINE_PLANE_OPEN"
else:
    outcome = "HOMOGENEOUS_2X2_NUMERICALLY_OPEN"

allowed_outcomes = {
    "HOMOGENEOUS_2X2_CONTROL_FAILED",
    "HOMOGENEOUS_2X2_SCHEDULE_DEPENDENT",
    "HOMOGENEOUS_2X2_CURVATURE_KERNEL_OPEN",
    "HOMOGENEOUS_2X2_KERNEL_NOT_EIGENLINE",
    "HOMOGENEOUS_2X2_EIGENLINE_PLANE_LEAKS",
    "HOMOGENEOUS_2X2_EIGENLINE_PLANE_OPEN",
    "HOMOGENEOUS_2X2_INVARIANT_EIGENLINE",
    "HOMOGENEOUS_2X2_NUMERICALLY_OPEN",
}
check("the preregistered hierarchy assigns the 2x2 outcome",
      outcome in allowed_outcomes, f"outcome={outcome}")

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "target_disclosed": True,
    "uniform_basis": "(1_30,0_30)/sqrt(30), (0_30,1_30)/sqrt(30)",
    "outcome": outcome,
    "literal_even_to_odd_orbit_map": even_to_odd,
    "schedule_comparison": {
        "kernel_line_distance": sf(kernel_schedule_distance),
        "kernel_line_epsilon": sf(kernel_schedule_epsilon),
        "kernel_line_label": kernel_schedule_label,
        "compressed_matrix_distance": sf(a_schedule_distance),
        "compressed_matrix_epsilon": sf(a_schedule_epsilon),
        "compressed_matrix_label": a_schedule_label,
    },
    "parities": {
        parity: {
            "controls_ok": record["controls_ok"],
            "trivial_sector_index": record["sector_index"],
            "boundary_mapping": list(record["mapping"]),
            "full_response_rank": record["full_rank"]["resolved_rank"],
            "full_response_zero_count": record["full_rank"]["zero_count"],
            "compressed_tangent": sm(record["A"]),
            "epsilon_compressed_tangent": sf(record["epsilon_A"]),
            "trace": sc(record["trace"]),
            "determinant": sc(record["determinant"]),
            "eigenvalues": [sc(value) for value in record["eigenvalues"]],
            "compressed_symplectic_defect": sf(record["symplectic_defect"]),
            "plane_leakage": sf(record["plane_leakage"]),
            "epsilon_plane_leakage": sf(record["epsilon_plane_leakage"]),
            "plane_label": record["plane_label"],
            "curvature_singular_values": [sf(value) for value in record["B_singular"]],
            "epsilon_curvature_singular": sf(record["epsilon_B"]),
            "curvature_rank_one": record["B_rank_one"],
            "kernel_vector_qp": [sc(record["kernel"][index]) for index in range(2)],
            "kernel_angular_error": sf(record["kernel_error"]),
            "kernel_variant_distance": sf(record["kernel_variant_distance"]),
            "full_eigenline_residual": sf(record["residual"]),
            "epsilon_full_eigenline_residual": sf(record["epsilon_residual"]),
            "eigenline_label": record["eigenline_label"],
            "rayleigh_multiplier": sc(record["mu"]),
            "epsilon_multiplier": sf(record["epsilon_mu"]),
            "distance_from_minus_one": sf(record["minus_distance"]),
            "minus_one_label": record["minus_label"],
        }
        for parity, record in records.items()
    },
    "classification": {
        "finite_linear_algebra": "DERIVED COMPUTATIONAL",
        "minus_one_test": "CONFIRMATORY TARGET-DISCLOSED",
        "gauge_lapse_or_time": "OPEN",
        "nonlinear_integrability": "OPEN",
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
