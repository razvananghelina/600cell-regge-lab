#!/usr/bin/env python3
"""Focused diagnostic for the first complete refined internal Hessian run.

Protocol commit: cdd0f69.
Only lexicographic schedule 0 is assembled/factorized.  No physical target,
spatial mode, continuum value, or desired nullity is loaded.
"""

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


SOURCE = HERE / "verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py"
FROZEN = HERE / "gravity_600cell_refined_nonhomogeneous_internal_hessian.json"
CURVATURE = HERE / "gravity_600cell_refined_local_curvature_mass.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_refined_h4_stationary_fill.py"
PROTOCOL = ROOT / "docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_diagnostic_protocol.md"
OUTPUT = HERE / "gravity_600cell_refined_nonhomogeneous_internal_hessian_diagnostic.json"

PROTOCOL_COMMIT = "cdd0f69"
EXPECTED_HASHES = {
    "reproducible/verify_gravity_600cell_refined_nonhomogeneous_internal_hessian.py":
        "2a96c8ce466d6d4e9be2cadc4ebf932b4e42eff16fc64fbfb08cd580d680879e",
    "reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian.json":
        "4a05968c68f8e6a35a1308ddf6114bb19b7106f214bfdcf798e7af2387bddec1",
    "reproducible/gravity_600cell_refined_local_curvature_mass.json":
        "180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091",
    "docs/gravity/gravity_600cell_refined_nonhomogeneous_internal_hessian_diagnostic_protocol.md":
        "6781e48da86cf79b853956477388d180374c023c2ee07652fb010fe1595587cf",
}

PAIR4 = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR4)}
LOCAL_EDGES = tuple(combinations(range(5), 2))
LOCAL_TRIANGLES = tuple(combinations(range(5), 3))
TRIANGLE_EDGES = tuple(combinations(range(3), 2))
TAU_TEXT = "0.0102"
DECIMAL_PRECISIONS = (100, 140)
DIFFERENCE_STEPS = ("1e-10", "5e-11", "2.5e-11")
EXPECTED_F = (2640, 17040, 28800, 14400)
EXPECTED_COUNTS = {
    "pentachora": 57600,
    "triangles": 149280,
    "boundary_edges": 34080,
    "internal_edges": 19680,
    "cross_edges": 17040,
    "vertical_edges": 2640,
    "all_edges": 53760,
}
FORWARD_INDICES = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14)
EXPECTED_OLD_ZERO = (6, 7, 12, 14)
NONZERO_GRADIENT_THRESHOLD = mp.mpf("1e-20")

tests = 0
passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return ok


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return format(float(value), ".17e")


def ldtext(value):
    return np.format_float_scientific(
        np.longdouble(value), precision=20, unique=False, trim="k"
    )


def clong(value):
    return (
        np.clongdouble(np.longdouble(mp.nstr(mp.re(value), 80)))
        + np.clongdouble(1j)
        * np.clongdouble(np.longdouble(mp.nstr(mp.im(value), 80)))
    )


def load_source_definitions():
    """Load definitions only; never execute the frozen verifier top level."""
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = {
        "ast": ast,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "sha256": sha256,
        "combinations": combinations,
        "permutations": permutations,
        "json": json,
        "Path": Path,
        "sys": sys,
        "mp": mp,
        "np": np,
        "sp": sp,
        "spla": spla,
        "build_600cell": build_600cell,
        "HERE": HERE,
        "ROOT": ROOT,
        "ACTION_SOURCE": ACTION_SOURCE,
        "PAIR4": PAIR4,
        "PAIR_INDEX": PAIR_INDEX,
        "LOCAL_EDGES": LOCAL_EDGES,
        "LOCAL_TRIANGLES": LOCAL_TRIANGLES,
        "TRIANGLE_EDGES": TRIANGLE_EDGES,
        "TAU_TEXT": TAU_TEXT,
        "DECIMAL_PRECISIONS": DECIMAL_PRECISIONS,
        "DIFFERENCE_STEPS": DIFFERENCE_STEPS,
        "EXPECTED_F": EXPECTED_F,
        "EXPECTED_COUNTS": EXPECTED_COUNTS,
        "tests": 0,
        "passed": 0,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(SOURCE), "exec"), namespace)
    return namespace


def angle_caches(source, actions, kinds, geometries):
    """Build the frozen binary cache and a long-double cache with envelopes."""
    binary = []
    extended = []
    for index, pattern in enumerate(kinds):
        print(f"  local angle pattern {index + 1}/{len(kinds)}", flush=True)
        low = source["angle_stencil_at_precision"](
            actions, pattern, geometries[100], 100
        )
        high = source["angle_stencil_at_precision"](
            actions, pattern, geometries[140], 140
        )
        selected = np.empty((10, 10), dtype=np.complex128)
        derivative_ld = np.empty((10, 10), dtype=np.clongdouble)
        derivative_error = np.empty((10, 10), dtype=np.longdouble)
        envelope = np.empty((10, 10), dtype=np.float64)
        for row in range(10):
            for column in range(10):
                value = high["fine"][row, column]
                error = 100 * (
                    abs(high["fine"][row, column]
                        - high["coarse"][row, column])
                    + abs(high["fine"][row, column]
                          - low["fine"][row, column])
                    + mp.mpf("1e-70")
                )
                selected[row, column] = complex(value)
                derivative_ld[row, column] = clong(value)
                derivative_error[row, column] = np.longdouble(mp.nstr(error, 80))
                envelope[row, column] = float(error)
        base = np.asarray([complex(value) for value in high["base"]])
        base_ld = np.asarray([clong(value) for value in high["base"]])
        base_error = np.asarray([
            np.longdouble(mp.nstr(
                100 * (abs(high["base"][row] - low["base"][row])
                       + mp.mpf("1e-70")), 80
            ))
            for row in range(10)
        ])
        binary.append({
            "base": base,
            "derivative": selected,
            "envelope": envelope,
        })
        extended.append({
            "base": base_ld,
            "base_error": base_error,
            "derivative": derivative_ld,
            "derivative_error": derivative_error,
        })
    return tuple(binary), tuple(extended)


def area_caches(source, kinds, geometries):
    binary = []
    extended = []
    for pattern in kinds:
        with mp.workdps(100):
            _, low_gradient, low_hessian = source["area_data"]([
                source["kind_value"](kind, geometries[100]) for kind in pattern
            ])
        with mp.workdps(140):
            _, high_gradient, high_hessian = source["area_data"]([
                source["kind_value"](kind, geometries[140]) for kind in pattern
            ])
        gradient = np.asarray([complex(value) for value in high_gradient])
        hessian = np.asarray([
            [complex(value) for value in row] for row in high_hessian
        ])
        gradient_ld = np.asarray([clong(value) for value in high_gradient])
        hessian_ld = np.asarray([
            [clong(value) for value in row] for row in high_hessian
        ])
        gradient_error = np.asarray([
            np.longdouble(mp.nstr(
                100 * (abs(high_gradient[row] - low_gradient[row])
                       + mp.mpf("1e-70")), 80
            )) for row in range(3)
        ])
        hessian_error = np.asarray([[
            np.longdouble(mp.nstr(
                100 * (abs(high_hessian[row][column]
                           - low_hessian[row][column]) + mp.mpf("1e-70")), 80
            )) for column in range(3)
        ] for row in range(3)])
        binary.append({"gradient": gradient, "hessian": hessian})
        extended.append({
            "gradient": gradient_ld,
            "gradient_error": gradient_error,
            "hessian": hessian_ld,
            "hessian_error": hessian_error,
        })
    return tuple(binary), tuple(extended)


def grouped_add(matrix, errors, absolute, counts, row_orbits, column_orbits,
                values, value_errors):
    codes = row_orbits.astype(np.int64) * 10 + column_orbits.astype(np.int64)
    flat = matrix.ravel()
    flat_error = errors.ravel()
    flat_absolute = absolute.ravel()
    flat_counts = counts.ravel()
    np.add.at(flat, codes, values)
    np.add.at(flat_error, codes, value_errors)
    np.add.at(flat_absolute, codes, np.abs(values))
    np.add.at(flat_counts, codes, 1)


def direct_longdouble_pullback(schedule, angle_cache, area_cache,
                               per_vertex_masses):
    """Sum actual incidences directly into orbit entries, not via sparse C."""
    simplex_pattern = schedule["simplex_pattern_global"]
    triangle_pattern = schedule["triangle_pattern_global"]
    simplex_triangles = schedule["simplex_triangles"]
    triangle_internal = schedule["triangle_internal"]
    simplex_internal = schedule["simplex_internal"]
    internal_orbits = schedule["internal_orbits"]

    angle_base = np.asarray([record["base"] for record in angle_cache])
    angle_base_error = np.asarray([
        record["base_error"] for record in angle_cache
    ])
    angle_derivative = np.asarray([
        record["derivative"] for record in angle_cache
    ])
    angle_derivative_error = np.asarray([
        record["derivative_error"] for record in angle_cache
    ])
    area_gradient = np.asarray([record["gradient"] for record in area_cache])
    area_gradient_error = np.asarray([
        record["gradient_error"] for record in area_cache
    ])
    area_hessian = np.asarray([record["hessian"] for record in area_cache])
    area_hessian_error = np.asarray([
        record["hessian_error"] for record in area_cache
    ])

    pi_ld = clong(mp.pi)
    curvature = np.where(
        schedule["triangle_boundary"], pi_ld, 2 * pi_ld
    ).astype(np.clongdouble)
    curvature_error = np.full(
        len(curvature), np.longdouble(4) * np.finfo(np.longdouble).eps
        * np.abs(pi_ld), dtype=np.longdouble
    )
    simplex_angles = angle_base[simplex_pattern]
    simplex_angle_errors = angle_base_error[simplex_pattern]
    np.add.at(curvature, simplex_triangles.ravel(), simplex_angles.ravel())
    np.add.at(
        curvature_error,
        simplex_triangles.ravel(),
        simplex_angle_errors.ravel(),
    )

    matrix = np.zeros((10, 10), dtype=np.clongdouble)
    errors = np.zeros((10, 10), dtype=np.longdouble)
    absolute = np.zeros((10, 10), dtype=np.longdouble)
    counts = np.zeros((10, 10), dtype=np.int64)

    coefficient = -np.clongdouble(1j) * curvature
    coefficient_error = curvature_error
    tri_grad = area_gradient[triangle_pattern]
    tri_grad_error = area_gradient_error[triangle_pattern]
    tri_hess = area_hessian[triangle_pattern]
    tri_hess_error = area_hessian_error[triangle_pattern]

    for row_position in range(3):
        rows = triangle_internal[:, row_position]
        for column_position in range(3):
            columns = triangle_internal[:, column_position]
            valid = (rows >= 0) & (columns >= 0)
            values = coefficient[valid] * tri_hess[
                valid, row_position, column_position
            ]
            local_error = (
                coefficient_error[valid]
                * np.abs(tri_hess[valid, row_position, column_position])
                + np.abs(coefficient[valid])
                * tri_hess_error[valid, row_position, column_position]
                + coefficient_error[valid]
                * tri_hess_error[valid, row_position, column_position]
            )
            grouped_add(
                matrix, errors, absolute, counts,
                internal_orbits[rows[valid]], internal_orbits[columns[valid]],
                values, local_error,
            )

    for hinge in range(10):
        triangle_indices = simplex_triangles[:, hinge]
        rows_all = triangle_internal[triangle_indices]
        derivative = angle_derivative[simplex_pattern, hinge]
        derivative_error = angle_derivative_error[simplex_pattern, hinge]
        for row_position in range(3):
            rows = rows_all[:, row_position]
            gradient = tri_grad[triangle_indices, row_position]
            gradient_error = tri_grad_error[triangle_indices, row_position]
            for column_position in range(10):
                columns = simplex_internal[:, column_position]
                valid = (rows >= 0) & (columns >= 0)
                values = (
                    -np.clongdouble(1j) * gradient[valid]
                    * derivative[valid, column_position]
                )
                local_error = (
                    gradient_error[valid]
                    * np.abs(derivative[valid, column_position])
                    + np.abs(gradient[valid])
                    * derivative_error[valid, column_position]
                    + gradient_error[valid]
                    * derivative_error[valid, column_position]
                )
                grouped_add(
                    matrix, errors, absolute, counts,
                    internal_orbits[rows[valid]],
                    internal_orbits[columns[valid]],
                    values, local_error,
                )

    tau = np.longdouble(TAU_TEXT)
    vertical_indices = np.flatnonzero(schedule["internal_is_vertical"])
    ranks = schedule["internal_orbits"][vertical_indices] - 6
    for rank in range(4):
        population = int(np.count_nonzero(ranks == rank))
        mass = np.longdouble(mp.nstr(per_vertex_masses[rank], 80))
        value = -2 * np.longdouble(mp.nstr(mp.pi, 80)) * mass * tau
        matrix[6 + rank, 6 + rank] += population * value
        absolute[6 + rank, 6 + rank] += population * abs(value)
        counts[6 + rank, 6 + rank] += population

    unit_roundoff = np.finfo(np.longdouble).eps / 2
    gamma = np.zeros((10, 10), dtype=np.longdouble)
    valid_counts = counts > 0
    operations = counts[valid_counts] + 32
    gamma[valid_counts] = (
        operations * unit_roundoff / (1 - operations * unit_roundoff)
    )
    errors += (gamma + 4 * unit_roundoff) * absolute

    real = matrix.real
    real = (real + real.T) / 2
    symmetric_error = (errors + errors.T) / 2
    return real, symmetric_error, counts, {
        "maximum_imaginary": ldtext(np.max(np.abs(matrix.imag))),
        "maximum_entry_error": ldtext(np.max(symmetric_error)),
        "maximum_term_count": int(counts.max()),
    }


def matrix_comparison(left, right, left_error, right_error):
    difference = np.asarray(left, dtype=np.longdouble) - np.asarray(
        right, dtype=np.longdouble
    )
    absolute = np.abs(difference)
    position = tuple(map(int, np.unravel_index(np.argmax(absolute), absolute.shape)))
    maximum = absolute[position]
    if np.ndim(left_error) == 0:
        left_bounds = np.full(absolute.shape, np.longdouble(left_error))
    else:
        left_bounds = np.asarray(left_error, dtype=np.longdouble)
    if np.ndim(right_error) == 0:
        right_bounds = np.full(absolute.shape, np.longdouble(right_error))
    else:
        right_bounds = np.asarray(right_error, dtype=np.longdouble)
    gates = np.longdouble(100) * (left_bounds + right_bounds)
    fractions = np.divide(
        absolute,
        gates,
        out=np.full_like(absolute, np.longdouble(np.inf)),
        where=gates > 0,
    )
    fraction_position = tuple(map(
        int, np.unravel_index(np.argmax(fractions), fractions.shape)
    ))
    return {
        "maximum_absolute_difference": ldtext(maximum),
        "maximum_position": list(position),
        "signed_difference_at_maximum": ldtext(difference[position]),
        "gate_at_maximum": ldtext(gates[position]),
        "maximum_gate_fraction": ldtext(fractions[fraction_position]),
        "maximum_gate_fraction_position": list(fraction_position),
        "matches": bool(np.all(absolute <= gates)),
    }


def select_corruption(schedule, area_cache, threshold=None):
    for simplex in range(len(schedule["slab"])):
        for hinge in range(10):
            triangle = int(schedule["simplex_triangles"][simplex, hinge])
            tri_pattern = int(schedule["triangle_pattern_global"][triangle])
            for row_position in range(3):
                row = int(schedule["triangle_internal"][triangle, row_position])
                if row < 0:
                    continue
                gradient = area_cache[tri_pattern]["gradient"][row_position]
                if threshold is not None and abs(gradient) <= threshold:
                    continue
                for column_position in range(10):
                    column = int(
                        schedule["simplex_internal"][simplex, column_position]
                    )
                    if column < 0:
                        continue
                    return {
                        "simplex": simplex,
                        "hinge": hinge,
                        "triangle": triangle,
                        "row_position": row_position,
                        "column_position": column_position,
                        "gradient_absolute": float(abs(gradient)),
                        "matrix_entry_change": float(abs(gradient) * 1e-4),
                    }
    raise RuntimeError("no admissible corruption incidence found")


def componentwise_metrics(matrix, solution, rhs):
    residual = rhs - matrix @ solution
    raw = np.linalg.norm(residual, ord=np.inf) / max(
        1.0, np.linalg.norm(rhs, ord=np.inf)
    )
    denominator = abs(matrix) @ np.abs(solution) + np.abs(rhs)
    ratios = np.divide(
        np.abs(residual), denominator,
        out=np.zeros_like(denominator), where=denominator > 0,
    )
    return residual, float(raw), float(np.max(ratios))


def spectral_diagnostic(source, matrix, tangent, operator_error, frozen_pair):
    size = matrix.shape[0]
    column = sp.csr_matrix(tangent.reshape(-1, 1))
    bordered = sp.bmat(
        [[matrix, column], [column.T, sp.csr_matrix((1, 1))]], format="csc"
    )
    factor = spla.splu(bordered, permc_spec="COLAMD", diag_pivot_thresh=1.0)
    inverse = spla.LinearOperator(
        bordered.shape,
        matvec=factor.solve,
        rmatvec=lambda value: factor.solve(value, trans="T"),
        dtype=np.float64,
    )

    max_row_nnz = int(np.diff(bordered.tocsr().indptr).max())
    unit_roundoff = np.finfo(np.float64).eps / 2
    gamma_q = max_row_nnz * unit_roundoff / (
        1 - max_row_nnz * unit_roundoff
    )
    backward_gate = 100 * gamma_q
    solves = []
    all_backward_stable = True
    for probe_index in range(8):
        rhs = source["fixed_probe"](size + 1, probe_index)
        solution = factor.solve(rhs)
        ladder = []
        for refinement in range(4):
            residual, raw, eta = componentwise_metrics(
                bordered, solution, rhs
            )
            ladder.append({
                "refinement": refinement,
                "raw_relative_residual": ftext(raw),
                "componentwise_backward_error": ftext(eta),
            })
            if refinement < 3:
                solution += factor.solve(residual)
        all_backward_stable &= float(ladder[0]["componentwise_backward_error"]) <= backward_gate
        solves.append({"probe_index": probe_index, "ladder": ladder})

    eigen_runs = []
    for run_index, tolerance in enumerate((1e-10, 1e-12)):
        values, vectors = spla.eigsh(
            bordered,
            k=32,
            sigma=0.0,
            which="LM",
            OPinv=inverse,
            v0=source["fixed_probe"](size + 1, run_index + 1),
            tol=tolerance,
            maxiter=30000,
        )
        order = np.argsort(np.abs(values))
        values = values[order]
        vectors = vectors[:, order]
        residuals = np.asarray([
            np.linalg.norm(bordered @ vectors[:, index]
                           - values[index] * vectors[:, index])
            for index in range(32)
        ])
        separated = np.asarray([
            abs(values[index])
            > 100 * (operator_error + residuals[index])
            for index in range(32)
        ])
        eigen_runs.append({
            "tolerance": ftext(tolerance),
            "eigenvalues_nearest_zero": [ftext(value) for value in values],
            "ritz_residuals": [ftext(value) for value in residuals],
            "separated_from_zero": [bool(value) for value in separated],
            "all_observed_values_separated": bool(np.all(separated)),
        })

    first = np.asarray([
        float(value) for value in eigen_runs[0]["eigenvalues_nearest_zero"]
    ])
    second = np.asarray([
        float(value) for value in eigen_runs[1]["eigenvalues_nearest_zero"]
    ])
    max_ritz = max(
        max(float(value) for value in eigen_runs[0]["ritz_residuals"]),
        max(float(value) for value in eigen_runs[1]["ritz_residuals"]),
    )
    run_gate = 100 * (operator_error + max_ritz)
    paired = bool(np.max(np.abs(first - second)) <= run_gate)

    old_values = np.asarray([
        float(value) for value in
        frozen_pair["bordered_spectrum"]["eigen_runs"][0][
            "eigenvalues_nearest_zero"
        ]
    ])
    old_ritz = float(
        frozen_pair["bordered_spectrum"]["maximum_ritz_residual"]
    )
    reproduction_gate = 100 * (operator_error + max_ritz + old_ritz)
    old_reproduced = bool(
        np.max(np.abs(first[:8] - old_values)) <= reproduction_gate
    )
    return {
        "matrix_infinity_norm": ftext(source["sparse_row_norm"](bordered)),
        "maximum_row_nnz": max_row_nnz,
        "backward_stability_gate": ftext(backward_gate),
        "solves": solves,
        "all_initial_solves_backward_stable": bool(all_backward_stable),
        "eigen_runs": eigen_runs,
        "paired_run_gate": ftext(run_gate),
        "paired_eigensolves": paired,
        "old_eigenvalue_reproduction_gate": ftext(reproduction_gate),
        "old_eigenvalues_reproduced": old_reproduced,
        "observed_32_exhaustive_for_kernel": False,
        "factor_permutation_row_sha256": sha256(
            np.asarray(factor.perm_r, dtype="<i4").tobytes()
        ).hexdigest(),
        "factor_permutation_column_sha256": sha256(
            np.asarray(factor.perm_c, dtype="<i4").tobytes()
        ).hexdigest(),
    }


print("=" * 78)
print("FOCUSED NONHOMOGENEOUS INTERNAL-HESSIAN DIAGNOSTIC")
print("=" * 78)

actual_hashes = {name: digest(ROOT / name) for name in EXPECTED_HASHES}
provenance_ok = check(
    "the frozen failed run and diagnostic protocol have exact provenance",
    actual_hashes == EXPECTED_HASHES and PROTOCOL_COMMIT == "cdd0f69",
    str(actual_hashes),
)

frozen = json.loads(FROZEN.read_text())
curvature = json.loads(CURVATURE.read_text())
source = load_source_definitions()
actions = source["load_action_definitions"]()
definitions_ok = check(
    "only definitions are loaded from both frozen verifiers",
    "OUTPUT" not in actions
    and "OUTPUT" not in source
    and {"schedule_geometry", "assemble_internal", "aggregate_hessian"}
        <= set(source),
)

_, adjacency, _ = build_600cell()
coarse_top = actions["tetrahedra_from_adjacency"](adjacency)
_, top, colours = actions["barycentric_chambers"](coarse_top)
orders = tuple(permutations(range(4)))

print("Building the 24 combinatorial schedules (no full matrices)...", flush=True)
all_schedules = tuple(
    source["schedule_geometry"](actions, top, colours, order)
    for order in orders
)
_, all_triangle_kinds = source["global_pattern_catalogue"](all_schedules)
topology_ok = check(
    "the declared carrier, lexicographic schedule and 12 representatives are fixed",
    len(colours) == EXPECTED_F[0]
    and len(all_schedules) == 24
    and all(schedule["counts"] == EXPECTED_COUNTS for schedule in all_schedules)
    and all_schedules[0]["order"] == (0, 1, 2, 3)
    and tuple(
        pair["forward_schedule"]
        for pair in frozen["census"]["time_reversal_pairs"]
    ) == FORWARD_INDICES,
)

geometries = {dps: actions["exact_geometry"](dps) for dps in (100, 140)}
for geometry in geometries.values():
    geometry["mass"] = mp.mpf(0)
masses_total = tuple(
    mp.mpf(value)
    for value in curvature["selected_rank_matter"]["total_masses"]
)
masses_per_vertex = tuple(
    mp.mpf(value)
    for value in curvature["selected_rank_matter"]["per_vertex_masses"]
)

# Rebuild schedule 0 with a local pattern catalogue, avoiding the 96-pattern
# all-schedule stencil census.  Pattern values depend only on their kind.
schedule0 = source["schedule_geometry"](
    actions, top, colours, orders[0]
)
simplex_kinds0, triangle_kinds0 = source["global_pattern_catalogue"]((schedule0,))
print("Building only schedule-0 high-precision stencils...", flush=True)
angle_binary, angle_extended = angle_caches(
    source, actions, simplex_kinds0, geometries
)
area_binary0, area_extended0 = area_caches(
    source, triangle_kinds0, geometries
)

matrix0, matrix_diagnostics = source["assemble_internal"](
    schedule0, angle_binary, area_binary0, masses_per_vertex, geometries[140]
)
matrix_digest = source["csr_digest"](matrix0)
expected_digest = frozen["census"]["schedules"][0]["csr_sha256"]
digest_ok = check(
    "schedule 0 reproduces the exact frozen CSR matrix",
    matrix_digest == expected_digest,
    f"actual={matrix_digest}, expected={expected_digest}",
)

p64 = source["orbit_pullback"](matrix0, schedule0["internal_orbits"])
pld, pld_error, pld_counts, pld_diagnostics = direct_longdouble_pullback(
    schedule0, angle_extended, area_extended0, masses_per_vertex
)
aggregate_combinatorics = actions["schedule_combinatorics"](
    top, colours, orders[0]
)
pagg, pagg_error = source["aggregate_hessian"](
    actions, aggregate_combinatorics, geometries, masses_total
)
p64_error = np.longdouble(matrix_diagnostics["operator_error_row_bound"])

comparisons = {
    "P64_vs_PLD": matrix_comparison(p64, pld, p64_error, pld_error),
    "P64_vs_PAGG": matrix_comparison(p64, pagg, p64_error, pagg_error),
    "PLD_vs_PAGG": matrix_comparison(pld, pagg, pld_error, pagg_error),
}
frozen_difference = float(
    frozen["census"]["schedules"][0][
        "aggregate_pullback_maximum_difference"
    ]
)
frozen_pullback_reproduced = check(
    "the known P64-PAGG discrepancy is independently reproduced",
    abs(float(comparisons["P64_vs_PAGG"]["maximum_absolute_difference"])
        - frozen_difference)
    <= 100 * (matrix_diagnostics["operator_error_row_bound"] + pagg_error),
    str(comparisons["P64_vs_PAGG"]),
)

match_pattern = tuple(
    comparisons[name]["matches"]
    for name in ("P64_vs_PLD", "P64_vs_PAGG", "PLD_vs_PAGG")
)
if match_pattern == (False, False, True):
    pullback_outcome = "BINARY64_ACCUMULATION_OR_PULLBACK_ENVELOPE_LOCALIZED"
elif match_pattern == (True, False, False):
    pullback_outcome = "AGGREGATE_FINITE_DIFFERENCE_CONTROL_LOCALIZED"
elif match_pattern == (False, True, False):
    pullback_outcome = "DIRECT_LONGDOUBLE_SUM_INVALID"
elif match_pattern == (True, True, True):
    pullback_outcome = "ALL_THREE_COMPATIBLE_OLD_ENVELOPE_TOO_TIGHT"
else:
    pullback_outcome = "PULLBACK_FORMULA_OPEN"
pullback_hierarchy_ok = check(
    "the three-way pullback result is classified by the frozen hierarchy",
    pullback_outcome in {
        "BINARY64_ACCUMULATION_OR_PULLBACK_ENVELOPE_LOCALIZED",
        "AGGREGATE_FINITE_DIFFERENCE_CONTROL_LOCALIZED",
        "DIRECT_LONGDOUBLE_SUM_INVALID",
        "ALL_THREE_COMPATIBLE_OLD_ENVELOPE_TOO_TIGHT",
        "PULLBACK_FORMULA_OPEN",
    },
    f"pattern={match_pattern}, outcome={pullback_outcome}",
)

print("Building exact area gradients for corruption controls...", flush=True)
all_area_binary = source["build_area_cache"](
    all_triangle_kinds, geometries[140]
)
corruption_records = []
old_zero = []
corrected_all_detected = True
for forward_index in FORWARD_INDICES:
    schedule = all_schedules[forward_index]
    old = select_corruption(schedule, all_area_binary)
    corrected = select_corruption(
        schedule, all_area_binary, float(NONZERO_GRADIENT_THRESHOLD)
    )
    operator_error = float(
        frozen["census"]["schedules"][forward_index]["diagnostics"][
            "operator_error_row_bound"
        ]
    )
    detected = corrected["matrix_entry_change"] > 100 * operator_error
    corrected_all_detected &= detected
    if old["matrix_entry_change"] == 0:
        old_zero.append(forward_index)
    corruption_records.append({
        "forward_schedule": forward_index,
        "old_selection": old,
        "corrected_selection": corrected,
        "detection_gate": ftext(100 * operator_error),
        "corrected_detected": bool(detected),
    })
old_control_ok = check(
    "the old corruption rule reproduces exactly its four known null choices",
    tuple(old_zero) == EXPECTED_OLD_ZERO,
    f"zero schedules={old_zero}",
)
corrected_control_ok = check(
    "the preregistered nonzero-gradient corruption is detected in all 12 pairs",
    corrected_all_detected,
)

tangent0 = source["product_tangent"](schedule0, geometries[140])
print("Factoring only bordered schedule 0 and requesting 32 Ritz values...", flush=True)
frozen_pair0 = frozen["census"]["time_reversal_pairs"][0]
try:
    spectral = spectral_diagnostic(
        source,
        matrix0,
        tangent0,
        matrix_diagnostics["operator_error_row_bound"],
        frozen_pair0,
    )
    factorization_ok = True
except Exception as error:
    spectral = {"error": repr(error)}
    factorization_ok = False

check("the single preregistered bordered factorization completes", factorization_ok)
backward_ok = check(
    "all eight initial LU solves pass the componentwise backward-error gate",
    factorization_ok and spectral["all_initial_solves_backward_stable"],
    spectral.get("backward_stability_gate", "factorization failed"),
)
eigen_pair_ok = check(
    "the two 32-vector Ritz runs agree within the operator envelope",
    factorization_ok and spectral["paired_eigensolves"],
)
old_eigen_ok = check(
    "the first eight frozen Ritz values are reproduced",
    factorization_ok and spectral["old_eigenvalues_reproduced"],
)
observed_separated = bool(
    factorization_ok
    and all(run["all_observed_values_separated"]
            for run in spectral["eigen_runs"])
)
check(
    "the diagnostic does not promote 32 observed Ritz values to an exhaustive kernel proof",
    factorization_ok
    and observed_separated
    and not spectral["observed_32_exhaustive_for_kernel"],
    "observed cluster may be resolved; completeness remains OPEN",
)

controls_ok = all((
    provenance_ok,
    definitions_ok,
    topology_ok,
    digest_ok,
    frozen_pullback_reproduced,
    pullback_hierarchy_ok,
    old_control_ok,
    corrected_control_ok,
    factorization_ok,
    backward_ok,
    eigen_pair_ok,
    old_eigen_ok,
    observed_separated,
))
if not controls_ok:
    outcome = "DIAGNOSTIC_INVALID"
elif pullback_outcome == "PULLBACK_FORMULA_OPEN":
    outcome = "VALID_DIAGNOSTIC_PULLBACK_OPEN"
else:
    outcome = "VALID_DIAGNOSTIC_NO_KERNEL_VERDICT"

artifact = {
    "title": "Focused diagnostic of the first complete refined internal Hessian run",
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "scope": {
        "factorized_schedule_indices": [0],
        "full_census_rerun": False,
        "full_suite_run": False,
        "physical_target_loaded": False,
        "kernel_verdict_permitted": False,
    },
    "schedule0": {
        "order": list(schedule0["order"]),
        "csr_sha256": matrix_digest,
        "operator_error_row_bound": ftext(
            matrix_diagnostics["operator_error_row_bound"]
        ),
        "local_simplex_pattern_count": len(simplex_kinds0),
        "local_triangle_pattern_count": len(triangle_kinds0),
    },
    "pullback_diagnostic": {
        "P64": [[ftext(value) for value in row] for row in p64],
        "PLD": [[ldtext(value) for value in row] for row in pld],
        "PAGG": [[ftext(value) for value in row] for row in pagg],
        "PLD_error": [[ldtext(value) for value in row] for row in pld_error],
        "PLD_term_counts": pld_counts.tolist(),
        "PLD_diagnostics": pld_diagnostics,
        "PAGG_error": ftext(pagg_error),
        "comparisons": comparisons,
        "match_pattern_P64_PLD__P64_PAGG__PLD_PAGG": list(match_pattern),
        "outcome": pullback_outcome,
    },
    "corruption_diagnostic": {
        "threshold": mp.nstr(NONZERO_GRADIENT_THRESHOLD, 10),
        "expected_old_zero_schedules": list(EXPECTED_OLD_ZERO),
        "observed_old_zero_schedules": old_zero,
        "corrected_all_detected": bool(corrected_all_detected),
        "records": corruption_records,
    },
    "spectral_diagnostic": spectral,
    "status_labels": {
        "pullback": pullback_outcome,
        "old_raw_solve_residual_as_eigenvalue_uncertainty": (
            "CATEGORY_ERROR" if backward_ok else "OPEN"
        ),
        "observed_soft_cluster": (
            "DERIVED_COMPUTATIONAL_FOR_32_OBSERVED_RITZ_PAIRS"
            if observed_separated else "OPEN"
        ),
        "complete_internal_kernel": "OPEN",
        "physics": "NOT_TESTED",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"Pullback outcome: {pullback_outcome}")
print(f"Observed 32 Ritz values separated: {observed_separated}")
print("Complete kernel: OPEN (exhaustiveness not certified)")
print(f"Outcome: {outcome}")
print(f"Tests: {passed}/{tests} passed")
print(f"Artifact: {OUTPUT.relative_to(ROOT)}")
sys.exit(0 if passed == tests else 1)
