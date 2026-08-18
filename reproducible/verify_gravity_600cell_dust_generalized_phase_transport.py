#!/usr/bin/env python3
"""Canonical phase transport audit for the generalized mode fibers.

Prior-art/framing commit: 60aafe1.
Preregistered protocol commit: 3794418.
"""

from collections import Counter
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy

import mpmath as mp
from flint import ctx


HERE = Path(__file__).resolve().parent
RESIDUAL_SOURCE = (
    HERE / "verify_gravity_600cell_dust_generalized_bundle_residual.py"
)
RESIDUAL_ARTIFACT = (
    HERE / "gravity_600cell_dust_generalized_bundle_residual.json"
)
OUTPUT = HERE / "gravity_600cell_dust_generalized_phase_transport.json"

PRIOR_ART_COMMIT = "60aafe1"
PROTOCOL_COMMIT = "3794418"
EXPECTED_HASHES = {
    "residual_source": (
        "ccf2ebe03c6e39c3d6e6b538d1c02d278804553987d65db0eeb67fce7936ca5a"
    ),
    "residual_artifact": (
        "3244185127aecf7c9a44261cced0be521c9dc42bf8e44f909d8a0ce10a96eadf"
    ),
}

DPS = 100
BALL_DPS = 80
MP_FLOOR = mp.mpf("1e-70")
PARITIES = ("even", "odd")
TARGET_SECTORS = (4, 5)
VARIANTS = (
    "operational_primary",
    "operational_shadow",
    "validation_primary",
    "validation_shadow",
)

mp.mp.dps = DPS
ctx.dps = BALL_DPS
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


def smp(value, digits=30):
    value = mp.mpf(value)
    if mp.isinf(value):
        return "inf" if value > 0 else "-inf"
    if mp.isnan(value):
        return "nan"
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def block_diagonal_phase(projector):
    n = projector.rows
    result = mp.matrix(2 * n, 2 * n)
    conjugate = projector.conjugate()
    for row in range(n):
        for column in range(n):
            result[row, column] = projector[row, column]
            result[n + row, n + column] = conjugate[row, column]
    return result


def assemble_blocks(blocks):
    n = blocks["A"].rows
    result = mp.matrix(2 * n, 2 * n)
    for row in range(n):
        for column in range(n):
            result[row, column] = blocks["A"][row, column]
            result[row, n + column] = blocks["B"][row, column]
            result[n + row, column] = blocks["C"][row, column]
            result[n + row, n + column] = blocks["D"][row, column]
    return result


def leakage_label(value, error):
    if not mp.isfinite(value) or not mp.isfinite(error):
        return "LEAKAGE_OPEN"
    if value <= 10 * error:
        return "LEAKAGE_ZERO_CONSISTENT"
    if value > 100 * error:
        return "LEAKAGE_NONZERO_RESOLVED"
    return "LEAKAGE_OPEN"


hashes = {
    "residual_source": sha256(RESIDUAL_SOURCE),
    "residual_artifact": sha256(RESIDUAL_ARTIFACT),
}
provenance_ok = hashes == EXPECTED_HASHES
check("the residual-certified bundle inputs retain exact provenance",
      provenance_ok, str(hashes))

print("[setup] replaying residual-certified high-precision projectors",
      flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    residual = runpy.run_path(str(RESIDUAL_SOURCE))
residual_replay_ok = bool(
    residual["passed"] == residual["tests"] == 10
    and residual["outcome"]
    == "RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED"
    and sha256(RESIDUAL_ARTIFACT) == EXPECTED_HASHES["residual_artifact"]
    and len(residual["projectors"]) == 32
)
check("the residual-certified projectors replay byte-identically",
      residual_replay_ok)

prior = residual["prior"]
mp_operator_norm = residual["mp_operator_norm"]
acb_midpoint_to_mp = residual["acb_midpoint_to_mp"]
acb_radius_frobenius = residual["acb_radius_frobenius"]

print("[setup] reconstructing the exact second-slab tangent balls", flush=True)
tangent_cells = {}
reconstruction_records = []
reconstruction_ok = True
for parity in PARITIES:
    a1, _ = [
        mp.mpf(value)
        for value in prior["first_tick"]["solutions"][parity]["state"]
    ]
    a2, r2 = [
        mp.mpf(value)
        for value in prior["second_tick"]["solutions"][parity][
            "state_absolute"
        ]
    ]
    index_data, kind_values = prior["slab_index_data"](
        prior["models"][parity], a1, a2, r2
    )
    geometry = prior["prepare_geometry"](
        prior["models"][parity], index_data
    )
    mapping = prior["boundary_mapping"](index_data)
    sectors, sector_control = prior["high_precision_sector_bases"](index_data)
    pattern_cache, branch = prior["high_precision_pattern_cache"](
        geometry["patterns"], kind_values
    )
    kernels, kernel_control = prior["assemble_full_representative_kernels"](
        index_data, geometry, pattern_cache
    )
    parity_control = bool(
        tuple(mapping) == tuple(range(30))
        and branch["entry_pass"]
        and kernel_control["maximum_imaginary"] < prior["ARITHMETIC_FLOOR"]
        and tuple(item["dimension"] for item in sectors)
        == (3, 2, 2, 2, 1, 1, 1)
        and sector_control["isotypic_dimensions"] == [1, 1, 1, 4, 4, 4, 9]
    )
    reconstruction_ok &= parity_control

    for sector_index in TARGET_SECTORS:
        sector = sectors[sector_index]
        reconstruction_ok &= sector["dimension"] == 1
        for variant in VARIANTS:
            raw = prior["project_full_kernel"](kernels[variant], sector)
            block = (raw + raw.H) / 2
            _, determinant, tangent_ball, symplectic_defect = (
                prior["build_tangent_ball"](block, 1, mapping)
            )
            determinant_ok = not determinant.contains(0)
            symplectic_ok = all(
                symplectic_defect[row, column].contains(0)
                for row in range(symplectic_defect.nrows())
                for column in range(symplectic_defect.ncols())
            )
            principal, identities = prior["reconstruct_principal"](
                prior["split_tangent"](tangent_ball, 30)
            )
            principal_ok = all(
                identity[row, column].contains(0)
                for identity in identities.values()
                for row in range(identity.nrows())
                for column in range(identity.ncols())
            )
            complete = bool(
                parity_control and determinant_ok and symplectic_ok
                and principal_ok
            )
            reconstruction_ok &= complete
            key = (parity, sector_index, variant)
            tangent_cells[key] = tangent_ball
            reconstruction_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "identity_boundary_mapping": tuple(mapping) == tuple(range(30)),
                "internal_determinant_excludes_zero": determinant_ok,
                "symplectic_ball_contains_zero": symplectic_ok,
                "principal_identities_contain_zero": principal_ok,
                "complete": complete,
            })

check("all 16 exact second-slab tangent balls reconstruct canonically",
      reconstruction_ok
      and len(tangent_cells) == len(reconstruction_records) == 16)

block_records = []
full_records = []
block_counts = Counter()
full_counts = Counter()
assembly_ok = True
compatibility_ok = True
all_finite = True
for parity in PARITIES:
    for sector_index in TARGET_SECTORS:
        for variant in VARIANTS:
            key = (parity, sector_index, variant)
            old = residual["projectors"][
                ("old", parity, sector_index, variant)
            ]
            shifted = residual["projectors"][
                ("shifted", parity, sector_index, variant)
            ]
            p0 = old["projector"]
            p1 = shifted["projector"]
            cp0 = p0.conjugate()
            cp1 = p1.conjugate()
            one = mp.eye(30)
            tangent_ball = tangent_cells[key]
            tangent = acb_midpoint_to_mp(tangent_ball)
            rows0 = tuple(range(30))
            rows1 = tuple(range(30, 60))
            block_balls = {
                "A": prior["submatrix"](tangent_ball, rows0, rows0),
                "B": prior["submatrix"](tangent_ball, rows0, rows1),
                "C": prior["submatrix"](tangent_ball, rows1, rows0),
                "D": prior["submatrix"](tangent_ball, rows1, rows1),
            }
            blocks = {
                name: acb_midpoint_to_mp(ball)
                for name, ball in block_balls.items()
            }
            residual_blocks = {
                "A": (one - p1) * blocks["A"] * p0,
                "B": (one - p1) * blocks["B"] * cp0,
                "C": (one - cp1) * blocks["C"] * p0,
                "D": (one - cp1) * blocks["D"] * cp0,
            }
            cell_block_labels = []
            for name in ("A", "B", "C", "D"):
                matrix = blocks[name]
                matrix_norm = mp_operator_norm(matrix)
                epsilon_x = (
                    acb_radius_frobenius(block_balls[name])
                    + MP_FLOOR * max(mp.mpf(1), matrix_norm)
                )
                residual_norm = mp_operator_norm(residual_blocks[name])
                residual_error = (
                    epsilon_x
                    + (old["eta"] + shifted["eta"]
                       + old["eta"] * shifted["eta"])
                    * (matrix_norm + epsilon_x)
                    + MP_FLOOR * max(mp.mpf(1), matrix_norm)
                )
                label = leakage_label(residual_norm, residual_error)
                cell_block_labels.append(label)
                block_counts[(name, label)] += 1
                all_finite &= bool(
                    mp.isfinite(matrix_norm)
                    and mp.isfinite(residual_norm)
                    and mp.isfinite(residual_error)
                )
                block_records.append({
                    "parity": parity,
                    "sector_index": sector_index,
                    "variant": variant,
                    "block": name,
                    "operator_norm": smp(matrix_norm),
                    "operator_error": smp(epsilon_x),
                    "residual_norm": smp(residual_norm),
                    "residual_error": smp(residual_error),
                    "error_units": smp(residual_norm / residual_error),
                    "label": label,
                })

            q0 = block_diagonal_phase(p0)
            q1 = block_diagonal_phase(p1)
            full_residual = (mp.eye(60) - q1) * tangent * q0
            assembled = assemble_blocks(residual_blocks)
            assembly_error = mp_operator_norm(full_residual - assembled)
            assembly_scale = max(mp.mpf(1), mp_operator_norm(full_residual))
            cell_assembly_ok = assembly_error / assembly_scale < mp.mpf("1e-60")
            assembly_ok &= cell_assembly_ok
            tangent_norm = mp_operator_norm(tangent)
            epsilon_t = (
                acb_radius_frobenius(tangent_ball)
                + MP_FLOOR * max(mp.mpf(1), tangent_norm)
            )
            residual_norm = mp_operator_norm(full_residual)
            residual_error = (
                epsilon_t
                + (old["eta"] + shifted["eta"]
                   + old["eta"] * shifted["eta"])
                * (tangent_norm + epsilon_t)
                + MP_FLOOR * max(mp.mpf(1), tangent_norm)
            )
            full_label = leakage_label(residual_norm, residual_error)
            full_counts[full_label] += 1
            if "LEAKAGE_NONZERO_RESOLVED" in cell_block_labels:
                cell_compatible = full_label == "LEAKAGE_NONZERO_RESOLVED"
            elif all(
                label == "LEAKAGE_ZERO_CONSISTENT"
                for label in cell_block_labels
            ):
                cell_compatible = full_label == "LEAKAGE_ZERO_CONSISTENT"
            else:
                cell_compatible = full_label == "LEAKAGE_OPEN"
            compatibility_ok &= cell_compatible
            all_finite &= bool(
                mp.isfinite(tangent_norm)
                and mp.isfinite(residual_norm)
                and mp.isfinite(residual_error)
            )
            full_records.append({
                "parity": parity,
                "sector_index": sector_index,
                "variant": variant,
                "tangent_norm": smp(tangent_norm),
                "tangent_error": smp(epsilon_t),
                "residual_norm": smp(residual_norm),
                "residual_error": smp(residual_error),
                "error_units": smp(residual_norm / residual_error),
                "assembly_relative_error": smp(assembly_error / assembly_scale),
                "assembly_ok": bool(cell_assembly_ok),
                "block_labels": cell_block_labels,
                "label_compatible": bool(cell_compatible),
                "label": full_label,
            })

block_complete = len(block_records) == sum(block_counts.values()) == 64
full_complete = len(full_records) == sum(full_counts.values()) == 16
check("all 64 block and 16 full phase residuals receive labels",
      block_complete and full_complete and all_finite,
      f"blocks={dict(block_counts)}, full={dict(full_counts)}")
check("the direct full residuals equal their four-block assemblies",
      assembly_ok)
check("all full and block leakage labels are compatible",
      compatibility_ok)

flat_block_counts = Counter(item["label"] for item in block_records)
controls_ok = bool(
    provenance_ok and residual_replay_ok and reconstruction_ok
    and len(tangent_cells) == 16 and block_complete and full_complete
    and all_finite and assembly_ok and compatibility_ok
)
if not controls_ok:
    outcome = "GENERALIZED_PHASE_TRANSPORT_CONTROL_FAILED"
elif (
    flat_block_counts["LEAKAGE_NONZERO_RESOLVED"]
    or full_counts["LEAKAGE_NONZERO_RESOLVED"]
):
    outcome = "GENERALIZED_PHASE_TRANSPORT_REFUTED"
elif (
    flat_block_counts["LEAKAGE_OPEN"]
    or full_counts["LEAKAGE_OPEN"]
):
    outcome = "GENERALIZED_PHASE_TRANSPORT_OPEN"
else:
    outcome = "GENERALIZED_PHASE_TRANSPORT_CERTIFIED"

allowed = {
    "GENERALIZED_PHASE_TRANSPORT_CONTROL_FAILED",
    "GENERALIZED_PHASE_TRANSPORT_REFUTED",
    "GENERALIZED_PHASE_TRANSPORT_OPEN",
    "GENERALIZED_PHASE_TRANSPORT_CERTIFIED",
}
check("the preregistered phase-transport hierarchy assigns one outcome",
      outcome in allowed, outcome)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": hashes,
    "outcome": outcome,
    "controls_ok": controls_ok,
    "fitted_alignment_used": False,
    "reconstruction_records": reconstruction_records,
    "block_counts": {
        f"{block}:{label}": count
        for (block, label), count in sorted(block_counts.items())
    },
    "block_records": block_records,
    "full_counts": dict(full_counts),
    "full_records": full_records,
    "classification": {
        "generalized_configuration_rotation": "DERIVED COMPUTATIONAL",
        "full_cotangent_phase_transport": (
            "DERIVED COMPUTATIONAL"
            if outcome == "GENERALIZED_PHASE_TRANSPORT_CERTIFIED"
            else "DERIVED COMPUTATIONAL REFUTATION"
            if outcome == "GENERALIZED_PHASE_TRANSPORT_REFUTED"
            else "OPEN"
        ),
        "transported_intersection_or_lagrangian_graph": "NOT COMPUTED",
        "reduced_propagator": "NOT COMPUTED",
        "dispersion_mass_inertia_or_speed": "NOT COMPUTED",
    },
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

print("-" * 78)
print(f"SCIENTIFIC OUTCOME: {outcome}")
print(f"block labels: {dict(block_counts)}")
print(f"full labels: {dict(full_counts)}")
print(f"Tests: {passed}/{tests}")
print(f"Artifact: {OUTPUT}")
if passed != tests:
    raise SystemExit(1)
