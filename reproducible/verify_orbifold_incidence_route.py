#!/usr/bin/env python3
"""Corrected comparison for the canonical binary (2,3,5) incidence route."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import runpy

import sympy as sy


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "orbifold_incidence_preregistered.json").read_text(
    encoding="utf-8"
))
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


def quiet_run(path):
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        return runpy.run_path(str(path))


print("=" * 78)
print("CORRECTED (2,3,5) ORBIFOLD INCIDENCE ROUTE")
print("=" * 78)

# ---------------------------------------------------------------------------
# 1. Derive the comparison modules independently from the exact 2I table.
# ---------------------------------------------------------------------------

harmonic = quiet_run(HERE / "verify_fibonacci_nonbinary_dynamics.py")
irreps = harmonic["irreps_2i"]
dims = tuple(harmonic["irrep_dims"])
class_sizes = (1, 1, 30, 20, 20, 12, 12, 12, 12)

# Repository order is rho1,...,rho9.  The derived diagonal homomorphism uses
# color rho5=3' and weak rho2=2.  Compute their tensor product rather than
# inserting the McKay edge as an assumption.
color_weak_product = tuple(
    sy.simplify(sum(
        class_sizes[c] * irreps[4][c] * irreps[1][c] * irreps[i][c]
        for c in range(9)
    ) / 120)
    for i in range(9)
)
check("[DERIVED] rho5 tensor rho2 is exactly rho9",
      color_weak_product == (0, 0, 0, 0, 0, 0, 0, 0, 1))

m16 = tuple(
    int(color_weak_product[i])
    + (2 if i == 4 else 0)
    + (1 if i == 1 else 0)
    + (2 if i == 0 else 0)
    for i in range(9)
)
m15 = tuple(value-(1 if i == 0 else 0) for i, value in enumerate(m16))
three_m16 = tuple(3*value for value in m16)
three_m15 = tuple(3*value for value in m15)
targets = {
    "M15": m15,
    "M16": m16,
    "3*M15": three_m15,
    "3*M16": three_m16,
}

central_parity = tuple(int(character[1]/character[0]) for character in irreps)
check("[DERIVED] diagonal M16 has the audited character, dimension, and center trace",
      m16 == (2, 1, 0, 0, 2, 0, 0, 0, 1)
      and sum(a*b for a, b in zip(m16, dims)) == 16
      and sum(a*character[1] for a, character in zip(m16, irreps)) == 0
      and central_parity == (1, -1, -1, 1, 1, 1, -1, 1, -1))
check("[DERIVED] M15 is M16 with one trivial sterile singlet removed",
      m15 == (1, 1, 0, 0, 2, 0, 0, 0, 1)
      and sum(a*b for a, b in zip(m15, dims)) == 15)

# ---------------------------------------------------------------------------
# 2. Quantify why a virtual induction identity carries no selection evidence.
# ---------------------------------------------------------------------------

induction = sy.Matrix([
    record["irrep_multiplicities"] for record in DATA["modules"]
])
particular = sy.Matrix([
    1, -1, 0, 0, -2, 0, 0, 0, 0, 0,  # C10/V; generator is 10B-labelled
    0, 0, 0, 0,                         # C4/E harmonics
    1, 1, 0, 0, 0, 0,                   # C6/F harmonics
])
check("[DERIVED] the displayed virtual identity is one integer solution",
      induction.T * particular == sy.Matrix(m16))
check("[DERIVED] all virtual solutions form an affine lattice of rank eleven",
      induction.rank() == 9 and len(induction.T.nullspace()) == 11,
      "20 induction generators minus rank 9 = 11 free integral directions")

# ---------------------------------------------------------------------------
# 3. Compare only operator-dependent outputs: ker, coker, and middle H1.
# ---------------------------------------------------------------------------

operators = DATA["operators"]
complexes = DATA["short_complexes"]

operator_outputs = []
object_hits = {name: [] for name in targets}
slot_hits = {name: [] for name in targets}
for record in operators:
    outputs = {
        "kernel": tuple(record["kernel_irrep_multiplicities"]),
        "cokernel": tuple(record["cokernel_irrep_multiplicities"]),
    }
    operator_outputs.extend(outputs.values())
    for target_name, target in targets.items():
        matching_slots = [name for name, character in outputs.items()
                          if character == target]
        if matching_slots:
            object_hits[target_name].append(record["id"])
            slot_hits[target_name].extend(
                f"{record['id']}:{slot}" for slot in matching_slots
            )

middle_outputs = []
for record in complexes:
    middle = tuple(record["H1_irrep_multiplicities"])
    middle_outputs.append(middle)
    for target_name, target in targets.items():
        if middle == target:
            object_hits[target_name].append(record["id"])
            slot_hits[target_name].append(f"{record['id']}:H1")

all_outputs = operator_outputs + middle_outputs
output_dimensions = [sum(a*b for a, b in zip(character, dims))
                     for character in all_outputs]
check("[DERIVED NEGATIVE] no kernel, cokernel, or middle cohomology is a target",
      all(not hits for hits in object_hits.values()),
      f"object hits={object_hits}")
check("[DERIVED NEGATIVE] hit fraction is 0/63 objects and 0/125 output slots",
      len(operators) + len(complexes) == 63
      and len(operator_outputs) + len(middle_outputs) == 125
      and all(not hits for hits in slot_hits.values()))
check("[DERIVED NEGATIVE] output dimensions already exclude every target",
      not set(output_dimensions) & {15, 16, 45, 48},
      f"output dimensions={sorted(set(output_dimensions))}")

def nonzero_parities(character):
    return {central_parity[i] for i, multiplicity in enumerate(character)
            if multiplicity}


check("[DERIVED] every actual output has pure central parity",
      all(len(nonzero_parities(character)) <= 1 for character in all_outputs))
check("[DERIVED NEGATIVE] every comparison module has mixed central parity",
      all(nonzero_parities(character) == {-1, 1}
          for character in targets.values()),
      "a submodule of one induced line module cannot have the required mix")

# The sole genuine complex is the ordinary oriented cellular complex.  Its
# Euler character is endpoint data; H1 is the operator-dependent object.
check("[DERIVED] the sole complex is F0 -> E2 -> V0",
      len(complexes) == 1
      and complexes[0]["modules"] == ["F0", "E2", "V0"]
      and complexes[0]["composition_zero_exactly"])
check("[DERIVED NEGATIVE] its middle cohomology vanishes",
      complexes[0]["H1_dimension"] == 0
      and complexes[0]["H1_irrep_multiplicities"] == [0]*9)

# Virtual indices are included only to close the requested look-elsewhere
# ledger.  Their non-hit is not used as evidence about an operator.
virtual_indices = [tuple(record["index_irrep_multiplicities"])
                   for record in operators]
virtual_indices += [tuple(record["euler_index_irrep_multiplicities"])
                    for record in complexes]
virtual_hits = {
    name: sum(character == target for character in virtual_indices)
    for name, target in targets.items()
}
check("[DERIVED inventory] none of 63 endpoint indices equals a comparison module",
      virtual_hits == {name: 0 for name in targets},
      f"virtual hits={virtual_hits}; this is endpoint arithmetic only")

# The three untwisted ranks in the old run were real, although its 4/4' label
# order differs from the repository order.  The missing twists and complex,
# not these ranks, invalidate that run.
untwisted = {(record["source"], record["target"]): record
             for record in operators
             if record["source"] in {"E0", "F0"}
             and record["target"] in {"V0", "E0"}}
check("[DERIVED] old untwisted ranks survive the corrected calculation",
      untwisted[("E0", "V0")]["matrix_rank"] == 12
      and untwisted[("F0", "V0")]["matrix_rank"] == 12
      and untwisted[("F0", "E0")]["matrix_rank"] == 20)

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("TARGET_HITS_BY_OBJECT=0/63")
print("TARGET_HITS_BY_OUTPUT_SLOT=0/125")
print("VERDICT_OLD_N_EQUALS_3=REFUTED")
print("VERDICT_OLD_NO_SHORT_COMPLEX=REFUTED")
print("VERDICT_CORRECTED_LINE_INCIDENCE_ROUTE=KILLED")
if passed != tests:
    raise SystemExit(1)
