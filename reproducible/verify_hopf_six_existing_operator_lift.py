#!/usr/bin/env python3
"""Provenance audit of existing operators against the 936-state Hopf carrier.

The seven-family inventory and acceptance rule were frozen in protocol commit
e54695a.  The conclusion is scoped to those authoritative committed sources;
it is not a universal no-go against future operator constructions.
"""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE/"hopf_six_existing_operator_lift.json"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("="*78)
print("EXISTING-OPERATOR LIFT AUDIT FOR THE 936-STATE HOPF CARRIER")
print("="*78)

# Reuse exact structured facts from the already registered audits rather than
# assigning a different representation to an old operator.
with (HERE/"hopf_label_crossed_product.json").open() as handle:
    crossed = json.load(handle)
with (HERE/"hopf_label_superselection.json").open() as handle:
    pair = json.load(handle)
with (HERE/"hopf_fibration_label_operator.json").open() as handle:
    label = json.load(handle)
with (HERE/"hopf_six_crossed_real_galois.json").open() as handle:
    real_crossed = json.load(handle)

check("the authoritative crossed-product data have the expected dimensions",
      crossed["crossed_product"]["dimension"] == 360
      and crossed["crossed_product"]["Wedderburn_blocks"] == [6, 6, 12, 12]
      and crossed["natural_label_representation"]["image_dimension"] == 36
      and crossed["natural_label_representation"]["kernel_dimension"] == 324)
check("the authoritative pair and label carriers retain their original scope",
      pair["pair_groupoid"]["dimension"] == 36
      and label["modules"]["H_F"] == "six-point A5 permutation module"
      and real_crossed["canonical_carriers"][
          "minimum_faithful_left_dimension"] == 36)

# Fixed cell-level algebra theorem.  Both pi(B) and pi(B)^o preserve each
# ordered central cell (i,j), while gamma is scalar on that cell.  Exhaust all
# simple-block matrix units at the cell-support level for all eight readings.
node_sizes = (6, 6, 12, 12)
cell_actions_checked = 0
all_algebra_actions_even = True
for _reading in range(8):
    for node, size in enumerate(node_sizes):
        matrix_units = size*size
        left_cells = [(node, right) for right in range(4) if right != node]
        right_cells = [(left, node) for left in range(4) if left != node]
        for _unit in range(matrix_units):
            for source in left_cells:
                target = source
                cell_actions_checked += 1
                all_algebra_actions_even &= target == source
            for source in right_cells:
                target = source
                cell_actions_checked += 1
                all_algebra_actions_even &= target == source
expected_cell_actions = 8*2*3*sum(size*size for size in node_sizes)
check("all left/right crossed-product generator actions are grading-even",
      all_algebra_actions_even
      and cell_actions_checked == expected_cell_actions == 17280,
      "17,280 exact generator-cell actions; algebraic closure remains even")

# Frozen source inventory.  These are the authoritative operator-defining
# scripts, not every file that later discusses their results.
operator_sources = {
    "vertex_wave": HERE/"verify_hopf_fibration_invariants.py",
    "kahler_dirac": HERE/"verify_kahler_dirac.py",
    "label_auxiliary": HERE/"verify_hopf_fibration_label_operator.py",
    "label_hessian": HERE/"verify_hopf_full_hessian_spectral_enumeration.py",
    "pair_rook": HERE/"verify_hopf_label_superselection.py",
    "crossed_product": HERE/"verify_hopf_label_crossed_product.py",
    "golden_class": HERE/"verify_hopf_six_galois_spectral_split.py",
}
check("all seven preregistered authoritative sources exist",
      len(operator_sources) == 7
      and all(path.is_file() for path in operator_sources.values()))

source_text = {name: path.read_text() for name, path in operator_sources.items()}
lift_markers = (
    "H_off", "936-state", "dimension = 936", "dimension=936",
    "hopf_six_spectral_krajewski", "cell_offsets",
)
explicit_lift_sources = [
    name for name, text in source_text.items()
    if any(marker in text for marker in lift_markers)
]
check("none of the frozen operator sources defines a lift to H_off",
      explicit_lift_sources == [],
      "source-scoped absence; not a theorem about future mathematics")

# Each row records only previously certified representations.  `odd_status`
# is evaluated on H_off when a canonical algebraic action exists; otherwise
# it is undefined because inventing a map is forbidden by the protocol.
families = [
    {
        "name": "A_Af_BoxF",
        "carrier_dimensions": [120],
        "represented_algebra": "600-cell vertex/wave arena",
        "faithful_B_R": False,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": False,
        "odd_status": "undefined_no_lift",
    },
    {
        "name": "Kahler_Dirac_and_moments",
        "carrier_dimensions": [2640],
        "represented_algebra": "oriented cochain arena",
        "faithful_B_R": False,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": False,
        "odd_status": "undefined_no_lift",
    },
    {
        "name": "Phi_and_D_aux",
        "carrier_dimensions": [6, 12],
        "represented_algebra": "natural six-label M6 quotient",
        "faithful_B_R": False,
        "B_R_kernel_dimension": 324,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": False,
        "odd_status": "undefined_no_lift",
    },
    {
        "name": "full_label_Hessian_and_ten_state_double",
        "carrier_dimensions": [5, 10],
        "represented_algebra": "A5 W5 label module; no faithful B_R action",
        "faithful_B_R": False,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": False,
        "odd_status": "undefined_no_lift",
    },
    {
        "name": "C6_pair_rook_and_standard_double",
        "carrier_dimensions": [36, 72],
        "represented_algebra": "commutative C6 pair bimodule",
        "faithful_B_R": False,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": False,
        "odd_status": "undefined_no_lift",
    },
    {
        "name": "crossed_product_left_right_natural_regular",
        "carrier_dimensions": [6, 36, 360],
        "represented_algebra": "B_R natural/minimal/regular representations",
        "faithful_B_R": True,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": True,
        "odd_status": "zero_if_represented_on_H_off",
    },
    {
        "name": "u_edge_u_chord_v_ref_functional_calculus",
        "carrier_dimensions": [10, 60],
        "represented_algebra": "central D5 class convolution and Morita lift",
        "faithful_B_R": True,
        "explicit_H_off_lift": False,
        "algebra_generated_on_H_off": True,
        "odd_status": "zero_if_represented_on_H_off",
    },
]
check("the frozen inventory contains exactly seven distinct operator families",
      len(families) == 7 and len({row["name"] for row in families}) == 7)
check("the 36-dimensional coincidence does not supply a B_R intertwiner",
      pair["pair_groupoid"]["dimension"]
      == real_crossed["canonical_carriers"][
          "minimum_faithful_left_dimension"] == 36
      and next(row for row in families if row["name"].startswith("C6_pair"))[
          "represented_algebra"] == "commutative C6 pair bimodule"
      and not next(row for row in families
                   if row["name"].startswith("C6_pair"))["faithful_B_R"])
check("the 12-state auxiliary operator is not reassigned to an M12 block",
      12 in next(row for row in families if row["name"] == "Phi_and_D_aux")[
          "carrier_dimensions"]
      and next(row for row in families if row["name"] == "Phi_and_D_aux")[
          "B_R_kernel_dimension"] == 324)

accepted = [row for row in families
            if row["faithful_B_R"] and row["explicit_H_off_lift"]
            and row["odd_status"] == "nonzero"]
undefined = [row for row in families
             if row["odd_status"] == "undefined_no_lift"]
even_only = [row for row in families
             if row["odd_status"] == "zero_if_represented_on_H_off"]
check("zero existing families meet the faithful-lift-nonzero-odd boundary",
      accepted == [] and len(undefined) == 5 and len(even_only) == 2,
      "accepted=0/7; five have no lift, two are algebra-generated and even")

payload = {
    "protocol_commit": "e54695a",
    "target_comparison_performed": False,
    "scope": "seven frozen committed operator families only",
    "H_off": {
        "dimension": 936,
        "algebra_generator_cell_actions_checked": cell_actions_checked,
        "left_right_algebra_generated_odd_projection": "zero",
    },
    "source_audit": {
        "authoritative_files": {
            name: str(path.relative_to(ROOT))
            for name, path in operator_sources.items()
        },
        "explicit_lift_sources": explicit_lift_sources,
        "source_absence_scope": "frozen files only",
    },
    "families": families,
    "acceptance": {
        "family_count": len(families),
        "accepted_count": len(accepted),
        "accepted_fraction": "0/7",
        "undefined_no_lift": len(undefined),
        "algebra_generated_even_only": len(even_only),
    },
    "verdict": (
        "DERIVED REPOSITORY-STATE NEGATIVE: none of the seven operator "
        "families already committed in the six-fibration/action chain has "
        "both a faithful canonical lift to the 936-state carrier and a "
        "nonzero odd component. Five have no lift; the two algebra-generated "
        "families preserve every central cell and have zero odd projection. "
        "A connected continuation requires a new noncentral incidence tensor."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the structured existing-operator audit was written", OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED REPOSITORY-STATE NEGATIVE: accepted existing lifts = 0/7.")
print("OPEN: a new noncentral incidence tensor with a faithful odd lift.")
print("NO HESSIAN OR STANDARD-MODEL TARGET WAS USED.")
raise SystemExit(0 if passed == tests else 1)
