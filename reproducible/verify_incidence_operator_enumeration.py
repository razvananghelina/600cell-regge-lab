#!/usr/bin/env python3
"""Target-independent verifier for the corrected (2,3,5) incidence census.

This supersedes the three-operator calculation in commit 36bd682.  That
calculation averaged a twisted kernel over both full stabilizers, so the sum
factorized and killed every nontrivial character.  A Mackey orbit kernel is
instead a covariant section over (H x K)/C2.  It is nonzero precisely when the
two characters have the same central parity.

No external comparison module is defined or loaded here.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from preregister_orbifold_incidence import enumerate_preregistration


HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / "orbifold_incidence_preregistered.json"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    ok = bool(condition)
    passed += int(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("CORRECTED TARGET-INDEPENDENT (2,3,5) INCIDENCE ENUMERATION")
print("=" * 78)

committed = json.loads(DATA_PATH.read_text(encoding="utf-8"))
regenerated = enumerate_preregistration()
check("[DERIVED] committed census regenerates byte-for-object exactly",
      regenerated == committed)
check("[STRUCTURAL] census declares no external-module comparison",
      committed["status"] == "TARGET_INDEPENDENT_NO_EXTERNAL_MODULE_COMPARISON")

geometry = committed["icosahedral_cell_certificate"]
check("[DERIVED] exact Q(phi) cell construction has f-vector (12,30,20)",
      geometry["cell_counts"] == {"vertices": 12, "edges": 30, "faces": 20})
check("[DERIVED] exact vertex dot products are the icosahedral values",
      geometry["normalized_vertex_dot_products_exact"]
      == ["-1", "-sqrt(5)/5", "sqrt(5)/5"])

relations = geometry["relations"]
relation_signature = [
    (record["source_type"], record["target_type"],
     record["double_coset_count"], record["double_coset_sizes"])
    for record in relations
]
check("[DERIVED] cross-cell double-coset counts and sizes are exact",
      relation_signature == [
          ("E", "V", 6, [20]*6),
          ("F", "V", 4, [30]*4),
          ("F", "E", 10, [12]*10),
      ])
check("[DERIVED] each relation has one pure 60-pair incidence orbit",
      all(sorted(record["pair_orbit_incidence_counts"])
              == [0]*(record["double_coset_count"]-1) + [60]
              for record in relations))
check("[DERIVED] every cross-cell double-coset intersection is the center C2",
      all(record["intersection_orders"] == [2]*record["double_coset_count"]
          for record in relations))

check("[DERIVED] all 20 induced rows span the rank-nine representation ring",
      committed["induction_matrix_rank"] == 9
      and committed["induction_relation_lattice_rank"] == 11)
check("[DERIVED] exact Hom diagonal histogram matches the Gram census",
      committed["hom_diagonal_histogram"]
      == {"3": 8, "4": 2, "7": 4, "8": 2, "15": 2, "16": 2})
check("[DERIVED] exact ordered off-diagonal Hom histogram has no dimension one",
      committed["hom_off_diagonal_ordered_histogram"] == {
          "0": 200, "2": 32, "3": 8, "4": 60, "6": 48,
          "7": 4, "10": 24, "14": 2, "15": 2,
      } and committed["off_diagonal_pairs_with_hom_dimension_one"] == 0)

operators = committed["operators"]
by_relation = Counter((record["source"][0], record["target"][0])
                      for record in operators)
check("[DERIVED] central parity leaves 20+30+12=62 incidence-map lines",
      by_relation == Counter({("E", "V"): 20,
                              ("F", "V"): 30,
                              ("F", "E"): 12})
      and len(operators) == 62)
check("[DERIVED] every map has support only on its incidence double coset",
      all(sum(record["double_coset_basis_support"]) == 1
          and record["double_coset_basis_support"][
              record["incidence_double_coset_index"]
          ] == 1
          and record["incidence_supported_hom_dimension"] == 1
          for record in operators))
check("[DERIVED] two good-prime rank certificates attain every exact upper bound",
      all(len(record["modular_channel_rank_certificates"]) == 2
          and all(certificate == record["channel_ranks"]
                  for certificate in record["modular_channel_rank_certificates"])
          and record["channel_ranks"] == record["channel_upper_bounds"]
          for record in operators))

dims = committed["irrep_dimensions"]
check("[DERIVED] all kernel and cokernel dimensions match their exact characters",
      all(sum(a*b for a, b in zip(record["kernel_irrep_multiplicities"], dims))
              == record["kernel_dimension"]
              and sum(a*b for a, b in zip(record["cokernel_irrep_multiplicities"], dims))
              == record["cokernel_dimension"]
              and record["matrix_rank"] + record["kernel_dimension"]
              == record["matrix_shape"][1]
              and record["matrix_rank"] + record["cokernel_dimension"]
              == record["matrix_shape"][0]
              for record in operators))

complexes = committed["short_complexes"]
check("[DERIVED] exactly one of 60 normalized boundary pairs is a complex",
      committed["counts"]["composable_incidence_pairs_tested"] == 60
      and len(complexes) == 1
      and complexes[0]["modules"] == ["F0", "E2", "V0"]
      and complexes[0]["composition_zero_exactly"])
check("[DERIVED] the cellular complex has H2=1, H1=0, H0=1",
      (complexes[0]["H2_irrep_multiplicities"],
       complexes[0]["H1_irrep_multiplicities"],
       complexes[0]["H0_irrep_multiplicities"])
      == ([1, 0, 0, 0, 0, 0, 0, 0, 0],
          [0, 0, 0, 0, 0, 0, 0, 0, 0],
          [1, 0, 0, 0, 0, 0, 0, 0, 0]))

counts = committed["counts"]
check("[DERIVED] preregistered look-elsewhere census is N=63, 28 indices",
      counts["N_canonical_objects_up_to_adjoint"] == 63
      and counts["distinct_index_characters_up_to_adjoint"] == 28)

print("-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT_OLD_THREE_OPERATOR_ENUMERATION=REFUTED")
print("VERDICT_CANONICAL_INCIDENCE_MAP_PAIRS=62")
print("VERDICT_CANONICAL_SHORT_COMPLEX_PAIRS=1")
print("NO_EXTERNAL_MODULE_COMPARISON_PERFORMED")
if passed != tests:
    raise SystemExit(1)
