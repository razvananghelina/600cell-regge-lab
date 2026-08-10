#!/usr/bin/env python3
"""Comparison of preregistered refined Hopf/Whitney gaps with the seed.

Definition commit: be5914a.  Blind refined spectrum commit: 973b512.
The committed JSON is the sole numerical input.  No coefficient is varied.
"""

import hashlib
import json
from pathlib import Path

import sympy as sp


DATA = Path(__file__).with_name("whitney_hopf_refinement_blind.json")
DEFINITION_COMMIT = "be5914a"
BLIND_RESULT_COMMIT = "973b512"
BLIND_SHA256 = "18bb2956b1a9d9acbe46fde8355afbba2b38404ac7273cb142b5ddeaa1b6c6d2"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("=" * 78)
print("REFINED WHITNEY/HOPF GAP COMPARISON")
print("=" * 78)

raw_data = DATA.read_bytes()
check("blind artifact hash is the preregistered one",
      hashlib.sha256(raw_data).hexdigest() == BLIND_SHA256,
      f"sha256={hashlib.sha256(raw_data).hexdigest()}")
data = json.loads(raw_data)
seed = 5
phi = (1+sp.sqrt(5))/2
check("independent bootstrap seed remains exact",
      sp.simplify(1+2*sp.cos(2*sp.pi/seed)-phi) == 0)
check("refinement artifact declares target-blind protocol",
      data["protocol"]
      == "target-blind tensor extension fixed before refined comparison",
      f"definition={DEFINITION_COMMIT}, blind result={BLIND_RESULT_COMMIT}")

coarse = data["coarse"]
fine = data["fine"]
coarse_ratio = coarse["gap_ratio_cross_over_fiber"]
fine_ratio = fine["gap_ratio_cross_over_fiber"]
ratio_drift = fine_ratio/coarse_ratio-1
fiber_drift = fine["fiber_gap"]/coarse["fiber_gap"]-1
cross_drift = fine["cross_gap"]/coarse["cross_gap"]-1

check("coarse generalized ratio agrees with the seed before refinement",
      abs(coarse_ratio-seed) < 1e-8,
      f"coarse ratio={coarse_ratio:.10f}")
check("first-refinement generalized ratio does not equal the seed",
      abs(fine_ratio-seed) > 0.3,
      f"fine ratio={fine_ratio:.10f}")
check("relative ratio drift is larger than six percent",
      ratio_drift > 0.06,
      f"relative drift={100*ratio_drift:.6f}%")
check("fiber and cross gaps renormalize by different amounts",
      abs(fiber_drift-cross_drift) > 0.05,
      f"fiber drift={100*fiber_drift:.6f}%, "
      f"cross drift={100*cross_drift:.6f}%")

compression = data["compression_relative_residuals"]
check("drift occurs despite machine-precision Galerkin compression",
      max(compression.values()) <= 1e-14,
      f"compression residuals={compression}")
check("fiber and cross kernel counts remain resolved in the eigen windows",
      coarse["fiber_kernel_in_window"] == fine["fiber_kernel_in_window"] == 12
      and coarse["cross_kernel_in_window"] == fine["cross_kernel_in_window"] == 1,
      "fiber kernels 12->12; cross kernels 1->1")

coarse_fiber_first = coarse["fiber_low_spectrum"][1]
fine_fiber_first = fine["fiber_low_spectrum"][1]
coarse_cross_first = coarse["cross_low_spectrum"][1]
fine_cross_first = fine["cross_low_spectrum"][1]
check("both first positive clusters retain multiplicity four",
      coarse_fiber_first["multiplicity_in_window"]
      == fine_fiber_first["multiplicity_in_window"]
      == coarse_cross_first["multiplicity_in_window"]
      == fine_cross_first["multiplicity_in_window"] == 4)

check("local tensor rule was positive and coefficient-free",
      data["definition"] == {
          "cross_tensor_minimum_tangent_eigenvalue": 0.499999999963,
          "fiber_tensor_positive_eigenvalue": 0.5,
          "fiber_tensor_rank": 1,
          "local_fiber_edges_per_parent": 1,
      })

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED NEGATIVE: the exact coarse ratio is not stable under the preregistered")
print("                  coefficient-free local-tensor refinement.")
print("DERIVED: the failure is dynamical, not a compression or kernel-count failure.")
print("KILL: within this extension, identifying the first-gap ratio with fixed c^2 is closed.")
print("OPEN: a different observable could survive, but must be preregistered independently.")
raise SystemExit(0 if passed == tests else 1)
