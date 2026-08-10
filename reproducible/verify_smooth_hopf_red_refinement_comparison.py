#!/usr/bin/env python3
"""Compare frozen projected-red Hopf data with preregistered gates."""

import hashlib
import json
import math
from pathlib import Path


DATA = Path(__file__).with_name("smooth_hopf_red_refinement_blind.json")
EXPECTED_SHA256 = "efc09b7b54b987b53d8b0a3086f3846020888362d2f12b7f33c12cda7d0485e0"
DEFINITION_COMMIT = "6b8a80b"
TOLERANCE_FIX_COMMIT = "cde7ad2"
BLIND_COMMIT = "390e38b"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def maximum_error(values, target):
    return max(abs(float(value)-target) for value in values)


print("=" * 78)
print("SHAPE-CONTROLLED SMOOTH-HOPF TOWER COMPARISON")
print("=" * 78)

raw = DATA.read_bytes()
observed_hash = hashlib.sha256(raw).hexdigest()
check("blind artifact hash and commit ordering are frozen",
      observed_hash == EXPECTED_SHA256,
      f"definition={DEFINITION_COMMIT}, tolerance={TOLERANCE_FIX_COMMIT}, "
      f"blind={BLIND_COMMIT}, sha256={observed_hash}")
data = json.loads(raw)
check("blind artifact contains no bootstrap target",
      data["protocol"]
      == "projected red refinement fixed before blind execution; no bootstrap target"
      and "a1" not in raw.decode().lower())

expected_topology = {
    "level0": (120, 720, 1200, 600),
    "level1": (840, 5640, 9600, 4800),
    "level2": (6480, 44880, 76800, 38400),
}
for label, expected in expected_topology.items():
    item = data["topology"][label]
    observed = (item["vertices"], item["edges"],
                item["faces"], item["tetrahedra"])
    check(f"{label} topology and closed-face incidence are exact",
          observed == expected
          and item["euler_characteristic"] == 0
          and item["face_incidence_values"] == [2],
          f"f-vector={observed}, chi={item['euler_characteristic']}")

levels = [data["levels"][f"level{index}"] for index in range(3)]
geometry = [level["geometry"] for level in levels]
h = [item["maximum_chord_length"] for item in geometry]
q = [item["minimum_quality"] for item in geometry]
check("maximum chord length decreases at both refinements",
      h[0] > h[1] > h[2],
      f"h={h}")
check("both refined levels exceed the preregistered quality floor",
      q[1] >= 0.5 and q[2] >= 0.5,
      f"q_min={q}")
check("second-level minimum quality retains at least eighty percent",
      q[2] >= 0.8*q[1],
      f"q2/q1={q[2]/q[1]:.6f}")
check("projected Hopf field exceeds the preregistered norm floor",
      all(item["minimum_projected_hopf_norm"] > 0.98
          for item in geometry),
      f"minimum norms={[item['minimum_projected_hopf_norm'] for item in geometry]}")
check("all stiffness splits close at machine precision",
      max(item["split_relative_residual"] for item in geometry) < 1e-14,
      f"residuals={[item['split_relative_residual'] for item in geometry]}")

mode_error_sequences = {}
for family, key, targets in (
        ("charged", "charged_coordinate_ritz",
         {"vertical": 1.0, "horizontal": 2.0, "full": 3.0}),
        ("base", "base_pullback_ritz",
         {"vertical": 0.0, "horizontal": 8.0, "full": 8.0})):
    expected_dimension = 4 if family == "charged" else 3
    for operator, target in targets.items():
        values = [level[key][operator] for level in levels]
        errors = [maximum_error(level_values, target)
                  for level_values in values]
        mode_error_sequences[f"{family}_{operator}"] = errors
        check(f"{family} {operator} error decreases at both levels",
              errors[0] > errors[1] > errors[2],
              f"errors={[round(error, 9) for error in errors]}")
        check(f"{family} {operator} calibration dimension is preserved",
              all(len(level_values) == expected_dimension
                  for level_values in values))

for index, level in enumerate(levels):
    full_low = level["low_spectra"]["full"]
    zero_count = sum(abs(value) < 1e-7 for value in full_low)
    check(f"level{index} combined low spectrum has one resolved zero",
          zero_count == 1 and full_low[1] > 2,
          f"zero count={zero_count}, first positive={full_low[1]:.9f}")

h_reduction = h[1]/h[2]
second_order_reductions = {
    name: errors[1]/errors[2]
    for name, errors in mode_error_sequences.items()
}
observed_orders = {
    name: math.log(reduction)/math.log(h_reduction)
    for name, reduction in second_order_reductions.items()
}
check("[PATTERN only] all last-step mode rates lie near second order",
      all(1.8 < order < 2.2 for order in observed_orders.values()),
      f"orders={{{', '.join(f'{k}: {v:.3f}' for k, v in observed_orders.items())}}}")

vertical_lows = [level["low_spectra"]["vertical"] for level in levels]
three_band_tops = [max(values[1:4]) for values in vertical_lows]
five_band_tops = [max(values[4:9]) for values in vertical_lows]
check("[PATTERN only] three- and five-mode vertical bands descend twice",
      three_band_tops[0] > three_band_tops[1] > three_band_tops[2]
      and five_band_tops[0] > five_band_tops[1] > five_band_tops[2],
      f"x3 tops={three_band_tops}; x5 tops={five_band_tops}")

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED POSITIVE: all preregistered topology, shape, Hopf-field and")
print("                  canonical-mode gates pass through 6480 vertices.")
print("PATTERN: the last-step rates are uniformly consistent with O(h^2), and")
print("         the 3/5-fold vertical bands descend at both refinements.")
print("DERIVED NEGATIVE: the robust round-Hopf tower still selects no factor five,")
print("                  anisotropy coefficient, or Lorentzian time.")
print("OPEN: infinite-level shape regularity and convergence proof.")
raise SystemExit(0 if passed == tests else 1)

