#!/usr/bin/env python3
"""Exact continuum audit of the scalar Hopf vertical/horizontal split.

This verifier checks the representation-theoretic spectrum before any new
finite-element observable is chosen.  It also records, without extrapolating,
the kernel counts already committed by the blind first-refinement run.
"""

import json
from pathlib import Path


DATA = Path(__file__).with_name("whitney_hopf_refinement_blind.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def weights(N):
    return tuple(range(-(N-1), N, 2))


def sectors(cutoff):
    for N in range(1, cutoff+1):
        for k in weights(N):
            yield N, k, N  # fixed-weight complex multiplicity


print("=" * 78)
print("CONTINUUM HOPF VERTICAL-INFINITY AUDIT")
print("=" * 78)

cutoff = 40
check("each S3 harmonic level has N Hopf weights",
      all(len(weights(N)) == N for N in range(1, cutoff+1)))
check("fixed-weight multiplicities reconstruct the N^2 scalar multiplicity",
      all(sum(mult for NN, _, mult in sectors(N) if NN == N) == N*N
          for N in range(1, cutoff+1)))

split_exact = True
horizontal_nonnegative = True
for N, k, _ in sectors(cutoff):
    full = N*N-1
    vertical = k*k
    horizontal = full-vertical
    split_exact &= full == vertical+horizontal
    horizontal_nonnegative &= horizontal >= 0
check("Delta=Delta_H+Delta_V on every simultaneous sector", split_exact)
check("horizontal spectrum is non-negative", horizontal_nonnegative)

for level in (1, 2, 3, 8, 17, 40):
    observed = sum(mult for N, k, mult in sectors(level) if k == 0)
    expected = ((level+1)//2)**2
    check(f"vertical kernel count through N<={level} is ceil(N/2)^2",
          observed == expected,
          f"observed={observed}, expected={expected}")

kernel_counts = [
    sum(mult for N, k, mult in sectors(level) if k == 0)
    for level in (10, 20, 30, 40)
]
check("vertical kernel count grows without a finite cutoff plateau",
      kernel_counts == [25, 100, 225, 400],
      f"counts={kernel_counts}")

unit_vertical_multiplicities = [
    sum(mult for N, k, mult in sectors(level) if k*k == 1)
    for level in (10, 20, 30, 40)
]
check("positive vertical eigenvalue one also has growing multiplicity",
      unit_vertical_multiplicities == [60, 220, 480, 840],
      f"counts={unit_vertical_multiplicities}")

horizontal_zeros = [
    (N, k, mult) for N, k, mult in sectors(cutoff)
    if N*N-1-k*k == 0
]
check("horizontal kernel consists only of the constant sector",
      horizontal_zeros == [(1, 0, 1)],
      f"zero sectors={horizontal_zeros}")

horizontal_positive = [
    (N*N-1-k*k, mult, N, k)
    for N, k, mult in sectors(cutoff) if N*N-1-k*k > 0
]
horizontal_gap = min(value for value, _, _, _ in horizontal_positive)
horizontal_gap_mult = sum(
    mult for value, mult, _, _ in horizontal_positive
    if value == horizontal_gap)
check("horizontal gap is two with multiplicity four",
      horizontal_gap == 2 and horizontal_gap_mult == 4,
      f"gap={horizontal_gap}, multiplicity={horizontal_gap_mult}")

for numerator, denominator in ((1, 7), (1, 1), (5, 1), (100, 1)):
    combined_zeros = []
    for N, k, mult in sectors(cutoff):
        # denominator * lambda_H + numerator * lambda_V
        value = denominator*(N*N-1-k*k)+numerator*k*k
        if value == 0:
            combined_zeros.append((N, k, mult))
    check(f"positive combined operator r={numerator}/{denominator} has one zero mode",
          combined_zeros == [(1, 0, 1)])

round_spectrum_ok = all(
    (N*N-1-k*k)+k*k == N*N-1 for N, k, _ in sectors(cutoff)
)
check("round coefficient r=1 reconstructs the scalar S3 Laplacian",
      round_spectrum_ok)

data = json.loads(DATA.read_text())
coarse = data["coarse"]
fine = data["fine"]
check("committed finite vertical kernel stays 12 at the first refinement",
      coarse["fiber_kernel_in_window"] == fine["fiber_kernel_in_window"] == 12)
check("committed finite cross-edge complement keeps only the constant zero mode",
      coarse["cross_kernel_in_window"] == fine["cross_kernel_in_window"] == 1)

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("DERIVED: the continuum vertical operator has infinite kernel and")
print("         infinite positive multiplicities; it is not elliptic.")
print("DERIVED: every positive horizontal+vertical combination has only")
print("         the constant zero mode.")
print("STRUCTURAL: separated component gaps are not physical speeds.")
print("CORRECTION: the old cross-edge tensor is rank three, not the true")
print("            rank-two horizontal Hopf projector.")
print("OPEN: a smooth fiber-adapted 600-cell refinement remains to be built.")
raise SystemExit(0 if passed == tests else 1)
