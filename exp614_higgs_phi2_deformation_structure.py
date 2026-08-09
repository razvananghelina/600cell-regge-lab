"""
exp614_higgs_phi2_deformation_structure.py
==========================================

Follow-up to exp612 and exp613.

exp612 found three simple A5-invariant operators on the 12-node icosahedral
base that preserve the ratio

    lambda_3' / lambda_3 = phi^2.

exp613 showed that positivity excludes the two nonlocal ones, leaving the
nearest-neighbor Laplacian as the unique positive solution.

This script asks a structural question:

    What do the nonlocal phi^2-preserving deformations actually do?

Answer:
  They do not mix the 3 and 3' sectors in any complicated way.
  The exact phi^2 condition is

      w3 = -5 w2

  and along this line the eigenvalues become

      lambda_1  = 0
      lambda_3  = (1 - w2) (5 - sqrt(5))
      lambda_3' = (1 - w2) (5 + sqrt(5))
      lambda_5  = (1 + w2) * 6

So every phi^2-preserving deformation does exactly two things:
  1. rescales the Galois pair (3,3') by a common factor (1 - w2),
  2. independently rescales the 5-sector by (1 + w2).

Therefore the Higgs/W tree ratio is protected by the equal-dimensional
Galois pair itself; the only additional freedom acts through the orthogonal
5-sector. The local operator is then selected by positivity (exp613).
"""

from __future__ import annotations

import numpy as np


PHI = (1 + np.sqrt(5)) / 2
SQRT5 = np.sqrt(5)


def spectrum_for_weights(w1: float, w2: float, w3: float) -> dict[str, float]:
    row_sum = 5 * w1 + 5 * w2 + w3
    return {
        "1": row_sum - 5 * w1 - 5 * w2 - w3,
        "3": row_sum - SQRT5 * w1 + SQRT5 * w2 + w3,
        "5": row_sum + w1 + w2 - w3,
        "3p": row_sum + SQRT5 * w1 - SQRT5 * w2 + w3,
    }


def main() -> None:
    print("=" * 72)
    print("EXP614: STRUCTURE OF THE phi^2-PRESERVING DEFORMATION LINE")
    print("=" * 72)

    print("\nSECTION 1: General family")
    print("-" * 72)
    print("  L(w1,w2,w3) = (5 w1 + 5 w2 + w3) I - w1 A1 - w2 A2 - w3 A3")
    print("  Fix scale by w1 = 1.")

    print("\nSECTION 2: Exact phi^2 condition")
    print("-" * 72)
    print("  Imposing lambda_3'/lambda_3 = phi^2 gives:")
    print("    w3 = -5 w2")
    print("  Substitute this back into the spectrum.")

    sample_values = [-0.5, -0.25, 0.0, 0.25, 0.5]

    print("\nSECTION 3: Closed-form spectrum on the phi^2 line")
    print("-" * 72)
    print("  Along w3 = -5 w2:")
    print("    lambda_1  = 0")
    print("    lambda_3  = (1 - w2)(5 - sqrt(5))")
    print("    lambda_3' = (1 - w2)(5 + sqrt(5))")
    print("    lambda_5  = (1 + w2) * 6")
    print()
    print("  Sample points:")
    for w2 in sample_values:
        w3 = -5.0 * w2
        spec = spectrum_for_weights(1.0, w2, w3)
        ratio = spec["3p"] / spec["3"]
        print(
            f"    w2={w2:+.2f}, w3={w3:+.2f}  "
            f"lambda_3={spec['3']:.10f}, lambda_5={spec['5']:.10f}, "
            f"lambda_3'={spec['3p']:.10f}, ratio={ratio:.12f}"
        )

    print("\nSECTION 4: Interpretation")
    print("-" * 72)
    print("  The entire phi^2-preserving family acts diagonally on the A5 sectors.")
    print("  It never distorts the 3/3' ratio internally; it only rescales")
    print("  the whole Galois pair together, while independently moving the 5-sector.")
    print("  So the tree-level electroweak ratio is robust inside the equal-dimensional")
    print("  pair itself. The residual ambiguity lives in the orthogonal 5-sector.")
    print("  exp613 then removes that ambiguity by positivity, which forces w2 = 0.")

    print("\nVERDICT")
    print("-" * 72)
    print("  The phi^2 ratio is structurally protected by the 3/3' Galois pair.")
    print("  Nonlocal deformations do not create a new ratio; they only rescale the")
    print("  pair and shift the 5-sector. Positivity selects the local operator.")


if __name__ == "__main__":
    main()
