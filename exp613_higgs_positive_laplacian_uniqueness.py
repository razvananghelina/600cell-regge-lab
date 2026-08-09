"""
exp613_higgs_positive_laplacian_uniqueness.py
=============================================

Strengthening exp612.

exp612 showed:
  - the tree-level Higgs/W ratio phi^2 is realized by the local icosahedral
    Laplacian on the 12-fiber Hopf base;
  - two additional A5-invariant, but nonlocal, distance-mixed operators also
    preserve the same 3'/3 ratio.

This script asks the sharper question:

  If we require a physically admissible Laplacian with nonnegative couplings
  on every A5 distance class, is the local operator still unique?

Answer:
  Yes.

On the 12-node icosahedral base, any A5-invariant Laplacian is

    L(w1,w2,w3) = (5 w1 + 5 w2 + w3) I - w1 A1 - w2 A2 - w3 A3

where
    A1 = nearest-neighbor adjacency,
    A2 = distance-2 adjacency,
    A3 = antipodal matching,
and w_i are the couplings on each relation class.

After fixing the overall scale by w1 = 1, the 3 and 3' eigenvalues are

    lambda_3  = (5-sqrt(5)) + (5+sqrt(5)) w2 + 2 w3
    lambda_3' = (5+sqrt(5)) + (5-sqrt(5)) w2 + 2 w3

Imposing lambda_3'/lambda_3 = phi^2 gives the exact linear condition

    w3 = -5 w2.

Therefore, if one requires nonnegative couplings w2 >= 0 and w3 >= 0,
the only solution is

    w2 = w3 = 0,

namely the unique local nearest-neighbor Laplacian.
"""

from __future__ import annotations

import numpy as np


PHI = (1 + np.sqrt(5)) / 2
SQRT5 = np.sqrt(5)


def spectrum_for_weights(w1: float, w2: float, w3: float) -> dict[str, float]:
    """
    Eigenvalues on the A5 irreducible sectors of the 12-node icosahedron.

    A1 eigenvalues on (1,3,5,3') are (5, +sqrt(5), -1, -sqrt(5)).
    A2 eigenvalues on (1,3,5,3') are (5, -sqrt(5), -1, +sqrt(5)).
    A3 eigenvalues on (1,3,5,3') are (1, -1, +1, -1).
    """
    row_sum = 5 * w1 + 5 * w2 + w3
    return {
        "1": row_sum - 5 * w1 - 5 * w2 - w3,
        "3": row_sum - SQRT5 * w1 + SQRT5 * w2 + w3,
        "5": row_sum + w1 + w2 - w3,
        "3p": row_sum + SQRT5 * w1 - SQRT5 * w2 + w3,
    }


def main() -> None:
    print("=" * 72)
    print("EXP613: UNIQUE POSITIVE LAPLACIAN BEHIND THE HIGGS TREE RATIO")
    print("=" * 72)

    print("\nSECTION 1: General A5-invariant Laplacian on the icosahedral base")
    print("-" * 72)
    print("  L(w1,w2,w3) = (5 w1 + 5 w2 + w3) I - w1 A1 - w2 A2 - w3 A3")
    print("  Normalize the overall scale by w1 = 1.")

    spec_local = spectrum_for_weights(1.0, 0.0, 0.0)
    print("\nSECTION 2: Local operator")
    print("-" * 72)
    print(f"  lambda_3   = {spec_local['3']:.10f} = 5 - sqrt(5)")
    print(f"  lambda_5   = {spec_local['5']:.10f} = 6")
    print(f"  lambda_3'  = {spec_local['3p']:.10f} = 5 + sqrt(5)")
    print(f"  ratio      = {spec_local['3p'] / spec_local['3']:.12f}")
    print(f"  phi^2      = {PHI**2:.12f}")

    print("\nSECTION 3: Analytic phi^2 condition")
    print("-" * 72)
    print("  For w1 = 1:")
    print("    lambda_3  = (5-sqrt(5)) + (5+sqrt(5)) w2 + 2 w3")
    print("    lambda_3' = (5+sqrt(5)) + (5-sqrt(5)) w2 + 2 w3")
    print("  Solving lambda_3'/lambda_3 = phi^2 gives exactly:")
    print("    w3 = -5 w2")

    sample_w2 = 0.25
    spec_line = spectrum_for_weights(1.0, sample_w2, -5.0 * sample_w2)
    print("\n  Check on the phi^2 line w3 = -5 w2:")
    print(f"    sample w2 = {sample_w2:.2f}, w3 = {-5.0 * sample_w2:.2f}")
    print(f"    lambda_3   = {spec_line['3']:.10f}")
    print(f"    lambda_3'  = {spec_line['3p']:.10f}")
    print(f"    ratio      = {spec_line['3p'] / spec_line['3']:.12f}")

    print("\nSECTION 4: Positivity and uniqueness")
    print("-" * 72)
    print("  Physical Laplacian requirement: relation couplings are nonnegative.")
    print("  So with w1 = 1 we require w2 >= 0 and w3 >= 0.")
    print("  But phi^2 forces w3 = -5 w2.")
    print("  Therefore:")
    print("    w2 >= 0 and w3 >= 0  and  w3 = -5 w2")
    print("    => w2 = w3 = 0")
    print("  Hence the only positive phi^2-preserving operator is the local one.")

    print("\nSECTION 5: Numerical scan")
    print("-" * 72)
    print("  Scan domain: w2,w3 in {0, 0.5, 1.0, ..., 3.0}")
    hits: list[tuple[float, float, float, float, float]] = []
    for w2 in [k / 2.0 for k in range(0, 7)]:
        for w3 in [k / 2.0 for k in range(0, 7)]:
            spec = spectrum_for_weights(1.0, w2, w3)
            if spec["3"] <= 1e-9 or spec["3p"] <= 1e-9:
                continue
            ratio = spec["3p"] / spec["3"]
            if abs(ratio - PHI**2) < 1e-9:
                hits.append((w2, w3, spec["3"], spec["5"], spec["3p"]))

    print(f"  Positive hits found = {len(hits)}")
    for w2, w3, l3, l5, l3p in hits:
        print(
            f"    (w2,w3)=({w2:.1f},{w3:.1f})  "
            f"spectrum = 0[1], {l3:.10f}[3], {l5:.10f}[5], {l3p:.10f}[3]"
        )

    assert hits == [(0.0, 0.0, spec_local["3"], spec_local["5"], spec_local["3p"])]

    print("\nVERDICT")
    print("-" * 72)
    print("  exp612 showed that phi^2 survives in three simple A5-invariant operators.")
    print("  exp613 shows that two of them necessarily use negative couplings on a")
    print("  distance class. Once positivity is imposed, they are excluded.")
    print("  So the standard nearest-neighbor icosahedral Laplacian is the unique")
    print("  positive A5-invariant Laplacian that preserves the Higgs tree ratio phi^2.")


if __name__ == "__main__":
    main()
