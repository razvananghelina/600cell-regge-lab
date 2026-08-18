# Result: the unique curvature kernel is localized to the uniform near-minus plane

Date: 2026-08-17

## Provenance

```text
prior-art gate                                      9177531
target-disclosed protocol                           53dc168
angular-resolution guard                            c2cbcd3
registered implementation                           21d0451
preserved 11/14 control-failed first run             c048518
control-correction note                             58d9590
non-duplicated safety-factor clarification           652e86c
corrected implementation                            3e80582
passing artifact                                    f886fa8
```

The first run was retained rather than overwritten.  It failed because the
implementation used an undocumented `1e-6` direct-eigenvector/Schur cutoff
and because two geometrically identical boundary-orbit lists had different
schedule-dependent type orderings.  The correction replaced the accidental
cutoff by the preregistered condition/separation calibration and derived the
unique schedule permutation from equality of literal 24-edge orbit sets.
No kernel component or candidate distance selected that permutation.

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_homogeneous_curvature_kernel.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_homogeneous_curvature_kernel.json
SHA-256 b55887ff3905afd94e86821852d58f0d60c227b52dfbd945044874bfe87540e9
```

Only this targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run.  Corrected result:

```text
14/14 PASS
HOMOGENEOUS_CURVATURE_KERNEL_SUBSPACE_LOCALIZED
elapsed 75.56 s
```

## Exact schedule identification

The two internal sort orders differ by the target-independent orbit-type map

```text
even -> odd:
[0,1,2,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8,9,10,11,21,22,23,24,25,26,27,28,29].
```

Within each matched orbit all 24 literal edge labels have the same group
coordinate.  After this canonical reordering, the two curvature-kernel lines
are `IDENTIFIED`:

```text
principal-line distance   3.10e-12,
calibrated uncertainty    3.71e-5.
```

Thus the one-dimensional nullity is not schedule-dependent at the resolution
of this construction.

## Complete 20-comparison ledger

The preregistered candidate labels were:

```text
IDENTIFIED          4 / 20,
SEPARATED           8 / 20,
NUMERICALLY_OPEN    8 / 20.
```

The same pattern occurred independently for both schedules:

| candidate | even | odd |
|---|---|---|
| near-`-1` contracting eigenline | OPEN | OPEN |
| near-`-1` expanding eigenline | OPEN | OPEN |
| their near-`-1` two-plane | **IDENTIFIED** | **IDENTIFIED** |
| uniform position line | SEPARATED | SEPARATED |
| uniform momentum line | SEPARATED | SEPARATED |
| uniform position/momentum plane | **IDENTIFIED** | **IDENTIFIED** |
| full pure-position subspace | SEPARATED | SEPARATED |
| full pure-momentum subspace | SEPARATED | SEPARATED |
| transported canonical weak lift | OPEN | OPEN |
| transported geometric lapse | OPEN | OPEN |

The uniform-plane distances are `1.5e-12`--`3.1e-12` with uncertainty
`1.85e-5`.  The line is not exactly pure momentum: its distance from the
entire pure-momentum subspace is `0.00343138`, resolved above the same
uncertainty.  It is almost orthogonal to the pure-position subspace, at
`0.9999941`.

The near-`-1` two-plane is cleanly selected from the remainder of the tangent
spectrum: the selection-gap ratio is about 159.44.  Kernel-to-plane distances
are `7.35e-9` and `1.57e-8`; the calibrated plane uncertainties are
`0.00120` and `0.00121`, below the preregistered `0.01` resolution cap.

## What remains numerically open

The individual near-`-1` eigenvectors are badly conditioned in the full
60-dimensional problem.  Although their raw distances from the kernel are
only about `5.2e-6`, their calibrated angular errors are about `0.467`, so
neither individual line is identified.

Likewise, the raw tangent-invariance diagnostics look striking:

```text
distance(K, T K)               3.21e-8,
Rayleigh multiplier            -0.99999376799...,
relative eigenvector residual  7.14e-14.
```

But the conservative 60-dimensional perturbation bound is larger than one,
so the verifier correctly labels `T K = K` **NUMERICALLY_OPEN**.  The weak
Schur-lift and geometric-lapse memberships are also OPEN: their raw distances
are `9.54e-6` and `3.82e-6`, but transport uncertainties are about `0.212`.
Small raw numbers do not override those bounds.

## Scientific verdict

**DERIVED COMPUTATIONAL:** the unique internal-curvature-preserving
boundary-phase direction is the same for both schedule triangulations and is
contained in the exact uniform two-plane spanned by global position-scale and
global momentum-scale variations.

**DERIVED COMPUTATIONAL / target-disclosed:** it is also contained in the
well-separated invariant two-plane of the real reciprocal tangent pair near
`-1`.

**DERIVED NEGATIVE:** it is neither the pure position line nor the pure
momentum line, and it is not contained in either full pure-position or full
pure-momentum half of phase space.

**OPEN:** which individual near-`-1` eigenline it is; whether it is exactly
invariant under the tangent; whether its transported geometry is the
canonical weak lift or the geometric lapse; and whether it integrates to a
nonlinear family.

This does not yet derive time or restore a gauge symmetry.  It reduces the
unresolved object from an arbitrary line in `C^60` to one line in a canonical
real two-dimensional homogeneous phase plane.  That reduction is the actual
advance.

## Framing and prior-art reconciliation

The curved-background Regge literature still prevents calling this line a
diffeomorphism merely because all internal deficits are stationary.  Exact
vertex-displacement gauges are controlled on flat backgrounds, while curved
discretizations generically carry pseudo-constraints.  See
[Hoehn](https://arxiv.org/abs/1411.5672),
[Bahr--Dittrich](https://arxiv.org/abs/0905.1670), and
[Dittrich--Hoehn](https://arxiv.org/abs/0912.1817).

No located primary source identifies this exact 600-cell dust line.  External
novelty remains **OPEN**.

## Next falsifiable gate

The canonical next step is now only a `2 x 2` homogeneous calculation:

1. restrict `F` to the exact uniform position/momentum plane and obtain its
   one-dimensional kernel without the ill-conditioned 60D eigensystem;
2. prove whether the tangent preserves that plane;
3. diagonalize the restricted `2 x 2` tangent at high precision;
4. decide exactly whether the kernel is an eigenline and which reciprocal
   multiplier it carries;
5. only then test its nonlinear integrability and lapse meaning.

If the restricted tangent does not preserve the uniform plane, the attractive
near-`-1` interpretation dies.  If it does and the kernel is an eigenline, the
`119 + 1` split acquires a precise canonical meaning rather than a count-only
pattern.
