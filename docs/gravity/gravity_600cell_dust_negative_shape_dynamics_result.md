# Result: the negative-stiffness shape carrier is autonomous; expansion remains open

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, TARGET-DISCLOSED.**  The two `15`-dimensional
negative-stiffness shape spaces selected by the preceding blind census are
invariant under both actual centered operators `Gamma` and `Omega`, in both
staircase schedules and all four derivative variants.  They therefore form a
canonical `30`-position subsystem and a regular `60`-state doubled recurrence
on the complete real carrier.

**OPEN:** the midpoint companion spectrum contains `15` multipliers below and
`15` above unit modulus, but the nonnormal Bauer--Fike envelopes include the
unit circle.  The certified outcome is

```text
NEGATIVE_SHAPE_AUTONOMOUS_UNIT_CONSISTENT,
```

not a resolved instability and not a nonexpansion theorem.

## Provenance

| stage | commit |
|---|---|
| prior-art and constraint-framing gate | `db1fa95` |
| target-disclosed protocol | `dde2164` |
| registered verifier before first execution | `2a3cada` |
| disclosed post-result outcome-semantics correction | `af65c64` |
| certified artifact | `f91fcbe` |

Verifier:

```text
reproducible/verify_gravity_600cell_dust_negative_shape_dynamics.py
```

Artifact:

```text
reproducible/gravity_600cell_dust_negative_shape_dynamics.json
SHA-256 51ff46a529958491d7338c62d34612721352c502e926bb4c1df98aa477c9a854
```

Two complete executions after the semantic correction returned byte-identical
artifacts and

```text
12/12 PASS.
```

Only the targeted verifier and its direct geometry import were run.  The full
suite was not run.

## Complete hypotheses

1. the fixed regular 600-cell and its `720` logarithmic signed-squared spatial
   edge coordinates;
2. the first two accepted nonstationary fixed-total-mass dust-Regge slabs;
3. the literal adjacent-slice identification of the centered Jacobi equation;
4. both schedules, four derivative variants and the seven frozen minimal
   binary-tetrahedral sectors;
5. the canonical conformal incidence and action-selected shape complement;
6. sectors `4` and `5`, disclosed in advance because the preceding result had
   already selected their `15 negative + 10 positive` stiffness inertia;
7. the committed midpoint/radius enclosures and `10/100` bands;
8. no exact constraint quotient, dust perturbation carrier, proper-time unit,
   refinement or continuum target.

This is not a blind rediscovery of the sectors and not a result for later
ticks, nonlinear perturbations or a continuum limit.

## Why no constraint quotient was applied

The complete pre-Legendre Jacobian on this fixed finite carrier is already
certified as

```text
rank 1560/1560,
error-consistent nullity 0,
120 resolved-nonzero weak pseudo-constraint candidates.
```

Exact pre/post constraints in variational discrete systems arise from
degenerate directions of the Lagrangian two-form.  Here there is no such
kernel.  Removing the 120 weak directions by a chosen numerical threshold
would change the theory by fitting.

Thus the negative carrier was tested as part of the regular finite dynamics.
A physical quotient remains **OPEN** and would require an independently
derived exact symmetry, extended matter carrier or controlled refinement
limit.

## The selected recurrence

For the shape basis `W`, the preceding result defines

```text
Gamma_S = W* Gamma W,
Omega_S = W* Omega W.
```

The centered equation is

```text
q_2 - 2q_1 + q_0 + Gamma_S(q_2-q_0) + Omega_S q_1 = 0.
```

If `E` spans the negative stiffness eigenspace, the verifier tested

```text
||(I-EE*) Gamma_S E||,
||(I-EE*) Omega_S E||.
```

Only after both passed did it construct the projected recurrence.  No fitted
Schur block was used.

## Invariance is strongly supported

All

```text
2 schedules * 2 sectors * 4 variants * 2 operators = 32
```

residuals are `INVARIANT_CONSISTENT`.

Their raw norms and complete error-unit ratios are

```text
residual norm       3.49e-9 ... 3.43e-8,
residual/error      0.0768  ... 0.3202.
```

The negative/positive stiffness eigengap is more than `12,211` error units in
every cell.  Thus the result is not produced by an ill-defined negative
eigenspace or a permissive invariance threshold.

**DERIVED COMPUTATIONAL:** within the committed enclosure, each selected
`15`-position negative carrier is reducing for both recurrence coefficients.
Together the two dimension-one sectors give the complete `30`-position
autonomous carrier identified by the preceding inertia census.

## Companion regularity

The projected recurrence has the unique state update

```text
C_- = [0, I;
       -(I+Gamma_-)^-1(I-Gamma_-),
        (I+Gamma_-)^-1(2I-Omega_-)].
```

All `16/16` restricted and all `16/16` full-sector matrices have
`REGULAR_RESOLVED` forward coefficients.  The smallest singular value of
`I+Gamma` clears its complete error by

```text
2.05e7 ... 2.28e7 units.
```

The doubled recurrence is therefore not marginally defined.

## What the spectrum does and does not say

For every restricted `30 x 30` companion, the binary midpoint has exactly

```text
15 roots with |mu|<1,
15 roots with |mu|>1.
```

The extrema are schedule- and variant-stable:

```text
minimum midpoint modulus  0.997529105633,
maximum midpoint modulus  1.002481867416.
```

However, the companion is nonnormal.  Its right-eigenvector condition factors
are

```text
4.78e3 ... 9.90e3,
```

and the complete Bauer--Fike eigenvalue errors are

```text
1.61e-3 ... 3.34e-3.
```

Those errors are comparable to or larger than the `2.48e-3` maximum departure
from unit modulus.  Consequently all

```text
16 cells * 30 eigenvalues = 480
```

receive `UNIT_CONSISTENT`; none is resolved expanding or contracting under the
frozen `100`-error rule.

The full `25`-position sectors show the same limitation: all `800` companion
entries are `UNIT_CONSISTENT`.  Their raw midpoint moduli span
`0.997529255056 ... 1.002481717102`.

The midpoint pairing is a **PATTERN**, not a certified instability.

## Why singular amplification is not the missing proof

The restricted companion's largest singular value is

```text
2.414224639877.
```

It is far above one compared with the matrix error.  But the exactly free
second-difference recurrence

```text
q_2 - 2q_1 + q_0 = 0
```

has the companion `[0,I;-I,2I]`, whose largest singular value is exactly

```text
1 + sqrt(2) = 2.414213562373...
```

despite having only the unit eigenvalue in a Jordan block.  The observed
excess over this free baseline is only about `1.108e-5`, or `32.9` companion
matrix-error units: open under the frozen `100`-unit nonzero rule.

Therefore the raw norm `2.414` is primarily a kinematic feature of storing two
successive positions.  It is not evidence of physical transient growth.

## Schedule robustness

All `32/32` available comparisons are `SCHEDULE_ROBUST`:

```text
2 sectors * 4 variants *
(full/restricted companion) * (eigen/singular spectrum).
```

The largest distance/error ratio is below `9.30e-9`.  Both sectors and all
derivative variants reproduce the same midpoint pattern.

## Disclosed semantic correction

The first mechanical outcome was named
`NEGATIVE_SHAPE_AUTONOMOUS_NONEXPANDING_CENSUS`.  That name was too strong:
`UNIT_CONSISTENT` means the error set contains unit modulus, not that it
excludes a slightly expanding true multiplier.

Commit `af65c64` split the last outcome into
`UNIT_CONSISTENT` and an actually all-contracting branch.  No operator,
threshold, residual, eigenvalue label or prior outcome branch changed.  The
scientific interpretation is now the weaker correct one: expansion is OPEN.

## Status ledger

- **DERIVED COMPUTATIONAL:** the two negative-stiffness spaces are autonomous
  under both centered coefficients.
- **DERIVED COMPUTATIONAL:** their complete doubled recurrence is regular and
  schedule robust.
- **PATTERN:** the midpoint roots form `15` inside/`15` outside pairs.
- **OPEN:** whether any true multiplier lies off the unit circle under the
  inherited coefficient uncertainty.
- **REFUTED AS EVIDENCE:** the companion singular norm near `2.414` does not
  establish growth; the free recurrence already has it.
- **OPEN:** physical scalar/vector/tensor identity, dust coupling, long-time
  Lyapunov growth, refinement and continuum interpretation.

This is genuine progress: the negative stiffness is not washed away by mode
mixing.  It defines a clean subsystem.  But its instability has not yet been
certified.

## Post-result prior-art reconciliation

The calculation exposed a structured quadratic matrix-polynomial problem.
Relevant primary numerical work includes:

- Aurentz--Mach--Robol--Vandebril--Watkins on backward-stable matrix-polynomial
  eigenvalue algorithms: <https://arxiv.org/abs/1611.10142>;
- Adhikari--Alam on structured backward errors and pseudospectra for structured
  polynomial eigenproblems: <https://arxiv.org/abs/0907.2545>;
- Noferini--Robol--Vandebril on structured backward errors in companion
  linearizations: <https://arxiv.org/abs/1912.04157>.

These works explain why an unstructured companion Bauer--Fike bound can be
too pessimistic and how polynomial structure may improve the audit.  They do
not decide the present coefficient balls or publish this 600-cell subsystem.
External novelty remains **OPEN**.

## Next load-bearing gate

Do not merely diagonalize the same midpoint at more decimal digits: the
dominant uncertainty comes from the inherited coefficient balls, not LAPACK
roundoff.

The next test should operate directly on the quadratic polynomial

```text
Q(z) = (I+Gamma_-) z^2 + (-2I+Omega_-) z + (I-Gamma_-).
```

It should preregister one of two genuine improvements:

1. a structure-preserving polynomial-eigenvalue enclosure that propagates the
   coefficient balls without the full companion eigenvector condition; or
2. a tighter reconstruction of the source derivative balls sufficient to
   separate `|mu|-1` from zero.

Only a certified root count inside/on/outside the unit circle can settle local
hyperbolicity.  Long-time stability would still require more independently
solved slabs and a non-autonomous cocycle analysis.
