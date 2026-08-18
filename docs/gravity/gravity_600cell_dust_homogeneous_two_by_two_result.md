# Result: the homogeneous plane closes, but the curvature kernel is not a mode

Date: 2026-08-17

## Provenance

```text
prior-art gate                                  24eed99
target-disclosed protocol                       f139f60
registered implementation                       0667c48
mechanical missing-helper correction             d7ada93
passing artifact                                1051329
```

The first execution stopped before reconstructing either operator because the
audited upstream helper `mp_frobenius` was omitted from the explicit import
set.  Commit `d7ada93` adds that one dependency; it changes no operator,
threshold or outcome rule.

Targeted verifier:
`reproducible/verify_gravity_600cell_dust_homogeneous_two_by_two.py`.

Artifact:

```text
reproducible/gravity_600cell_dust_homogeneous_two_by_two.json
SHA-256 d0017d4cfdf3a8833cf19bfcd287b21ac91a7f631c803d5d67114fdf64b77622
```

Only this targeted verifier and its direct 43-control geometry import were
run.  The full suite was not run.  Two complete runs were byte-identical:

```text
17/17 PASS
HOMOGENEOUS_2X2_KERNEL_NOT_EIGENLINE
```

## The exact two-dimensional carrier

The basis was fixed before evaluation:

```text
u_q = (1_30,0_30)/sqrt(30),
u_p = (0_30,1_30)/sqrt(30).
```

For both independent schedule parities, the complete 60-dimensional tangent
maps this plane into itself.  The off-plane norm and its calibrated error are

```text
|| (I-UU*) T U ||_2 = 6.9733927e-30,
epsilon             = 6.9733927e-19.
```

Thus the leakage is eleven orders below its conservative zero threshold.
This is not a compression masquerading as a restriction: the homogeneous
plane is a resolved invariant symplectic subspace.

After the fixed normalization, both schedules give the same matrix to a
Frobenius distance `1.41e-68`, versus calibrated uncertainty `1.23e-19`:

```text
A = U* T U
  = [   -4.9999893222494013187    -0.0137255543024608860 ]
    [ 1165.6889961192268782         2.9999519297508158930 ].
```

Its calibrated invariants are

```text
det A             = 1 + 2.18e-34,
trace A           = -2.0000373924985854257,
symplectic defect = 3.09e-34,
lambda_small      = -0.9939037270641566758,
lambda_large      = -1.0061336654344287499.
```

**DERIVED COMPUTATIONAL:** the one-slab homogeneous linearized dynamics is a
closed real `2 x 2` symplectic hyperbolic map.  The negative sign depends on
the canonical layer identification and is not, by itself, a physical period
or time reversal.  Likewise, off-unit moduli are not yet a continuum
instability without a multi-slab physical norm and clock.

## The curvature-kernel line

The restricted curvature response

```text
B = F U : C^2 -> C^160
```

has calibrated singular values

```text
sigma_min = 1.458e-34,
sigma_max = 0.19202874541896045,
epsilon   = 8.270e-21.
```

It is therefore resolved rank one by more than `2.3e19` uncertainty units.
Its unit kernel vector in the preregistered `(q,p)` basis is

```text
k = (-0.0034313802072921571,
      0.9999941127976069258).
```

The two schedules identify the same line: principal-line distance
`3.90e-45`, versus uncertainty `8.61e-20`.

This independently strengthens the previous 60-dimensional localization.
The unique curvature-preserving direction is genuinely homogeneous and is a
specific mixture of uniform scale and uniform conjugate momentum.

## Decisive negative

The full, not merely compressed, eigenline residual is

```text
||(I-vv*) T v||_2 = 3.2076983283e-8,
epsilon           = 2.008e-16,
v                 = U k.
```

It is nonzero by approximately `1.60e8` uncertainty units.  Its Rayleigh
multiplier is

```text
mu                 = -0.9999937679909155572,
|mu+1|             = 6.2320090844e-6,
epsilon_mu         = 1.004e-16.
```

Hence `mu=-1` is rejected by approximately `6.21e10` uncertainty units.
Both results are identical between schedules.

**DERIVED COMPUTATIONAL NEGATIVE:** the internal-curvature kernel is not an
eigenline of the one-slab canonical tangent and does not carry multiplier
`-1`.  Its earlier raw near-invariance was caused by the small separation of
the two homogeneous reciprocal eigenvalues, not by an exact invariant line.

This breaks the attractive count-only interpretation of

```text
119 strong reciprocal pairs + 1 curvature-kernel line.
```

The `+1` is not a dynamically neutral or separately propagating eigenmode.
It is a line crossing the two-dimensional near-`-1` hyperbolic plane.  Thus
the count does not canonically split the dynamics into `119` physical modes
plus one gauge/time direction.

## What this does and does not say physically

- **DERIVED COMPUTATIONAL:** the uniform scale/momentum plane is an exact
  invariant minisuperspace at linear order for this slab.
- **DERIVED COMPUTATIONAL:** the curvature-preserving line inside it is
  unique and schedule-independent.
- **DERIVED NEGATIVE:** that line is not a one-step eigenmode and is not an
  exact `-1` mode.
- **STRUCTURAL:** this is evidence against identifying the line itself as an
  exact preserved gauge/lapse direction on the fixed curved background.
- **OPEN:** a curved Regge pseudo-constraint need not behave like a flat
  exact gauge generator.
- **OPEN:** the background evolves.  The physically correct propagation test
  is therefore not necessarily `T_n K_n = K_n`, but
  `T_n K_n = K_{n+1}`, where `K_{n+1}` is reconstructed from the next accepted
  slab without fitting.
- **OPEN:** nonlinear integrability, a physical clock, stability under
  refinement and continuum meaning.

The distinction matters.  This result kills the fixed-line/eigenmode story;
it does not yet kill a time-dependent constraint line bundle along the
already accepted multi-tick trajectory.

## Post-result prior-art reconciliation

Hoehn proves preservation of vertex-displacement generators and identifies
curvature observables for *flat-background* linearized Regge evolution:
<https://arxiv.org/abs/1411.5672>.  Bahr--Dittrich and Dittrich--Hoehn show
that curvature generically breaks exact discrete gauge symmetry and replaces
constraints by background-dependent pseudo-constraints:
<https://arxiv.org/abs/0905.1670> and
<https://arxiv.org/abs/0912.1817>.

These sources make the negative unsurprising but do not supply the present
matrix or line.  No located primary source computes this exact invariant
600-cell dust plane.  External novelty remains **OPEN**.

## Next falsifiable gate

Use the already accepted first and second homogeneous slabs.  Independently
reconstruct, before comparison,

```text
K_1 = ker F_1,
K_2 = ker F_2,
T_1 K_1.
```

Then test the target-disclosed covariance equation

```text
T_1 K_1 = K_2
```

with the same four derivative variants, ball radii and schedule controls.
If it passes, the object is a transported time-dependent line rather than an
eigenline.  If it is resolved false, the curvature kernel is not a propagated
constraint/gauge direction even along the accepted homogeneous trajectory.
