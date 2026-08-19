# Protocol: invariant correction of the static Lorentz-covariance classifier

Date: 2026-08-19

The original registered verifier returned `12/13` and
`RELATIVE_POINCARE_CONTROL_FAILED`.  Its artifact and source remain frozen.
This protocol changes no original acceptance criterion retroactively.  It
tests the suspected classifier error in a new registered verifier.

## Frozen failed inputs

| input | SHA-256 |
|---|---|
| original protocol | `f88599f2b23d3459f95cba1be12401cd404d98d56a8615f3fa103d1112cc9c7a` |
| original verifier | `308d97cc0b057d3ac79cbc8a4706a63fe1b3d76792d84d8576a84df7d7d63514` |
| failed artifact | `3ac5cce9db2b2f828e0ced2114f301f761dd9371847b712ad47119709396cf7d` |
| failure note | `de81078a7a44af6c49461b0d22d3a59a4edff0e2e4c29432529e00c9f4e39109` |

The new verifier must confirm literally that the failed artifact says
`12/13`, that only the static record has `lorentz_covariance=false`, and that
the original overall outcome remains a control failure.

## Category error under test

A chosen future-timelike unit normal `n` splits `so(3,1)` into spatial
rotations and boosts relative to that observer.  Under a Lorentz
transformation `L`, the normal becomes `n'=L n`.  The invariant statement for
the spatial stabilizer is

```text
so(3)_n     = { A in so(3,1) : A n = 0 },
so(3)_(Ln)  = L so(3)_n L^(-1).
```

It is not invariantly meaningful to require `so(3)_(Ln)` to remain inside
the first three coordinate generators labelled rotations relative to the
old `n`.  The original classifier did require this indirectly by comparing
the separate coordinate rotation/boost projection ranks.

## Exact reconstruction

Independently reconstruct the regular tetrahedron, the three rational
representatives, the six Lorentz generators, the Poincare displacement
matrix `U`, the four-strut constraint matrix `C`, and the rational boost

```text
L = boost_x(cosh=5/4,sinh=3/4).
```

Let `B=diag(L,L,L,L)` act on stacked vertex displacements.  Let `H` be the
exact ten-parameter change

```text
A' = L A L^(-1),
b' = L b.
```

At every representative require exactly

```text
U' H = B U,
C' H = C,
ker(C') = H ker(C).
```

These are the decisive covariance equations.  Total Lorentz-projection rank,
translation-block rank and pure-translation-kernel dimension must be
unchanged.

## Static stabilizer test

At `lambda=1`, compute the Lorentz image of `ker(C)` without naming its
coordinate components.  Require

```text
image_L ker(C)  = ker[A -> A n],
image_L ker(C') = ker[A -> A (L n)],
image_L ker(C') = Ad_L(image_L ker(C)).
```

All three spaces must have exact dimension three.  Additionally require that
the old coordinate rotation/boost rank tuple changes after the boost.  This
is a positive falsifying control for the original non-invariant comparison,
not a defect.

At `lambda!=1`, require both original and transformed Lorentz images to be
all of `so(3,1)` and both translation blocks to be invertible.

## Outcome hierarchy

1. `POINCARE_COVARIANCE_CORRECTION_CONTROL_FAILED` if frozen provenance,
   failed-artifact preservation, exact basis or boost controls fail.
2. `POINCARE_COVARIANCE_REAL_DISAGREEMENT` if either exact intertwining
   equation or kernel transport fails.
3. `STATIC_STABILIZER_AND_EXPANDING_LORENTZ_CHART_CORROBORATED` if the
   covariance equations pass, the static images equal the appropriate
   observer stabilizers and the expanding images are full Lorentz graphs.
4. `POINCARE_COVARIANCE_CORRECTION_OPEN` otherwise.

## Interpretation firewall

A positive correction would establish only:

- **DERIVED EXACT:** the six flexes are a constrained relative-Poincare
  kernel;
- **DERIVED EXACT:** on expanding representatives that kernel is a graph
  over `so(3,1)`;
- **DERIVED EXACT:** on the static representative it instead contains three
  pure translations and projects to the observer-normal stabilizer.

It would **refute**, not rescue, a uniform six-component Lorentz-frame
interpretation through `lambda=1`.  A physical connection, extrinsic
curvature, global gluing, closure, shape matching and dynamics would all
remain **OPEN** and would still require an adversarial local replication
before any global protocol.

Only the correction verifier and static registry guards may be run.
