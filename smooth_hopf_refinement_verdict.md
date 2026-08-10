# Smooth orthogonal Hopf refinement: comparison and verdict

Date: 2026-08-10

## Provenance

The geometry and all comparison gates were fixed before the mode data:

- definition commit: `5f78826`;
- blind-result commit: `5202966`;
- blind JSON SHA-256:
  `7258a4755ac32af9d32d2415c09bf65f7b1c4a064475a0e2da304f43eb362ba8`.

The blind verifier contains no bootstrap integer and printed no mode values
before the JSON was committed.

## Convention correction

The preregistration wrote the representative generator as
`cos(pi/5)+u sin(pi/5)`.  The deterministic group enumeration actually chose

`g=(-0.30901699,-0.80901699,-0.5,0)`,

whose angle is `3*pi/5` up to the sign convention for `u`.  This is another
primitive generator of the same order-ten subgroup.  It has the same axis and
the same twelve right-coset fibers.  The code extracts
`u=Im(g)/|Im(g)|`, so the smooth field `X(q)=q*u` is unaffected.  **DERIVED
CONVENTION CORRECTION; no numerical choice changes.**

## Geometric and algebraic gates

The projected-barycentric mesh changes as follows:

```text
                              coarse        first refinement
vertices                         120                    2640
tetrahedra                       600                   14400
maximum chord length        0.618034                0.385708
min |projected X|            1.000000                0.996772
split residual              1.85e-16                 7.68e-16
```

At both levels the vertical projector has spectrum `(0,0,0,1)` in ambient
four-space and the horizontal projector has `(0,0,1,1)`.  Thus their tangent
ranks are exactly one and two, and their sum is the full tangent projector.
All algebraic and geometric preregistration gates pass.

## Canonical continuum-mode comparison

### Fiber-charged coordinate modes

The exact round-`S3` target is `(lambda_V,lambda_H,lambda_full)=(1,2,3)`.

```text
operator       coarse Ritz       fine Ritz       |fine-target|
vertical       1.129420018       1.024886392        0.024886392
horizontal     2.258840035       2.049772784        0.049772784
full           3.388260053       3.074659176        0.074659176
```

Each value has multiplicity four.  Every error decreases by a factor of about
`5.20`.

### Fiber-invariant base modes

The three components of the Hopf map have exact target `(0,8,8)`.

```text
operator       coarse Ritz       fine Ritz
vertical       0.211145618       0.066509143
horizontal     9.788854382       8.361166971
full          10.000000000       8.427676114
```

The vertical leakage falls by a factor of about `3.17`.  Horizontal error
falls from `1.78885` to `0.36117`, and full error from `2` to `0.42768`.
All preregistered convergence-direction gates pass.

The low vertical spectrum also begins

`0; 0.009822 (x3); 0.111287 (x5); ...`.

The multiplicities `1,3,5` are those of the first scalar harmonics on the
Hopf base.  Because no clustering tolerance or multiplicity acceptance rule
was preregistered, this is labelled **PATTERN**, not used as a pass condition.
The separately preregistered base-mode Ritz test is the load-bearing evidence.

## What happened to the infinity?

It appears in the correct place at the level actually tested.  The three
explicit fiber-invariant base functions acquire substantially smaller
vertical energies.  Additional apparent base harmonics enter a near-zero band,
but that second statement is only the unregistered `1,3,5` pattern above.
The analytic continuum theorem, not this one refinement, proves the
infinite-dimensional kernel.

It does not contaminate the combined geometry.  The full operator has exactly
one resolved zero mode at both levels and its first four positive modes move
from `3.38826` to `3.06801`, toward the round-`S3` value `3`.

Thus:

- **DERIVED:** the continuum infinity is real, and the three preregistered
  base modes show decreasing vertical leakage at the first refinement;
- **DERIVED:** the positive combined operator removes it down to constants;
- **DERIVED:** the smooth projected `P1` implementation passes every
  preregistered first-refinement calibration gate;
- **PATTERN:** the unregistered `1,3,5` near-zero cluster multiplicities are
  encouraging but not acceptance evidence.

## Consequence for `a_1=5`

The old exact five is now fully explained.  It belongs to the combinatorial
split of one fiber edge against the other five tetrahedral edges.  Its
so-called cross tensor has tangent spectrum `(1/2,1,1)` and is not the
rank-two horizontal Hopf projector.

The continuum-faithful split instead recovers the standard round-sphere
numbers `(1,2,3)` and `(0,8,8)`.  It supplies no dynamically selected five.
Inserting a coefficient `r=5` into `K_H+r K_V` would be an external choice and
is forbidden.

**DERIVED NEGATIVE:** the proposed identification of `a_1=5` with a
vertical/horizontal propagation-speed ratio is closed.  The earlier equality
was exact but combinatorial, not a continuum kinetic constant.

**DERIVED POSITIVE:** the projected first refinement gives a
continuum-consistent discretization of the genuine Hopf vertical/horizontal
geometry on all preregistered calibration modes.  This is not yet a convergence
theorem for an infinite tower: repeated-refinement shape regularity remains a
separate gate.  The result connects to `a_1=5` only through the already-known
choice of the decagonal fibration, not through a new physical constant.

## Next admissible physics gate

The round metric selects `r=1`.  A non-round Berger coefficient requires an
independent action, symmetry breaking mechanism, or variational selector.
Until such a selector exists, the honest dynamical operator is the combined
round one, and Lorentzian time remains **OPEN**.  The compact Hopf circle must
not be relabelled as time: doing so would introduce closed timelike curves
rather than derive causality.
