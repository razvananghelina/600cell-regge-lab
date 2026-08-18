# Prior-art and framing gate: residual-certified generalized bundle

Date: 2026-08-18

## Exact object, carrier and hypotheses

The fixed physical input is the committed three-slab homogeneous
Lorentzian Regge--dust trajectory on the regular 600-cell, with conserved
total dust mass.  The fixed numerical family consists of the four already
disclosed complex-step derivative schedules

```text
operational_primary, operational_shadow,
validation_primary,  validation_shadow.
```

For each of two centered times, two schedule parities and two disclosed
one-dimensional symmetry sectors, the action supplies Hermitian matrices on
the literal 30-orbit edge carrier.  The action-selected shape space is

```text
S_t = ker(U* M_t),       dim S_t = 25,
```

where `U` is the rank-5 conformal incidence image.  On `S_t` the generalized
Hermitian-definite pencil is

```text
A_t = -V_t,       B_t = -M_t > 0,
A_t x = lambda B_t x.
```

The preceding direct audit found a rank-15 lower generalized spectral cluster
in every one of the `32` cells.  Its old/shifted Euclidean projector distance
is approximately `3.96e-7`, but the generic whole-pencil perturbation bound is
approximately `2.43e-4 ... 6.46e-4` per projector.  Error attribution shows
that this generalized-eigenspace term dominates all `32/32` cells.

This mission asks whether a high-precision, a posteriori residual calculation
classifies the old and shifted finite-family fibers as robustly rotated,
zero-consistent, or open.  It does not ask for roots, a propagator, a wave
speed, mass or a physical particle label.

## Primary literature

**KNOWN.** Davis and Kahan treat both perturbation and computable-residual
bounds for invariant subspaces of Hermitian operators, with the subspace angle
controlled by a separating gap:
<https://doi.org/10.1137/0707001>.

**KNOWN.** The block-off-diagonal form of the Davis--Kahan `tan(2 Theta)`
theorem gives the sharp projector estimate

```text
||P-P0|| <= sin(0.5 atan(2 ||R|| / d)),
```

when a Hermitian operator is split into two diagonal spectral blocks separated
by `d>0` and an off-diagonal coupling of norm `||R||`:
<https://arxiv.org/abs/math/0302020>.

**KNOWN.** A Hermitian-definite generalized problem with `B>0` reduces by a
positive Cholesky factor to an ordinary Hermitian problem; this is the standard
LAPACK construction:
<https://netlib.org/lapack/lug/node54.html> and
<https://doi.org/10.1137/1.9780898719581.ch5>.

**KNOWN.** Perturbation theory for partitioned Hermitian-definite pencils
distinguishes diagonal changes, which chiefly move eigenvalues, from
off-diagonal changes, which rotate invariant subspaces:
<https://doi.org/10.1137/100808356>.

**KNOWN.** Stewart's residual analysis takes an approximate invariant
subspace `X`, its Rayleigh compression `X* H X`, and residual
`R=H X-X(X* H X)` as the a posteriori objects; the cited paper bounds spectral
errors, while the Davis--Kahan paper above supplies the subspace-angle result:
<https://doi.org/10.1137/0612016>.

These sources fix the linear-algebraic certification.  They say nothing about
the present Regge operator, its two-time rotation or its physical meaning.
External novelty is **OPEN**.

## Framing correction

A numerical residual bound cannot prove that two analytically defined
subspaces are exactly equal.  At finite precision it can establish only one
of the following:

1. every admissible old/new projector pair is separated by more than its
   certified error, so finite-family rotation is resolved;
2. every observed separation is smaller than the zero-consistency threshold;
3. the result remains open.

Therefore the old mechanical label `COMMON_BUNDLE_RESOLVED` must not be read
as an equality theorem.  Exact equality would require an independent
structural identity, for example an exactly vanishing common reducing
projector relation for both pencils.  No such identity is currently known.

The four finite-difference schedules also do not constitute a rigorous
analytic enclosure of the continuum Hessian.  The new calculation may certify
the Arb source balls and solver residual for each schedule separately and may
show agreement across the complete frozen family.  Its result remains
conditional on that family until an analytic/automatic derivative enclosure
is built.

## Proposed difference from the previous bound

The previous estimate replaced the complete derivative-family variation by
one unstructured operator-norm ball.  Such a ball allows the entire error to
sit in the off-diagonal block and therefore pays roughly
`whole-pencil error / spectral gap`, even when the observed schedules rotate
the subspace far less.

The proposed audit instead:

- keeps all four schedules as four disclosed operators rather than collapsing
  them into an arbitrary unstructured perturbation;
- reconstructs the conformal and shape carriers at high precision;
- Cholesky-whitens each Hermitian-definite pencil independently;
- certifies each computed invariant subspace from its actual off-diagonal
  residual and block separation;
- compares all old/new schedule pairs, so a result cannot depend on choosing
  the most favorable derivative step.

This is an a posteriori structured-family test, not a fitted smaller error bar.

## KNOWN, CONTROL, OPEN

- **DERIVED UPSTREAM:** every one of the `32` finite-family pencils has a
  resolved `15+10` generalized spectral split.
- **DERIVED UPSTREAM:** the old/shifted midpoint projector displacement is
  approximately `3.96e-7` in every matched cell.
- **CONTROL:** all input files and source hashes remain frozen.
- **CONTROL:** the lexicographic conformal-image and shape-nullspace
  constructions must have the certified ranks `5` and `25` without looking at
  old/new projector distances.
- **CONTROL:** every restricted kinetic form must remain positive and every
  lower/upper block gap must remain positive before a subspace is classified.
- **CONTROL:** the high-precision projectors must overlap the previously
  committed binary64 projectors within the old broad errors.
- **OPEN:** whether all cross-schedule old/new pairs resolve a nonzero
  rotation.
- **OPEN:** whether the analytic Hessian, outside the disclosed finite family,
  has the same rotation.
- **OPEN:** whether the action supplies a connection between distinct fibers.
- **FORBIDDEN HERE:** fitted alignment, Procrustes/polar transport, reduced
  propagators, root counts, dispersion, `c`, mass or particle inertia.
