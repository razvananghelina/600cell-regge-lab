# Finite-height full-boundary canonical tangent

Date: 2026-08-22

## Verdict

**DERIVED COMPUTATIONAL, MECHANICALLY DIFFERENTLY REPLICATED, UNDER THE
HYPOTHESES BELOW:** the first accepted positive-height homogeneous 600-cell
slab has a unique linearized canonical boundary response

```text
(delta old geometry, delta incoming momentum)
    -> (delta new geometry, delta outgoing momentum)
```

on the complete 1440-dimensional boundary phase space.  The map is
symplectic inside the frozen numerical envelopes and no dependence on the
two H4 staircase schedules is resolved.

This is the first accepted full forced response after the internal-kernel
calculation.  It is not a free internal mode and it is not an eigenmode or
wave equation.

## Complete hypotheses

The statement applies only to:

- the fixed regular 600-cell spatial carrier;
- the zero-cosmological-constant Lorentzian Regge action with the accepted
  boundary term and conserved homogeneous global dust;
- the representative incoming state `v=3/2` and its first accepted
  positive-height slab, reconstructed from the exact homogeneous equations;
- positive `h`, `lambda`, `rho` and `lambda-rho` on the already certified
  Lorentzian branch;
- logarithmic signed-squared-edge variables;
- all 720 old, 840 internal and 720 new physical edge variables;
- canonical boundary momenta obtained from the action's pre- and
  post-Legendre transforms;
- the exact physical final-to-old edge identification
  `{u,v}->{u+120,v+120}`;
- the two frozen even and odd H4 staircase schedules;
- a first-order tangent at this one background slab.

The statement does not cover a generic incoming state, a nonlinear
nonhomogeneous evolution, another slab, refinement or a continuum limit.

## Preregistered chronology

```text
c3fc22d  prior-art gate
2f6a4a5  primary protocol
aae40c7  primary verifier registration
3f0dd26  primary implementation before first scientific run
ff6e45a  preserved primary control failure
7ad0bf1  scalar sign-control correction only
88833b0  primary artifacts: 21/21
3be6eb6  adversarial dense protocol
b3ab177  adversarial verifier registration
fb80285  adversarial implementation before first scientific run
2516b3d  preserved adversarial bookkeeping failure
9445452  odd hostile-control bookkeeping correction only
6139101  adversarial artifact: 22/22
```

The adversarial code did not open a primary tangent entry or singular value
until its dense rank, canonicality and parity labels had been frozen in
memory.

## Primary construction

The primary verifier assembled a 95-by-95 group-convolution Hessian kernel
at 180 decimal digits, projected it to the seven deterministic minimal
right-regular sectors of the binary tetrahedral group, and classified every
pre-Legendre block by both a scaled binary64 SVD gap and a 140-decimal Flint
determinant ball.  It then solved 42 complex ball systems and reconstructed
the canonical sector maps.

Result:

```text
21/21 PASS
even: 7/7 REGULAR, 7/7 CANONICAL
odd:  7/7 REGULAR, 7/7 CANONICAL
schedule: SCHEDULE_ROBUST
outcome: FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY
```

The primary numerical archive SHA-256 is
`0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf`.

## Mechanically different adversarial construction

The adversarial verifier did not use the primary sector assembly, sector
projection, Flint systems or representation-theoretic canonicality test.  It
instead:

1. assembled four complete real 2280-by-2280 Hessians for each schedule in
   physical edge coordinates;
2. formed three dense Richardson levels;
3. classified each complete real 1560-by-1560 pre-Legendre matrix;
4. solved each regular system with all 1440 right-hand sides;
5. built the complete real 1440-by-1440 tangent map;
6. checked the three real symplectic block identities directly;
7. compared the schedules in the common physical edge ordering;
8. only then compared its complete singular-value multiset with the union of
   the primary minimal-sector singular values.

Result:

```text
22/22 PASS
even: REGULAR, CANONICAL, PRIMARY_AGREES
odd:  REGULAR, CANONICAL, PRIMARY_AGREES
schedule: SCHEDULE_ROBUST
outcome:
  FINITE_HEIGHT_FULL_BOUNDARY_TANGENT_
  SCHEDULE_ROBUST_ADVERSARIALLY_REPLICATED
```

Artifact SHA-256:
`ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7`.

## Numerical margins

The complete dense pre-Legendre normalized smallest singular values are

```text
even  4.0880740195610725e-4
odd   4.0880745009547740e-4
gate  3.6510838476029892e-7
```

Thus both rank gaps exceed the frozen gate by about 1120.  The corresponding
condition estimates are 114.28 and 114.32.

The raw dense Hessian reciprocity residuals, evaluated before the licensed
symmetrization, are about `4.61e-14`, against the conservative Hessian
envelope `3.65e-9`.  All three dense Richardson levels became identical in
binary64 at the frozen very small local steps, so the reported envelope is
roundoff-dominated rather than a measured truncation extrapolation.

The direct full-map symplectic defects are

```text
even  2.7033e-10
odd   2.7062e-10
```

against the canonical acceptance threshold `10*e_sym = 5.37e-2`.  This is a
very conservative envelope.  The hostile omission of the actual even
`K_NO` term gives defect `5.919e1`, and a cyclicly corrupted output-edge map
differs by `1.246` in the normalized Frobenius comparison.

The direct even--odd distances are `4.668e-12` at all three levels, against
schedule uncertainty `3.651e-9`.  The primary/adversarial singular-spectrum
distances are `5.044e-12` and `4.940e-12`.  They are approximately 1.99 and
1.95 times their conservative post-classification comparison uncertainties,
but below the preregistered `10*uncertainty` agreement boundary.  Therefore
the honest label is **PRIMARY_AGREES under the frozen gate**, not bitwise
identity of the two arithmetic routes.

## What this proves

- **DERIVED COMPUTATIONAL:** the complete pre-Legendre system is regular at
  the stated finite-height slab.
- **DERIVED COMPUTATIONAL:** fixed incoming geometry and momentum determine
  one first-order outgoing geometry and momentum.
- **DERIVED COMPUTATIONAL:** the action-generated map obeys the canonical
  symplectic identities within the frozen envelopes.
- **DERIVED COMPUTATIONAL:** no staircase-schedule dependence is resolved.
- **DERIVED COMPUTATIONAL / STRUCTURAL:** a group-reduced complex-ball route
  and a full dense real-space route agree on the unitary-invariant singular
  spectrum under their preregistered comparison gate.

## What this does not prove

- **OPEN:** which tangent directions are physical after the appropriate
  discrete constraint or pseudo-constraint reduction;
- **OPEN:** a second nonhomogeneous tick or a map along the later homogeneous
  history;
- **OPEN:** a canonical identification of tangent fibres at different
  background scales suitable for physical eigenvalues;
- **OPEN:** a discrete graviton equation, dispersion relation or stability
  theorem;
- **NOT DERIVED:** a limiting speed, physical tick, `c`, `G`, Planck scale,
  particle masses or Standard-Model content;
- **OPEN:** external novelty beyond the repository prior-art search;
- **OPEN:** convergence, infinite total proper duration and completeness of
  the homogeneous history.

In particular, the invariant-region theorem proves by induction a unique
physical successor at every later finite step for the representative
homogeneous branch.  It is not a theorem of a unique infinite evolution.

## Next falsifiable gate

Do not diagonalize this isolated one-step map and call its eigenvalues
particles or waves.  It maps tangent spaces at different background scales,
and a one-step spectrum depends on how those fibres are identified.

The next calculation must first derive, from the scale covariance and the
canonical symplectic form, a co-moving canonical trivialization of the
boundary phase space.  It must then construct the complete tangent on the
next accepted homogeneous slab and compose the two maps.  Only quantities
invariant under the preregistered admissible trivializations may be assigned
physical mode labels.  Failure to obtain a unique such reduction keeps the
mode interpretation **OPEN**; a schedule dependence or loss of regularity on
the second slab refutes the proposed local dynamical continuation.

No full suite was run.  The targeted primary and adversarial verifiers are
registered exactly once; at the adversarial run the registry contained 241
entries and 241 distinct names.

