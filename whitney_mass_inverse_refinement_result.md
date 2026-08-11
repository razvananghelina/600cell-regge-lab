# Whitney inverse depth grows at the first barycentric refinement

Date: 2026-08-11

Preregistration commit: `366fe4a`

Targeted verifier:
`reproducible/verify_whitney_mass_inverse_refinement.py`

Targeted result: **10/10 PASS**.  The verifier is registered, but the full
suite was not run by explicit user request.

## Result

The exact coarse minimal-polynomial degrees are

\[
(9,22,27,1),
\]

corresponding to inverse-polynomial degrees `(8,21,26,0)`.  On the complete
first barycentric subdivision, deterministic modular Krylov sequences give
the exact lower bounds

\[
\deg m_{\rm fine}\geq(34,60,70,1).
\]

Therefore

\[
\deg M_{\rm fine}^{-1}\geq(33,59,69,0).
\]

> **DERIVED DEGREE GROWTH AT LEVEL 1:** the base inverse depths 8, 21 and 26
> are insufficient after one refinement.  Fixed depth at the base values is
> refuted in every nontrivial cochain degree.

The top-form mass remains a scalar identity and provides the frozen positive
control: its minimal degree stays exactly one and its inverse degree stays
zero.

## Exact carrier and calibration

The verifier constructs every barycentric flag

\[
v\subset e\subset f\subset t
\]

inside all 600 parent tetrahedra.  The exact refined f-vector is

\[
(2640,17040,28800,14400).
\]

The 24 flags of a regular parent tetrahedron have exactly identical
rank-ordered rational Whitney mass blocks.  Those blocks are assembled over
the complete refined complex, denominators are cleared, and each mass is
reduced to a primitive symmetric integer matrix.  No coordinates, mass
lumping or fitted weights enter the Krylov calculation.

The same fixed probe rule, sequence lengths and three primes were applied at
both levels.  On the coarse matrices the method recovers

| form degree | exact coarse degree | complexities at the three primes |
|---:|---:|---:|
| 0 | 9 | 9, 9, 9 |
| 1 | 22 | 22, 22, 22 |
| 2 | 27 | 27, 27, 27 |
| 3 | 1 | 1, 1, 1 |

This is the preregistered calibration.  The refined results are

| form degree | fine dimension | fine complexities | certified lower bound | ratio to coarse degree |
|---:|---:|---:|---:|---:|
| 0 | 2,640 | 34, 34, 34 | 34 | 3.778 |
| 1 | 17,040 | 60, 60, 60 | 60 | 2.727 |
| 2 | 28,800 | 70, 70, 70 | 70 | 2.593 |
| 3 | 14,400 | 1, 1, 1 | 1 | 1.000 |

For a scalar sequence modulo a prime, Berlekamp--Massey complexity (L) is a
rigorous lower bound on the rational matrix minimal degree: every rational
annihilating polynomial reduces to a modular recurrence for that sequence.
The result needs neither a numerical eigensolver nor a generic-probe
assumption.  A probe may miss spectral factors, but that can only make the
reported bound smaller, never falsely increase it.

## Important censoring of the result

The frozen sequence lengths were

\[
(68,120,140,36).
\]

Berlekamp--Massey can return at most half those lengths as a stable lower
bound.  The first three refined complexities are exactly

\[
(34,60,70),
\]

so all three hit the preregistered ceiling.  The actual refined degrees are
not determined here and may be substantially larger.  Reporting the ceiling
as an equality would be false.

## Consequence for the “infinity” question

The preceding fixed-complex result remains correct: every finite mass inverse
is a finite polynomial.  But the newly found short polynomials are not stable
under even one geometric refinement.

What is now established is:

- **DERIVED:** no infinite algebraic operation is needed on either finite
  carrier separately;
- **DERIVED NEGATIVE:** the original finite depths `(8,21,26)` cannot serve as
  a refinement-independent internal tick;
- **PATTERN:** the required spectral complexity grows strongly at the first
  refinement;
- **OPEN:** whether the exact degree diverges along repeated refinements;
- **OPEN:** whether a different coefficient-free unitary construction avoids
  explicit mass inversion with bounded depth.

It would be illegitimate to fit a scaling exponent from two levels.  An
unbounded-depth continuum limit is now plausible, not proved.  Conversely,
bounded depth is no longer supported by the only refinement comparison
available.

## Physical reading

The geometry does not currently provide one universal finite number of
internal microsteps per tick.  If exact Whitney assembly is kept, finer
resolution demands at least more algebraic stages in all propagating form
degrees.  This is evidence that instantaneous metric glue is an effective
coarse description rather than a fundamental one-step process.

It is not yet a derivation of physical time, inertia or (c).  It instead
isolates a concrete obstruction any such derivation must overcome: either
allow scale-dependent internal depth, find a new reversible representation
that bypasses the inverse, or weaken exact Whitney consistency.

## Status ledger

- **DERIVED:** refined f-vector `(2640,17040,28800,14400)`.
- **DERIVED:** exact rational congruence of all 24 rank-ordered child masses.
- **DERIVED:** calibrated modular lower bounds `(34,60,70,1)`.
- **DERIVED:** inverse-degree lower bounds `(33,59,69,0)`.
- **DERIVED NEGATIVE:** base inverse depths fail unchanged at level one.
- **PATTERN:** strong depth growth under first refinement.
- **OPEN:** full refined minimal degrees and repeated-refinement scaling.
- **NOT CLAIMED:** a divergence theorem, physical infinity, time, mass or the
  speed of light.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_mass_inverse_refinement.py
```

Expected result: `10/10`.
