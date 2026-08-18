# Raw commuting recovery is topologically impossible

Date: 2026-08-12

Preregistration commit: `39aa5c8`

Targeted verifier:
`reproducible/verify_whitney_commuting_recovery.py`

Targeted result: **8/8 PASS**.  The verifier is registered exactly once.  No
spectrum, fit or phenomenological target was used.  The full suite was not
run, by explicit user request.

## Headline

The proposed cochain-map axiom cannot select a broken-FEEC recovery.  It is
incompatible with being a left inverse of the occurrence injection.

The preregistered local support argument proves this for every
strict-occurrence-local recovery, allowing arbitrary signed weights.  A
post-result cohomology argument strengthens it substantially:

> **DERIVED GLOBAL NO-GO:** there is no recovery of any support radius that
> is simultaneously a cochain map from the raw direct-sum tetrahedron
> complex and a left inverse of the conforming occurrence injection.

This is not a failure of broken FEEC.  It explains why its differential is
defined using projection, `d_h=D_pw P`, instead of requiring recovery to
commute with `D_pw`.

## Complete hypotheses

For

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

let `V_p` be the direct sum of local cochains on all tetrahedra, `W_p` the
global conforming cochains, `J:W->V` the signed occurrence injection,
`D_pw` the elementwise coboundary and `d` the global coboundary.  The
injection is an exact cochain map:

\[
 D_{\rm pw}J=Jd.
\]

The tested recovery conditions are

\[
 LJ=I,
 \qquad LD_{\rm pw}=dL.
\]

No positivity, symmetry, metric weighting, equivariance or candidate formula
is assumed in the global theorem.

## Preregistered strict-support theorem

For a strict-occurrence-local recovery, a local copy `(T,s)` can contribute
only to the same global simplex `s`.  In the column belonging to this copy,
`L D_pw` can see only cofaces `u` satisfying

\[
 s\subset u\subset T.
\]

But `d L` has coefficient `plus_or_minus w_(T,s)` at every global coface of
`s`.  Any coface outside `T` therefore gives the exact equation

\[
 w_{T,s}=0.
\]

The exhaustive census is:

| `k` | degree | occurrences forced to zero | external cofaces per occurrence |
|---:|---:|---:|---:|
| 1 | 0 | 480 / 480 | 5--11 |
| 1 | 1 | 720 / 720 | 2--4 |
| 1 | 2 | 480 / 480 | 1 |
| 2 | 0 | 3,840 / 3,840 | 5--11 |
| 2 | 1 | 5,760 / 5,760 | 2--4 |
| 2 | 2 | 3,840 / 3,840 | 1 |
| 4 | 0 | 30,720 / 30,720 | 5--11 |
| 4 | 1 | 46,080 / 46,080 | 2--4 |
| 4 | 2 | 30,720 / 30,720 | 1 |

Thus every admissible weight is zero, while `LJ=I` requires the occurrence
weights above each global simplex to sum to one.  The contradiction holds in
each degree separately.  **DERIVED on all frozen controls.**

This proof covers arbitrary real weights, including negative and asymmetric
ones.  Failure of the two named candidates is not used as evidence for the
general claim.

## Exact candidate controls

For completeness, the verifier constructs the exact commutators

\[
 C^X_p=L^X_{p+1}D_{{\rm pw},p}-d_pL^X_p
\]

for equal counting and diagonal-Whitney recovery:

| `k` | degree | counting: nonzeros / maximum | diagonal: nonzeros / maximum |
|---:|---:|---:|---:|
| 1 | 0 | 5,280 / `5/24` | 5,280 / `5/24` |
| 1 | 1 | 3,600 / `1/3` | 3,600 / `1/3` |
| 1 | 2 | 960 / `1/2` | 960 / `1/2` |
| 2 | 0 | 50,400 / `5/24` | 50,400 / `5/24` |
| 2 | 1 | 30,240 / `1/3` | 30,240 / `701/1312` |
| 2 | 2 | 7,680 / `1/2` | 7,680 / `21/32` |

All entries are exact rational numbers.  Neither the size ordering nor the
maximum coefficient selects a candidate.

## Post-result global strengthening

The local proof initially left open a recovery with larger support.  That
opening is closed by cohomology, and this strengthening is labelled
post-result because it was not a frozen decision gate.

The raw broken complex is a direct sum of the cochain complexes of the
individual tetrahedra.  One tetrahedron has exact differential ranks

\[
 (3,3,1)
\]

and Betti vector

\[
 (1,0,0,0).
\]

At `k=1`, the 120-fold direct sum therefore has

\[
 H^3(V)=0.
\]

The conforming closed complex has dimensions `(30,150,240,120)`.  Exact
integer nilpotency and modular ranks at two primes give

\[
 \operatorname{rank}(d_0,d_1,d_2)=(29,121,119).
\]

These are exact rational ranks, not probabilistic estimates: modular ranks
are lower bounds over the rationals, while the constant zero-cochain and
successive nilpotency give the identical upper bounds `(29,121,119)`.  Hence

\[
 H^\bullet(W)=(1,0,0,1).
\]

If `L` and `J` were cochain maps with `LJ=I`, they would induce

\[
 L_*J_*=I
\]

on cohomology.  In degree three, however,

\[
 J_*:\mathbb C\cong H^3(W)\longrightarrow H^3(V)=0,
\]

so `L_*J_*` is zero and cannot be the identity.  Contradiction.

This argument is independent of recovery support, metrics and weights.  It
holds at every subdivision level because the conforming carrier remains a
closed `S^3` while every summand of the raw broken complex is a contractible
tetrahedron.

## What the no-go means

The failed condition was

\[
 LD_{\rm pw}=dL.
\]

It is too strong because the raw direct-sum complex has forgotten precisely
the global gluing information that creates the top cohomology class.  No
linear recovery can both reconstruct that information as a left inverse and
intertwine the unaugmented differential.

Broken FEEC avoids the contradiction by changing the operative differential:

\[
 d_h=D_{\rm pw}P,
 \qquad P=JL.
\]

This projected differential has the correct conforming kernel after
stabilization, as already certified, but its finite positive spectrum remains
projection-dependent.

The next mathematically coherent route is therefore not a larger-support
solution of the impossible equation.  It is an augmented local complex that
retains interface data explicitly—for example jump/trace or Cech-type
degrees of freedom—so its cohomology can contain the global class before any
recovery is applied.  Whether the existing geometry uniquely selects the
metric and first-order operator on such an augmentation is **OPEN**.

## Status ledger

- **DERIVED:** `J` is an exact signed cochain map.
- **DERIVED:** every strict-occurrence weight is forced to zero on every
  frozen level and degree `0,1,2`.
- **DERIVED:** this contradicts every left-inverse sum independently of
  positivity or candidate formula.
- **DERIVED:** both named local candidates have exactly nonzero cochain
  commutators.
- **DERIVED POST-RESULT:** the `H^3` mismatch forbids any raw-piecewise
  cochain retraction, even nonlocal.
- **STRUCTURAL CORRECTION:** increasing the recovery support cannot repair
  this equation; the broken complex itself must be augmented or its
  differential projected.
- **OPEN:** a canonical local interface-augmented Hilbert complex and its
  metric.
- **NOT CLAIMED:** failure of broken FEEC, all local discrete dynamics, time,
  causality, inertia, mass or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_commuting_recovery.py
```

Expected result: `8/8`.
