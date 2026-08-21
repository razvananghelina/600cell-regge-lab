# Result: standard canonical and causal conditions do not select a branch

Date: 2026-08-21.

## Provenance

```text
post-result prior-art and selector protocol       e3c8927
primary verifier registered before first run      70129cf
primary selector artifact                         b0c2fe9
adversarial protocol                              94ee376
adversarial verifier registered before first run  9bec051
adversarial selector artifact                     eb0eae5
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_finite_height_selector_audit.py
  10/10 PASS

reproducible/verify_gravity_600cell_finite_height_selector_audit_adversarial.py
  10/10 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_finite_height_selector_audit.json
SHA-256 956cd655b8b3a5106029fb852df74b85bb59f922a4984542bc2e089f54799676

reproducible/gravity_600cell_finite_height_selector_audit_adversarial.json
SHA-256 1fe11f006cd928dc5418c3171154a4b0e26db79225e69845f2edd9566f820f0e
```

Both targeted verifiers were rerun together and reproduced these hashes byte
for byte.  No full-suite run was performed.  Static registry audit: `409/409`
distinct registrations, two deliberate exclusions, zero duplicates, zero
unregistered verifiers and zero missing files.

## Headline

```text
STANDARD_CANONICAL_SELECTORS_DO_NOT_RESOLVE_BRANCH_
ADVERSARIALLY_CORROBORATED
```

> **DERIVED NEGATIVE, SELECTOR-SCOPED / ADVERSARIALLY CORROBORATED:** the two
> physical second slabs at the frozen state `v=3/2` are both future oriented,
> lie on the same connected real causal branch of the cellular action and are
> locally regular sheets of the constrained discrete Legendre relation.
> Neither causality, time orientation, action-branch membership nor local
> regularity selects one of them.

This does not prove that no selector exists.  It proves that the standard
conditions above, applied without a fitted threshold, leave exactly the same
two branches that the composition calculation found.

## Complete hypotheses

Use the fixed homogeneous tetrahedral-frustum 600-cell action at zero
cosmological constant, with conserved global dust, positive proper strut
square `rho=h^2`, positive endpoint scale and the committed canonical
pre/post momentum convention.  Fix the accepted first state `v=3/2` and its
exact outgoing normalized canonical data `(m1,pi1)`.

The result does not cover an added matter clock, a boundary wave function or
path-integral contour, nonhomogeneous edge data, carrier refinement, a
cosmological constant or an independently derived global extendibility
condition.

## Exact local-regularity theorem

For arbitrary normalized incoming data `(m,pi)`, write

```text
C(h,q)=8*pi[mu(q)-m]+4*pi*h*q*mu(q),
P(h,q)=p(q)-pi-2*pi*h*mu(q).
```

Direct expansion gives

```text
det partial(C,P)/partial(h,q)
  -8*pi^2*h*mu(q)^2
  =4*pi*mu(q)[4*pi*mu'(q)+q*p'(q)].
```

The independently established state-function identity

```text
4*pi*mu'(q)+q*p'(q)=0
```

therefore implies

```text
det partial(C,P)/partial(h,q)=8*pi^2*h*mu(q)^2.
```

Moreover, `mu(q)>0` for every real `q`: the inverse-cosine argument lies in
`[1/3,1/2)`, and

```text
1/3 > cos(2*pi/5)=(sqrt(5)-1)/4.
```

Hence the determinant is strictly positive at every positive-height root.
Both physical solutions are locally invertible branches.  Since the same
incoming `(m1,pi1)` has two such outputs, the constrained discrete Legendre
relation is not globally injective at this state.  Local regularity cannot be
promoted to global uniqueness.

## Exact causal correction

The present `h` is the proper length of a timelike same-vertex strut.  It is
not the central-coordinate height used in the superficially similar cubical
frustum inequality.  In the certified tetrahedral Minkowski embedding,

```text
R=phi*L,
Delta_R=phi*h*q,
T^2=h^2+Delta_R^2,
strut^2=-T^2+Delta_R^2=-h^2<0.
```

The coordinate speed is therefore

```text
beta(q)^2=Delta_R^2/T^2
         =phi^2*q^2/(1+phi^2*q^2)<1
```

for every finite real `q`.  The exact causal margin is

```text
1-beta(q)^2=1/(1+phi^2*q^2)>0.
```

The two frozen branches give

| Branch | `q` | `h` | `beta^2` | causal margin |
|---|---:|---:|---:|---:|
| A | `0.0212594...` | `7.2893013...` | `0.00118186...` | `0.998818...` |
| B | `31.2792236...` | `0.0689338...` | `0.99960975...` | `0.000390250...` |

Branch B is near-null in central coordinates but remains strictly timelike.
There is no exact causal boundary between A and B.  Falsely setting the
central time to `T=h` makes branch B appear to have speed
`phi*q=50.6108...`; the adversarial verifier requires this convention trap to
fail.

Both roots also satisfy

```text
1/3 <= (q^2+2)/(2(q^2+3)) < 1/2
```

and have finite real boost argument, so they lie on the same real analytic
cellular-angle branch.

## Independent replication

The primary verifier used the reduced exact equations and the closed causal
formula.  The adversarial verifier instead used:

- roots obtained previously by direct two-variable solves of the complete
  differentiated action at 80, 120 and 180 decimal digits;
- direct full-action Jacobians at all three precisions;
- explicit Minkowski coordinates `Delta_R` and `T` rather than the primary
  `beta` identity;
- hostile central-height and Euclidean-sign conventions.

Only after these classifications were complete did it read the primary
selector artifact.  The roots, Jacobians and squared speeds then agreed
beyond 55 decimal digits.

## Status ledger

| Claim | Status |
|---|---|
| Every positive-height root of the reduced canonical equations is locally regular | **DERIVED EXACT** |
| The two `v=3/2` continuations are both future oriented | **DERIVED COMPUTATIONAL** |
| Both continuations have timelike proper struts and `beta^2<1` | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| Both lie on the same real action branch | **DERIVED / ADVERSARIALLY CORROBORATED** |
| Local Legendre regularity selects one continuation | **DERIVED NEGATIVE** |
| Causality or time orientation selects one continuation | **DERIVED NEGATIVE, selector-scoped** |
| The action is globally injective at the frozen state | **DERIVED NEGATIVE** |
| No possible additional physical selector exists | **OPEN** |
| A global extendibility condition selects a unique infinite history | **OPEN** |
| Nonhomogeneous equations remove the branching | **OPEN** |
| Refinement removes the branching | **OPEN** |
| A deterministic or fundamental tick is derived | **NO** |

## Interpretation and next gate

The homogeneous action defines a locally regular, globally multivalued
canonical relation.  This is mathematically coherent but is not yet a
deterministic evolution law.  Selecting the slow or fast branch because of
its numerical appearance would be fitting.

The next target-free test is **future extendibility**: construct the complete
third-slab root set for each of the two branches before judging either one.
If exactly one branch extends while the other terminates, that is a structural
global distinction, although promoting it to a physical selection principle
would still require justification.  If both extend, the branching is genuine
at a longer horizon.  If neither extends, the finite update is a short-lived
boundary-value relation rather than an evolution.
