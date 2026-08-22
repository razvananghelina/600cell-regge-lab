# Local stability of the `v=3/2` branch signature

Date: 2026-08-22.

Status: **DERIVED / ADVERSARIALLY CORROBORATED, LOCAL AND STRUCTURAL**.

## Complete hypotheses

The theorem assumes only:

- the fixed homogeneous tetrahedral-frustum 600-cell action;
- zero cosmological constant;
- conserved global dust;
- the committed pre/post canonical-momentum convention;
- positive proper height and positive endpoint scale;
- initial data on the analytic one-parameter curve

  ```text
  (m,pi)=(mu(v),p(v));
  ```

- the accepted invariant theorem on

  ```text
  D={(m,x): 0<m<=2/5, x=m*q>=125}.
  ```

No explicit neighbourhood radius, empirical target, continuum Friedmann
parameter or nonhomogeneous degree of freedom is an input.

## Theorem

There exists an unspecified real `epsilon>0` such that every incoming state
with

```text
|v-3/2|<epsilon
```

has the same ordered physical tree through slab four:

```text
slab 1: exactly one nontrivial physical child;
slab 2: exactly two physical children A<B in increasing q;
branch A: no physical child at slab 3 (DEAD);
branch B: exactly one slab-3 and one slab-4 child, then ENTERED_D.
```

Once branch B enters `D`, the already accepted invariant theorem gives exactly
one physical successor at every later finite step.

Thus the accepted `DEAD+ENTERED_D` history is not confined to one isolated
incoming point on the special curve.

## Why local constancy follows

The exact scalar relation obeys

```text
E_q=p(q)-pi.
```

The initial diagonal `q=v` is an exact persistent zero-height tangency, and

```text
p'(3/2)<0
```

proves it is nondegenerate. Every other root in the five visited states is
simple. All stationary values, origin values, analytic-tail coefficients,
physical gates and outgoing denominators are strictly separated from zero.
The branch-B entry has certified margins

```text
2/5-m > 0.0042556251475,
m*q-125 > 0.3317932609404.
```

The real-analytic implicit-function theorem continues every simple root and
recursive state locally. The nondegenerate diagonal tangency persists
identically; the stationary, origin and tail signs prevent an unseen real
root from appearing or escaping from infinity; the strict gate signs preserve
the physical and terminal labels. These are precisely the hypotheses needed
for some nonzero open neighbourhood.

No maximal or physically meaningful `epsilon` is inferred.

## Primary certificate

The interval-Newton primary at corrected verifier commit `676eb5d` returned

```text
RESULT: 10/10 checks passed
OUTCOME: LOCAL_SIGNATURE_PRIMARY_CERTIFIED
```

Artifact SHA-256:

```text
9f524cc22df8cfb5083f372481b3efd19868252b85551d56378327eea7a6d613.
```

It used outward-rounded Arb balls, a complete stationary/tail partition and
one strict interval-Newton inclusion for each non-diagonal root. All five
physical edges also satisfied a separately redifferentiated complete action.

The first two executions exposed a SymPy radical-positivity omission and two
Arb textual-export issues. They are preserved in
`gravity_600cell_finite_height_local_signature_first_run_failure.md`; the
scientific criteria were not changed.

## Mechanically different adversarial chain

The first direct rational-bracket route returned a preserved `4/11 OPEN`:

```text
artifact SHA-256
139dcee2e9ee021c131aae1090433fe16bd70c9f2b10ec52d32b0c5ebd7748a7.
```

It had strict endpoint signs but raw wide-interval derivative balls contained
zero. This was interval dependency, not a contradictory root.

The separately frozen monotone-factor resolution used the already accepted
exact formula

```text
p'(q)=-720*K(q^2)/[
  (q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8)
]
```

and the unique zero of `K` at `v_star`. With no interval Newton, subdivision
or discovery seed, endpoint signs plus the three exact monotonicity regions of
`p` certified every fixed integer bracket and reproduced all five primary
states, roots, gate signs and terminals. It returned a preserved `8/10 OPEN`
only because the last numerical enclosure of `r-(1+h*q)` was wider than an
auxiliary preregistered threshold:

```text
artifact SHA-256
70448c78be2156ef84fbaa986c543c6063bcca8ca4395ee77bdbf657ab2760d1.
```

All five identity balls contained zero. The final exact resolver proved

```text
r-(1+h*q)=-E/(2*pi*mu(q)),
```

verified `mu(q)>0` on every physical root interval, and rejected the false
relation `r=1+2*h*q` on all five edges. It returned

```text
RESULT: 7/7 checks passed
OUTCOME: LOCAL_SIGNATURE_ENDPOINT_IDENTITY_EXACTLY_RESOLVED
```

Artifact SHA-256:

```text
ccedea4620f7cd485381f8002a8fa29b39a7842a94867e430a630024e6e7eb60.
```

This closes the only remaining adversarial failure without loosening its
post-result width threshold.

## Provenance ledger

```text
5bace00  prior-art gate
6512791  primary protocol
241a7d4  primary verifier registered before first run
0a276b4  first primary failure preserved
a9c68dd  explicit radical-positivity resolution
676eb5d  final Arb export correction before successful run
c5a202b  primary artifact
da65e6b  first adversarial protocol
a5e0c1d  pre-implementation width consistency correction
bc30c63  first adversarial verifier registered
3ee8dcb  first adversarial OPEN preserved
55c2ef0  monotone-factor resolution protocol
1aa915f  monotone-factor verifier registered
91bb6db  8/10 identity-width OPEN preserved
f2bf410  exact identity-resolution protocol
78fcb84  exact identity resolver registered
```

## Evidential ledger and framing attack

- **DERIVED / ADVERSARIALLY CORROBORATED:** existence of some open
  neighbourhood on the special incoming curve with the exact frozen tree.
- **STRUCTURAL:** complete forward extendibility selects branch B on that
  neighbourhood within this homogeneous model.
- **OPEN:** the size of the neighbourhood, the complete incoming basin and the
  full two-dimensional canonical state plane.
- **STRUCTURAL WARNING:** future extendibility remains global in discrete
  time. This theorem does not make it a local equation of motion.
- **NOT DERIVED:** that `v=3/2` is selected, that the neighbourhood is large or
  generic, or that the post-hoc half-strip `D` is sharp or physical.
- **NOT TESTED HERE:** nonhomogeneous modes, local gravity, an absolute tick,
  `c`, `G`, Planck scales or particles.

The complete-domain discovery remains essential context: it found 36
depth-four signatures and 1080 diagnostic inputs with a branch still live
outside `D`. The local theorem removes the literal single-point objection but
does not turn the complicated global diagram into one generic branch law.

## Next decision

Do not compute more isolated homogeneous slabs. A full 50-cell interval atlas
would rigorously classify only depth four and still not decide the 1080 live
branches. The next physically informative gate is a geometry-selected
nonhomogeneous perturbation of this locally stable background, reconciled with
the already recorded no-go results for previous anisotropic carriers. No new
carrier may be fitted after seeing a desired spectrum.

Only the targeted local verifiers and documentation guard were run. No
full-suite result is claimed.

