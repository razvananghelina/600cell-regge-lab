# Resolution protocol: global Lagrange certificate after direct-quotient OPEN

Date: 2026-08-22.

Parent adversarial protocol commit: `26ef9c3`.

First adversarial implementation commit: `87d9aac`.

Preserved first result commit: `bfa4db4`.

Status: frozen after preserving the first adversarial result and before the
resolution verifier exists or is executed.

## 1. Preserved result and diagnosed limitation

The parent adversarial protocol returned

```text
INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN
```

with artifact SHA-256

```text
f7d1f36e5ed679c39d1c38dbc21509ae52211f6735b38a2da46046fb798f54d5.
```

At both 160 and 256 decimal digits, the Taylor-treated axis leaf passed but
all 63 leaves using interval evaluation of the direct quotients failed to
separate the signs of `Bbar` and `-C'`.  The enclosures were essentially
unchanged with precision.  The failure is dependency loss in expressions
whose first-order terms cancel, not roundoff and not a certified wrong sign.

The complete direct-action reconstruction of all five slabs and every hostile
control passed at both precisions.  The first result remains `OPEN`; this
protocol does not overwrite or reinterpret it.

The aggregate symbolic redifferentiation flag also failed because a one-shot
SymPy simplification did not reduce the exact identity

```text
p'(q)=-4*pi*mu'(q)/q.
```

The action and post-momentum scaling identities simplified exactly, the
derived `R'` identity simplified exactly, and the direct action passed.  The
resolution verifier must replace only this incomplete simplifier call by the
explicit radical chain-rule factorization already stated below.

## 2. Unchanged theorem, domain and gates

Test the same post-hoc half-strip

```text
0<m<=2/5,
x>=125
```

with the same homogeneous 600-cell action, zero cosmological constant,
conserved global dust, momentum convention, physical inequalities, global
root argument, two precisions, fifth-slab control and outcome hierarchy.

No threshold, sign, branch, seed, target or accepted input may change.  Do
not use `x>=124` or `x>=126` for acceptance.  Do not subdivide the domain.

## 3. Mechanically distinct resolution

Keep

```text
u=m^2/x^2,
0<=u<=4/390625,
N=P-P0+2*pi*u*M,
H=u*P'-(P-P0).
```

The primary certificate used exact integral means of derivatives.  The first
adversarial attempt used direct quotient interval evaluation on 63 leaves.
The resolution must use neither method.  Use one origin-centred Maclaurin
polynomial with a global Lagrange derivative remainder over the complete
rational `u` interval.

At 160 and 256 decimal digits, obtain outward-rounded Arb coefficients at
`u=0` and derivative ranges on the full interval.  Through degree six use

```text
W(u)=sum_(k=1)^6 p_k*u^(k-1)
     + [P^(7)(I)/7!]*u^6,

Bbar(u)=sum_(k=2)^6 n_k*u^(k-2)
        + [N^(7)(I)/7!]*u^5,

W'(u)=sum_(k=2)^6 h_k*u^(k-2)
      + [H^(7)(I)/7!]*u^5,
```

where

```text
p_k=P^(k)(0)/k!,
n_k=N^(k)(0)/k!,
h_k=H^(k)(0)/k!.
```

Require exact Arb overlap of `h_k` with `(k-1)*p_k` for `2<=k<=6` as an
internal identity control.  Then set

```text
C=W+4*pi*M,
-C'=-W'-4*pi*M'.
```

Evaluate `M` and `M'` directly on the one full interval; they contain no
removable quotient.  Use the same conservative rectangle bounds as the
primary proof for `U`, `y_plus/z`, `partial_z Y` and the normalized same-`x`
gap.  Every load-bearing ball must serialize explicit lower and upper
endpoints.

This route is an interval theorem on the continuum domain, not a grid or a
finite trajectory scan.

## 4. Exact derivative correction

Re-differentiate the full action.  For the state derivative, do not ask a
generic simplifier to discover nested radical identities.  With

```text
r=sqrt(q^2+4),
s=sqrt(3*q^2+8),
K=10*r-(q^2+3)*s*epsilon,
```

verify exactly

```text
mu'=180*q*K/[pi*r^3*(q^2+3)*s],
p'=-720*K/[r^3*(q^2+3)*s],
p'+4*pi*mu'/q=0.
```

The radical identities used in this reduction must themselves simplify to
zero.  This is an exact implementation repair, not a changed physical claim.

## 5. Controls and comparison order

Require again:

- the naively divided zero-containing direct quotient is unresolved;
- the preserved 63-leaf direct-quotient artifact is `OPEN` and has the frozen
  hash above;
- reversing the outgoing momentum sign fails;
- omitting the boost changes `Bbar(0)`;
- `x=60` fails the strict positive-height gate;
- the fifth successor is reproduced by the complete action at both
  precisions;
- the two rigorous precision records overlap with the same signs.

Build the complete resolution verdict before reading the primary invariant
artifact.  Compare with the primary only afterward.

## 6. Outcomes

### `INVARIANT_HALF_STRIP_ADVERSARIALLY_CORROBORATED`

Both global Lagrange certificates, the exact derivative factorization, direct
action, hostile controls and post-construction primary comparison pass.  Then
the primary theorem is adversarially accepted despite the separately
preserved direct-quotient `OPEN` route.

### `INVARIANT_HALF_STRIP_ADVERSARIAL_REFUTED`

A strict wrong sign or direct-action disagreement is certified.  Preserve the
witness and retract the primary theorem.

### `INVARIANT_HALF_STRIP_ADVERSARIAL_OPEN`

Any coefficient, remainder, sign, provenance or comparison remains
unresolved.  Do not increase the degree, subdivide, shrink the domain or
change a threshold in this mission.

Even a positive outcome remains homogeneous, dimensionless and scoped to the
representative accepted seed.  Infinite extendibility remains a STRUCTURAL
global admissibility criterion, not a local dynamical law.  Generic incoming
states, nonhomogeneous modes and an absolute tick remain OPEN.
