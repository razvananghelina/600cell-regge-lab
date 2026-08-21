# Frozen classification protocol: finite-height tangent branch

Date: 2026-08-21

Original finite protocol commit: `79747ea`.

Unexecuted verifier registration commit: `483a733`.

First exact construction commit: `4895e5f`.

First artifact SHA-256:
`c386d5dc16630ac4915f3ff634a0eb53e28b5d9f9760cdfaba225fb81fa47d4b`.

Status: frozen after the exact affine construction, before any numerical root
scan, threshold value or finite-update classification.

## 1. Exact first-run reductions

The construction run derived

```text
C0(v,q)=8*pi[mu(q)-mu(v)],
P0(v,q)=p(q)-p(v),
Ch(q)=4*pi*q*mu(q),
Ph(q)=-2*pi*mu(q).
```

Since `epsilon(q)>0`, one has `mu(q)>0` for every finite real `q`.  Therefore

```text
D(v,q)=-4*pi*mu(q) E(v,q),

E(v,q)=4*pi[mu(q)-mu(v)]+q[p(q)-p(v)].
```

The verifier correction may certify these identities exactly but may not
change the action, state family, physical inequalities or determinant sign.

For `q!=0`, every non-boundary root reconstructs

```text
h=[p(q)-p(v)]/[2*pi*mu(q)],
L_plus=1+hq=2*mu(v)/mu(q)-1.
```

At `q=0`, the lapse equation requires `mu(v)=mu(0)` and the momentum equation
gives `h=-p(v)/[2*pi*mu(0)]`.

## 2. Exact derivative structure

Let

```text
x=q^2,
r=sqrt(x+4),
s=sqrt(3x+8),
epsilon(x)=2*pi-5*acos((x+2)/(2(x+3))),
K(x)=10*r-(x+3)s*epsilon(x).
```

Certify exactly

```text
mu'(q)=180*q*K(q^2)/[
  pi*(q^2+4)^(3/2)*(q^2+3)*sqrt(3q^2+8)
],

p'(q)=-4*pi*mu'(q)/q
     =-720*K(q^2)/[
       (q^2+4)^(3/2)*(q^2+3)*sqrt(3q^2+8)
     ],

partial_q E(v,q)=p(q)-p(v).
```

Use the already-corroborated proof that `K` has exactly one positive squared
root `x_star`, with `v_star=sqrt(x_star)`.  Consequently, on `q>0`:

- `mu` increases on `(0,v_star)` and decreases on `(v_star,infinity)`;
- `p` decreases on `(0,v_star)` and increases on `(v_star,infinity)`.

Use oddness of `p` and evenness of `mu` for the negative half-line.

## 3. Frozen threshold definitions

Define, before numerical evaluation,

```text
p_infinity=60*pi-300*sqrt(3)*log(2)<0.
```

The inequality may use `log(2)>2/3`, `sqrt(3)>5/3` and `pi<22/7`.

Define three positive thresholds only by their equations and monotone
intervals:

```text
v_A in (0,v_star):
  p(v_A)=p_infinity,

v_M in (v_star,infinity):
  mu(v_M)=mu(0),

(v_C,q_C), with v_C in (v_M,infinity), q_C<0:
  E(v_C,q_C)=0,
  mu(q_C)=2*mu(v_C).
```

`v_A` is the state where a nontrivial root enters from `q=+infinity`.
`v_M` is the state whose root crosses `q=0`.  `(v_C,q_C)` is the candidate
endpoint-positivity or causality boundary `L_plus=0`.

No decimal threshold may be used to define a bracket.  Brackets must be
fixed by rational endpoints whose signs are printed, and all roots must be
proved unique before a decimal is reported.

## 4. Registered all-real root-count theorem under test

Using `E_q=p(q)-p(v)` and the monotone branches of `p`, prove or refute:

### Nontrivial roots of `E`

- `v<=0`: any off-diagonal roots have `h<0`; no positive-height update;
- `0<v<=v_A`: no finite off-diagonal root;
- `v_A<v<v_star`: exactly one off-diagonal root, with `q>v_star`;
- `v=v_star`: no off-diagonal root;
- `v_star<v<v_M`: exactly one off-diagonal root, with `0<q<v_star`;
- `v=v_M`: exactly one root at `q=0`;
- `v>v_M`: exactly one off-diagonal root, with `q<0`.

The ubiquitous diagonal root `q=v` reconstructs only `h=0` and is never a
positive update.

### Physical endpoint filter

On every nontrivial positive-height root require

```text
L_plus=2*mu(v)/mu(q)-1>0.
```

Prove or refute that the endpoint ratio crosses zero exactly once on the
`v>v_M` branch, at `(v_C,q_C)`, and that the complete positive physical
state set is

```text
v in (v_A,v_star) union (v_star,v_C).
```

At `v_A` the root is at infinity; at `v_star` it merges with the zero-height
diagonal; at `v_C` the endpoint is zero.  All three endpoints are excluded.

If uniqueness of `v_C` or any global root count cannot be proved, the outcome
remains `FINITE_HEIGHT_OPEN` regardless of numerical agreement.

## 5. Deterministic numerical diagnostics after the proof objects exist

Only after all exact identities and monotonic intervals are printed:

1. isolate `v_A`, `v_M` and `(v_C,q_C)` with rational sign brackets;
2. evaluate one representative state in each predicted interval;
3. bracket every off-diagonal `q` using the monotonic intervals specified
   above;
4. reconstruct `h` and `L_plus` from both original residuals at 100 decimals;
5. require residuals below `1e-80` and a nonzero Jacobian at every isolated
   representative.

Diagnostics may falsify the theorem but cannot prove its global quantifiers.

## 6. Interpretation frozen before threshold values

If the theorem passes, the outcome is

```text
FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY.
```

Label it **DERIVED EXACT / STRUCTURAL**, not a fundamental tick:

- the height is a state-dependent pseudo-constraint update `h(v)`, not one
  universal interval;
- the branch has a zero-height puncture at `v_star` and terminates at
  `L_plus=0`;
- composition, perturbative stability and refinement are untested;
- published 600-cell Regge cosmologies already contain causality-breaking
  stopping points, so only the exact present coefficient may be externally
  new, and that remains **OPEN**.

No result derives seconds, `c`, `G` or Planck time.
