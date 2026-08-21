# Result: exact finite state-dependent updates with a causal endpoint

Date: 2026-08-21.

## Provenance

```text
initial prior-art gate                              c8edca2
all-real domain correction                         726e52f
finite-height construction protocol                79747ea
construction verifier registered before execution 483a733
exact affine construction                          4895e5f
classification frozen before threshold scan        4b24abf
primary classifier registered before execution     4176c3d
primary artifact                                   f0a4209
post-result literature audit                       820697b
adversarial tangent protocol                       ccc7a4e
adversarial verifier registered before execution   6934610
adversarial artifact                               7fa861c
```

Targeted verifiers:

```text
reproducible/verify_gravity_600cell_finite_height_classification.py
  12/12 PASS

reproducible/verify_gravity_600cell_finite_height_classification_adversarial.py
  11/11 PASS
```

Accepted artifacts:

```text
reproducible/gravity_600cell_finite_height_classification.json
SHA-256 9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03

reproducible/gravity_600cell_finite_height_classification_adversarial.json
SHA-256 da8d60e95b5196beaf93ea234fbf9dfb93e3d5e6bd00fb0a85ed2ef4ba388996
```

No full-suite run was performed.

## Headline

```text
FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY_
ADVERSARIALLY_CORROBORATED
```

> **DERIVED EXACT / STRUCTURAL / ADVERSARIALLY CORROBORATED:** the fixed
> homogeneous cellular 600-cell Regge-plus-conserved-dust action has exactly
> one nontrivial positive-height update on two explicitly classified state
> intervals.  The height is isolated by the two canonical equations.  The
> physical branch ends at a unique point where the next spatial scale is
> zero.

This is a state-dependent discrete pseudo-constraint update.  It is not a
universal tick and does not derive an absolute unit of time.

## Complete hypotheses

Use the certified positive-Lorentzian homogeneous cellular slab at zero
cosmological constant, conserved global dust and

```text
L_minus=1,
M=mu(v),
p_pre=p(v),
rho=h^2,
h>0,
L_plus=1+h*q>0,
v real.
```

No value of `v`, including `v=0` and `K(v^2)=0`, is removed from the finite
census.  The result does not cover a changed mass, changed incoming momentum,
nonhomogeneous edge data, another Regge action, a refined carrier or a
cosmological constant.

## Exact reduction

Define

```text
epsilon(t)=2*pi-5*acos((t^2+2)/(2(t^2+3))),

mu(t)=180*epsilon(t)/(pi*sqrt(t^2+4)),

p(t)=180*t*epsilon(t)/sqrt(t^2+4)
     -600*sqrt(3)*asinh(t/sqrt(8(t^2+3))).
```

After differentiating the complete action, the constraint and incoming
momentum residuals are exactly affine in `h`:

```text
C(v,q,h)=8*pi[mu(q)-mu(v)]+4*pi*h*q*mu(q),

P(v,q,h)=p(q)-p(v)-2*pi*h*mu(q).
```

Since `mu(q)>0`, simultaneous nontrivial roots are governed by

```text
E(v,q)=4*pi[mu(q)-mu(v)]+q[p(q)-p(v)]=0.
```

The original unknowns then reconstruct uniquely as

```text
h=[p(q)-p(v)]/[2*pi*mu(q)],
L_plus=1+h*q=2*mu(v)/mu(q)-1.
```

The diagonal `q=v` always gives only `h=0` and is not an update.

## Primary global proof

Let

```text
K(x)=10*sqrt(x+4)
     -(x+3)*sqrt(3*x+8)
       *[2*pi-5*acos((x+2)/(2(x+3)))].
```

The already-adversarially-corroborated theorem gives one positive zero
`x_star` of `K`.  Direct differentiation gives

```text
mu'(q)=180*q*K(q^2)/[
  pi*(q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8)],

p'(q)=-720*K(q^2)/[
  (q^2+4)^(3/2)*(q^2+3)*sqrt(3*q^2+8)],

E_q(v,q)=p(q)-p(v).
```

Thus `mu` rises and then falls on the positive axis, while `p` falls and then
rises.  Their common turning point is `v_star=sqrt(x_star)`.  Oddness of `p`,
evenness of `mu`, the exact negative tail

```text
p_infinity=60*pi-300*sqrt(3)*log(2)<0
```

and the signs of `E` at every stationary point give a complete all-real root
count.  No finite grid is used for the quantifier.

## Adversarial dual proof

The second route redifferentiated the complete action and did not reuse the
primary decisive identity.  For `q!=0` it classified the tangent-line error

```text
T(v,q)=p(v)-p(q)+(4*pi/q)[mu(v)-mu(q)].
```

Its independent decisive identity is

```text
T_q(v,q)=4*pi[mu(q)-mu(v)]/q^2.
```

The primary proof is therefore controlled by equal-`p` points, while the
adversarial proof is controlled by equal-`mu` points.  The point `q=0`, lost
by division, was checked directly in the original residuals.  The second
route also:

- used disjoint rational Arb sign brackets;
- repeated the threshold calculation at 80, 120 and 180 decimal digits;
- solved the two full-action equations directly in `(h,q)`;
- checked nonzero direct Jacobians;
- tested time reflection, a reversed tangent sign, a mass shift and a
  post-causal endpoint.

All six independently reconstructed threshold values agreed with the primary
artifact beyond 70 decimal digits.

## Complete root and physical classification

Define the thresholds intrinsically by

```text
p(v_A)=p_infinity,             0<v_A<v_star,
mu(v_M)=mu(0),                 v_M>v_star,
E(v_C,q_C)=0,
mu(q_C)=2*mu(v_C),             v_C>v_M, q_C<0.
```

Their diagnostic values are

```text
v_A    = 1.243256643819371910201962759595...
v_star = 2.376203781903790658629184442612...
v_M    = 16.01540335249024050451687117992...
v_C    = 31.46931050890825531254526180794...
q_C    = -0.1214673515316825288110275918484...
h_C    = 8.232664888055687601722677170121...
```

The exact nontrivial-root census is:

| Incoming state | Nontrivial root |
|---|---|
| `v<=0` | no positive-height update |
| `0<v<=v_A` | none |
| `v_A<v<v_star` | exactly one, `q>v_star` |
| `v=v_star` | none |
| `v_star<v<v_M` | exactly one, `0<q<v_star` |
| `v=v_M` | exactly one, `q=0` |
| `v>v_M` | exactly one, `q<0` |

Every nontrivial root for `v>0` has `h>0`.  Endpoint positivity reduces the
physical state set to

```text
v in (v_A,v_star) union (v_star,v_C).
```

For `v_A<v<v_M`, the update expands the scale (`L_plus>1`).  At `v_M` it has
`L_plus=1`.  For `v_M<v<v_C`, it contracts (`0<L_plus<1`).  At `v_C`,
`L_plus=0`; beyond it the algebraic root remains but violates the frozen
physical endpoint condition.

The causal endpoint is unique.  On the negative-`q` branch write `q=-u`.
The tangent proof gives `dq/dv<0`, hence `u(v)` increases while remaining
below `1<v_star`.  Therefore `mu(u(v))` strictly increases and `mu(v)`
strictly decreases, so

```text
L_plus(v)=2*mu(v)/mu(u(v))-1
```

strictly decreases from `+1` to `-1` and crosses zero once.

## What is and is not established

| Claim | Status |
|---|---|
| A finite positive update exists for the frozen homogeneous action | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| The update is isolated for a fixed admissible state | **DERIVED EXACT** |
| The complete all-real state domain is classified | **DERIVED EXACT** |
| The physical branch has one causal endpoint | **DERIVED EXACT / ADVERSARIALLY CORROBORATED** |
| The interval is universal or state independent | **DERIVED NEGATIVE** |
| This is a fundamental tick | **NOT DERIVED** |
| Consecutive updates compose into a unique evolution under the frozen action | **DERIVED NEGATIVE in the later representative-scoped composition test** |
| The branch is stable to nonhomogeneous perturbations | **OPEN** |
| The branch survives carrier/action refinement | **OPEN** |
| Seconds, `c`, `G` or Planck time follow | **NOT DERIVED** |
| External novelty of the exact coefficient theorem | **OPEN** |

The later composition result is
`docs/gravity/gravity_600cell_finite_height_composition_result.md`.  It does
not alter this one-slab classification; it refutes the stronger interpretation
of the classified update as a unique deterministic tick.

## Prior-art interpretation

[De Felice--Fabri, arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093)
and [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077) already evolve a
dust-filled 600-cell and identify a causality-breaking stopping point.
[Jercher--Steinhaus, arXiv:2312.11639](https://arxiv.org/html/2312.11639v2)
show in a different Lorentzian frustum cosmology that matter can make the
height dynamical and uniquely fixed subject to a causal inequality.

Therefore the physical mechanism is **KNOWN**.  The present contribution is
an exact internal root theorem for the frozen cellular action.  A dedicated
literature review would be required before any claim that its coefficient or
threshold theorem is externally new.

## Next gate

The next question is not whether another isolated one-slab root exists.  It
is whether the selected endpoint is admissible canonical data for a second
slab and whether the second selected update matches the first slab's outgoing
momentum without fitting.  That is the composition gate.  Only after it
passes is nonhomogeneous stability or refinement informative.
