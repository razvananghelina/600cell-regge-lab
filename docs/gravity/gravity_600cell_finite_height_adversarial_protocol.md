# Frozen adversarial protocol: finite-height tangent classification

Date: 2026-08-21.

Primary protocol commit: `4b24abf`.

Primary verifier registration commit: `4176c3d`.

Primary artifact commit: `f0a4209`.

Primary artifact SHA-256:
`9bf4cc33d42d540e137f620eaf952d44ac49105648c828efba0ac8bdf4762f03`.

Status: frozen before the adversarial verifier exists or is executed.

## 1. Purpose and independence boundary

The primary proof eliminated the two affine residuals to

```text
E(v,q)=4*pi[mu(q)-mu(v)]+q[p(q)-p(v)]
```

and used the decisive identity `E_q=p(q)-p(v)`.  Repeating that proof would
only test reproducibility.

The adversarial implementation must instead reconstruct the finite residuals
from the complete action and classify their tangent-line error

```text
T(v,q)=p(v)-p(q)+(4*pi/q)[mu(v)-mu(q)],  q!=0.
```

It may compare `T=-E/q` only after deriving and classifying `T` independently.
Its decisive identity is the dual formula

```text
T_q(v,q)=4*pi[mu(q)-mu(v)]/q^2,
```

so the global partition is controlled by the two branches of `mu`, not by the
two branches of `p`.  The lost point `q=0` must be evaluated in the original
two residuals and may not be recovered by continuity of `T` alone.

## 2. Independent all-real root proof

Use only the already-adversarially-corroborated theorem that `K(q^2)` has one
positive zero `v_star`, together with the independently rederived signs

```text
mu'(q) has sign q*K(q^2),
mu is even,
p is odd and negative for q>0.
```

For `q>0`, split by the critical points where `mu(q)=mu(v)`:

- if `0<v<v_star`, these are `q=v` and one outer point;
- if `v_star<v<v_M`, these are one inner point and `q=v`;
- if `v=v_M`, the inner point is `q=0`;
- if `v>v_M`, only `q=v` remains.

Use `T_q`, the exact diagonal value `T(v,v)=0`, and

```text
lim_{q->infinity} T(v,q)=p(v)-p_infinity
```

to recover or refute the primary positive-axis root counts.  For `q<0`, use
the reflected `mu` critical points and the exact equal-mass values

```text
T(v,-v)=2*p(v)<0,
T(v,-m(v))=p(v)+p(m(v))<0
```

to exclude hidden negative-axis roots below `v_M`.  Above `v_M`, require one
root between `-v` and zero.  Prove `v<=0` by direct time reflection, treating
`v=0` separately.

## 3. Independent endpoint proof

On the negative-`q` branch derive from `T`, rather than from the primary `F`,

```text
T_q>0,
T_v=4*pi*mu'(v)(1/q-1/v)>0,
dq/dv=-T_v/T_q<0.
```

Thus `u(v)=-q(v)` increases.  Certify with real-ball arithmetic that the
root remains in `0<u<1<v_star`.  Since `mu(v)` strictly decreases and
`mu(u(v))` strictly increases,

```text
L_plus(v)=2*mu(v)/mu(u(v))-1
```

strictly decreases from `+1` to `-1` and has one zero.  This is the
adversarial uniqueness proof for the causal endpoint.

## 4. Numerical and convention attacks

Only after the exact tangent proof objects exist:

1. use rigorous Arb balls at rational points to bracket
   `v_A`, `v_star`, `v_M`, `v_C` and `q_C`;
2. solve the original two residuals directly in `(h,q)` at 80, 120 and 180
   decimal digits and require nested agreement;
3. recompute the primary thresholds without reading their decimal strings,
   then compare all values beyond 60 decimal digits;
4. verify the exact time-reflection convention
   `(v,q,h)->(-v,-q,-h)`;
5. use `v=1` as a known no-update control and a state beyond `v_C` as a known
   endpoint-positivity failure;
6. reverse the tangent sign and shift the conserved mass by `1/11`; both
   hostile constructions must fail the original residuals;
7. require nonzero direct Jacobians at every isolated representative.

## 5. Outcome boundary

If the dual `mu`-controlled proof, direct residual solve and convention
attacks all agree with the primary result, emit

```text
FINITE_HEIGHT_ISOLATED_UPDATES_WITH_CAUSALITY_BOUNDARY_ADVERSARIALLY_CORROBORATED.
```

Otherwise the result returns to `FINITE_HEIGHT_OPEN`, and the disagreement is
the headline.

Even a pass is only **DERIVED EXACT / STRUCTURAL** for this fixed homogeneous
action.  Published work already establishes dynamical frustum heights and
600-cell causal stopping points.  Composition, nonhomogeneous stability,
refinement and an absolute time scale remain **OPEN**.
