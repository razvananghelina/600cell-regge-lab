# Preregistration: exact regular-lapse identity of the 600-cell dust slab

Date: 2026-08-16

Prior-art gate: `898d123`

Status: **frozen before any new evaluation of the candidate identity away
from the already committed continuation data**.

## 1. Carrier, action and domain

Use both derived order-24 schedule parities and exactly the complete
Lorentzian Regge plus De Felice--Fabri dust action, angle branch, edge signs,
orbit multiplicities and per-edge logarithmic momentum convention of
`verify_gravity_600cell_dust_canonical_continuation.py`.

Set `L=ARB_L0`, `rho0=ARB_TAU^2`, and define

```text
q_old = q_new = L^2,
x_diagonal = L^2-rho,
x_pole magnitude = rho,
physical pole square = -rho,
tau = sqrt(rho)>0.
```

The proof domain is the connected positive-lapse interval

```text
0 < rho <= rho0.
```

No assertion is made at `rho=0`, where the slab is degenerate, or for a
negative/orientation-reversed lapse.

## 2. Candidate exact formula

Let

```text
epsilon3 = 2*pi-5*acos(1/3),
M        = (90/pi)*epsilon3*L.
```

The static-prism factorization predicts, in the conventions of the code,

```text
S_grav(rho) = 720*epsilon3*L*tau,
S_dust(rho) = -8*pi*M*tau = -720*epsilon3*L*tau,
S_total(rho)=0,

g_internal(rho)=0 on all 35 internal orbit coordinates,
p_pre,e(rho) = -(epsilon3*L/4)*tau for every old boundary edge e,
(P p_post)_e(rho) = +(epsilon3*L/4)*tau.
```

The boundary formula is for the stored per-edge derivative with respect to
the logarithm of the squared edge length.  Its sign must be checked against
the implementation; no absolute-value rescue is allowed.

The action restriction alone does not imply the 35 unrestricted internal
equations.  Those equations and the boundary derivatives must be proved
separately.

## 3. Exact symbolic route

Normalize `L=1` and write `u=rho/L^2`.  Reconstruct the slab combinatorially
for each parity and deduplicate all regular simplex, hinge and edge-incidence
types without using floating-point geometry.

For every type:

1. construct its squared-distance Gram matrix over `Q(u)`;
2. compute signed simplex, facet and hinge volume squares exactly;
3. derive the corrected Lorentzian angle argument used by the existing code;
4. reduce every Schlaefli gradient component by exact symbolic algebra;
5. combine logarithms only after certifying their branch continuously from
   the committed base point throughout `0<u<=u0`.

The branch certificate must prove one timelike Gram direction, nonzero
leading minors, nonzero angle arguments and no logarithm-cut crossing on the
open interval.  Exact polynomial signs may be certified by factorization,
Sturm/root isolation or equivalent rational algebra.  Sampling is not a
branch proof.

An exact proof must establish all 35 internal components and all 30 pre/post
components in both parity carriers.  Symmetry may reduce identical
expressions, but the verifier must expand the result back to every stored
component and check its incidence multiplicity.

## 4. Independent numerical controls

Only after the symbolic expressions are fixed, evaluate the original
arbitrary-precision geometric implementation at the preregistered ratios

```text
rho/rho0 in {1, 3/4, 1/2, 1/4, 1/16, 1/256}.
```

Use 100 decimal digits.  At each point require:

- the same Lorentzian branch as the published datum;
- maximum imaginary contamination below `1e-70`;
- every internal gradient below `1e-60` absolutely;
- relative error below `1e-60` for each nonzero pre/post momentum formula;
- action-factorization error below `1e-60` relative to
  `max(1,|S_grav|,|S_dust|)`.

These controls may falsify a symbolic derivation or expose a convention
error.  They cannot promote a finite sample to a theorem.

## 5. Consequence for the committed homotopy

Independently establish from the exact time-reflection/orbit map that

```text
P p_post(rho0) = -p_pre(rho0).
```

Only if the candidate formula and this endpoint relation are both proved may
the continuation target be simplified to

```text
p(lambda)=(1-2 lambda)*p_pre(rho0).
```

For the positive-lapse family this would imply, without fitting,

```text
rho(lambda)=rho0*(1-2 lambda)^2,  0<=lambda<1/2.
```

The limiting point `lambda=1/2` is classified as a zero-lapse boundary, not
as an accepted Lorentzian root.  For `lambda>1/2`, a solution would require a
separately defined temporal-orientation branch or a different spatial
geometry; neither may be inferred from the positive-lapse identity.

## 6. Mechanical outcomes

Assign exactly one primary outcome:

1. `REGULAR_LAPSE_IDENTITY_PROVED` only if the exact symbolic route, branch
   certificate, complete-component expansion and every numerical control
   pass;
2. `REGULAR_LAPSE_IDENTITY_REFUTED` if any exact or numerical component is
   demonstrably nonzero or has the wrong coefficient/sign;
3. `REGULAR_LAPSE_PATTERN_ONLY` if all samples agree but exact symbolic or
   branch certification is incomplete;
4. `REGULAR_LAPSE_IDENTITY_NUMERICALLY_OPEN` if the exact and numerical
   routes cannot be evaluated reliably enough to decide.

If outcome 1 holds, additionally report whether the same-orientation
connected canonical homotopy is proved to meet the zero-lapse boundary at
`lambda=1/2`.  This consequence is forbidden for outcomes 2--4.

## 7. Acceptance and kill boundaries

Acceptance is an exact identity, not a fit to the 41 known continuation
points.  The coefficient `epsilon3*L/4`, the six control ratios and all
tolerances above are fixed before new evaluation.

A refutation kills the proposed analytic explanation but not the previously
certified local continuation.  A proof closes the same-orientation regular
branch: it would show that this branch changes only lapse and cannot itself
produce a nondegenerate forward spatial frame.  It would not establish a
physical clock, expansion, Friedmann dynamics or emergent time.

Only the targeted verifier for this mission will be run.  It must be
registered in `reproducible/run_all.py`; the full suite will not be run.

