# Preregistration: local conserved-mass homothetic slab test

Date: 2026-08-16

Prior-art gate: `8865346`.

Status: **frozen before evaluating the non-static homothetic states in the
Regge action**.

## 1. Carrier, action and coordinates

Use both derived order-24 staircase schedule parities and exactly the
Lorentzian Regge plus fixed-total-mass dust action, logarithmic gradients,
complex-angle branch and orbit normalization of
`verify_gravity_600cell_dust_canonical_continuation.py`.

Fix

```text
L_- = L0,
rho = tau0^2,       tau0 = 0.0102,
M   = (90/pi)*(2*pi-5*acos(1/3))*L0.
```

Let `s=log(L_+/L0)`.  The geometrically homothetic state is

```text
q_old(s)       = L0^2,
q_new(s)       = exp(2s)*L0^2,
pole magnitude = rho,
diagonal(s)    = exp(s)*L0^2-rho.
```

No internal edge, mass, coefficient or momentum is optimized.  Only `s=0`
is evaluated before this protocol through the upstream static theorem.

The verifier must independently derive the diagonal formula from exact
unit-600-cell dot products and the five-dimensional Lorentzian embedding.
Every evaluated state must have 2400 simplex representatives counted with
exactly one timelike Gram direction, positive nonzero leading-minor margins,
angle-argument modulus above `1e-6`, and maximum action/gradient imaginary
contamination below `1e-70`.

## 2. Two equations that must remain separate

Let `g_D[30]` and `g_P[5]` be the complete per-edge-normalized logarithmic
internal gradients for the staircase diagonals and positive pole
magnitudes.  The complete local residual is

```text
R_local(s) = (g_D(s), g_P(s)) in C^35.
```

The globally restricted lapse derivative, divided by the common orbit size
24, is the chain-rule combination

```text
E_lapse(s)
  = sum_j g_D,j(s) * [-rho/diagonal(s)] + sum_k g_P,k(s).
```

This follows from varying `log rho` at fixed `s`, because

```text
d log(diagonal)/d log(rho) = -rho/diagonal,
d log(pole magnitude)/d log(rho) = 1.
```

Independently validate this chain rule at the displaced points using
central differences of the complete restricted action in `log rho`.  The
action derivative must agree with `24*E_lapse` inside the calibrated
finite-difference error.  A zero of `E_lapse` is not promoted to a solution
unless all 35 components of `R_local` vanish.

For this independent check, at every frozen nonzero `s` use

```text
k = 1e-6,
Q(k)   = [S(log rho+k)-S(log rho-k)]/(2k),
Q(k/2) = [S(log rho+k/2)-S(log rho-k/2)]/k.
```

Here every displaced `rho'` is inserted simultaneously into
`diagonal=exp(s)L0^2-rho'` and the five pole magnitudes.  Require

```text
|Q(k/2)-24*E_lapse(s)|
  <= 10*|Q(k/2)-Q(k)| + 1e-60.
```

The check is validation only; neither `Q` value enters the derivative
verdict of section 3.

## 3. Frozen local derivative calculation

At fixed `rho=tau0^2`, evaluate both parities at

```text
s = +/- h, +/- h/2, +/- h/4,
h = 1e-4.
```

For each scalar or vector function `f`, form

```text
D(h)   = [f(h)-f(-h)]/(2h),
D(h/2) = [f(h/2)-f(-h/2)]/h,
D(h/4) = [f(h/4)-f(-h/4)]/(h/2),
R4(h)  = [4D(h/2)-D(h)]/3,
R4(h/2)= [4D(h/4)-D(h/2)]/3.
```

Use `R4(h/2)` as the operational derivative and
`|R4(h/2)-R4(h)|` componentwise as the truncation proxy.  Add an arithmetic
floor of `1e-70`.  A scalar is resolved nonzero only if its absolute value
exceeds 100 times its proxy.  A vector is resolved nonzero only if its
Euclidean norm exceeds 100 times the Euclidean proxy norm and at least one
component is individually resolved by the same rule.

Record:

- `dR_local/ds` with all 35 components;
- `dE_lapse/ds`;
- derivatives of the 30 pre- and 30 post-momenta as diagnostics;
- within-type spreads and even/odd differences;
- all branch margins at every displaced point.

No derivative sign or expected magnitude is preregistered.  No tolerance
may be relaxed after inspection.

## 4. Static and independent controls

Before using a derivative verdict, require:

1. the upstream exact artifact reports `REGULAR_LAPSE_IDENTITY_PROVED` and
   `13/13`;
2. a fresh 100-decimal evaluation at `s=0` has all 35 local residuals below
   `1e-60` and reproduces the exact pre/post momentum formula below
   `1e-60` relative error;
3. the exact identities
   `2(1-phi/2)=phi^-2` and `diagonal=L_-L_+-rho` pass symbolically;
4. the action-gradient chain rule in section 2 passes at all twelve
   displaced parity/state combinations;
5. the operational derivatives from both parities agree within ten times
   their combined proxy norms before a parity-independent physical reading
   is allowed.

Pre-evaluation completeness clarification: item 5 applies to
`dR_local/ds` and `dE_lapse/ds`, which are the two decision objects of this
mission.  The pre/post-momentum derivatives are reported diagnostics and do
not enter this schedule gate; comparing their components physically would
require a separately reconstructed boundary-orbit identification.

If the two schedules disagree beyond that gate, the primary outcome is
schedule-dependent and no homogeneous physical conclusion is allowed.

## 5. Mechanical outcome hierarchy

Assign exactly one primary outcome, in this order:

1. `HOMOTHETIC_CONTROL_FAILED` if a static, exact-geometry, chain-rule or
   Lorentzian-branch control fails;
2. `HOMOTHETIC_SCHEDULE_DEPENDENT` if the even/odd derivative gate fails;
3. `LOCAL_HOMOTHETIC_STATIC_ONLY_GLOBAL_AND_LOCAL` if both
   `dE_lapse/ds` and `dR_local/ds` are resolved nonzero;
4. `LOCAL_HOMOTHETIC_STATIC_ONLY_LOCAL` if `dR_local/ds` is resolved
   nonzero but `dE_lapse/ds` is not;
5. `LOCAL_HOMOTHETIC_GLOBAL_ONLY_LOCAL_TANGENT` if `dE_lapse/ds` is
   resolved nonzero but `dR_local/ds` is not;
6. `LOCAL_HOMOTHETIC_TANGENT_SURVIVES` if neither derivative is resolved
   nonzero;
7. `HOMOTHETIC_NUMERICALLY_OPEN` if precision/proxy bands cannot decide.

For outcomes 3 or 4, the exact upstream line `s=0` together with one
resolved transverse local derivative gives a **DERIVED COMPUTATIONAL LOCAL**
implicit-function obstruction: in a neighborhood of the published point,
the complete homothetic internal system has only the static family.  This is
not a global no-go for disconnected roots or nonhomothetic internal lengths.

For outcome 3, the same statement also applies to the global restricted
lapse equation.  It means that fixed mass locally forces `L_+=L_-`; it does
not select a nonzero tick.

For outcome 5, the global equation alone cannot establish a solution of the
complete carrier and is labelled **STRUCTURAL / minisuperspace-only**.  For
outcome 6, no tick is yet accepted; it only authorizes a separately
preregistered higher-order/branch continuation.

## 6. Acceptance and kill boundaries

This mission accepts no physical evolution from a nonzero derivative alone.
A physical tick would still require a non-static state that passes all 35
local equations and the independent canonical junction momentum.  The
present local test can kill the geometrically rigid homothetic route near the
published sandwich, but cannot kill the full 65-variable nonregular route.

If a global lapse equation survives while local equations fail, label it
**STRUCTURAL / minisuperspace-only**.  If a lapse value appears only because
the finite curved triangulation breaks gauge, label it **STRUCTURAL /
candidate pseudo-constraint**, not an emergent clock.

Only the targeted verifier for this mission will be run.  It must be
registered exactly once in `reproducible/run_all.py`; the full suite will not
be run.
