# Preregistration: published 600-cell dust-sandwich control

Date: 2026-08-13

Prior-art gate commit: `e7d8bd5`

Upstream certified commits: `d9fe159`, `b11185e`

Status: **protocol frozen before evaluating the repository action at the
published dust-sandwich point**.

This is an external reproduction/control, not a search.  No parameter is
adjusted to reduce a residual, and no `a1=5`, speed, Planck scale or measured
cosmological quantity is compared.

## 1. Published source and frozen data

Use De Felice and Fabri,
[The Friedmann universe of dust by Regge Calculus](https://arxiv.org/abs/gr-qc/0009093),
sections 3.1--3.2.3.  Reconstruct their displayed initial sandwich from the
unrounded formulas, not by fitting their printed decimals:

```text
M_star = 10
zeta = (pi^2*sqrt(2)/50)^(1/3)
R0 = 4*M_star/(3*pi)
l0 = zeta*R0
epsilon_3 = 2*pi - 5*acos(1/3)
M = (90/pi)*epsilon_3*l0
tau = 0.0102
```

The two boundary 600-cells have equal squared edge length `l0^2`, every pole
has squared length `-tau^2`, and every staircase slant has

```text
d^2 = l0^2-tau^2.
```

Before action evaluation, reproduce the paper's printed controls

```text
l0^2             = 7.69379990138304  (printed)
l0^2-tau^2       = 7.69369586138304  (printed)
d^2 solved       = 7.69369586138301  (printed)
new l^2 solved   = 7.69379990138297  (printed)
M approximately  = 10.202...
```

only to the precision justified by those displayed values.

## 2. Frozen action and source normalization

Use the exact 65-variable Lorentzian curvature action and angle branch from
`verify_gravity_global_boundary_legendre.py`.  No cosmological/volume term is
added.

The paper's action is

```text
S_paper = (1/(8*pi))*S_gravity - (M/120)*sum_poles tau_i.
```

Work with the equivalent action multiplied by `8*pi`:

```text
S_total = S_gravity - (8*pi*M/120)*sum_poles sqrt(rho_i).
```

Each of the five pole orbits contains 24 poles, so its dust contribution and
derivative with respect to the common squared pole magnitude are frozen as

```text
S_dust,k = -(8*pi*M/5)*sqrt(rho_k),
dS_dust,k/drho_k = -(4*pi*M)/(5*sqrt(rho_k)).
```

No dust derivative acts on a staircase diagonal or boundary edge.  Confirm
the pole normalization independently by the paper's length equation

```text
2*tau*(1/24)*dS_gravity/drho_k = pi*M/15
```

when the total pole residual vanishes.  Do not change a sign after seeing the
result; a sign mismatch is a control failure to diagnose.

## 3. Frozen complete-action controls

For both ordered-schedule parities:

1. reconstruct all 30 diagonal, five pole and 30 final-boundary variables;
2. require all 2400 four-simplices to be nondegenerate Lorentzian and on the
   real certified angle branch;
3. compare the 100-simplex-orbit evaluator with the direct 2400-simplex sum in
   the gravitational action and all 65 gravitational derivatives;
4. add the analytic dust action and five derivatives;
5. reproduce all 35 total internal derivatives with centered differences of
   the complete total action at relative step `3e-6`.

Frozen implementation tolerances are:

```text
full/reduced action and gradient relative error     3e-8
maximum imaginary residual                         3e-7
direct total-action derivative relative error      5e-5
minimum absolute Gram eigenvalue                    > 1e-8
minimum angle-argument modulus                      > 1e-6.
```

These gates validate the evaluator.  A nonzero physical residual is a result,
not a verifier failure.

## 4. Frozen full-equation audit

Divide every orbit derivative by 24 to obtain its local edge equation.  Report
for each parity:

- all five total pole residuals;
- all 30 staircase-diagonal residuals;
- maximum absolute and Euclidean norms of each block;
- the full 35-component norm;
- for each of the ten unordered phase pairs, the three separate diagonal
  residuals and their sum;
- spread within each phase-pair triple;
- differences between the two phase parities.

The phase-pair sum is the derivative of the globally restricted action when
its three stabilizer-orbit lengths are forced equal.  It is **not** to be
silently identified with the paper's local-variation choice: the later
global-versus-local Regge literature distinguishes those procedures.

Use one frozen stationarity threshold:

```text
absolute per-edge residual <= 1e-7.
```

No averaging can turn a failing individual edge equation into a full pass.

## 5. Frozen outcome labels

Assign exactly one primary outcome for each parity:

- **FULL_REPRODUCTION:** all 35 individual residuals pass;
- **POLE_ONLY_REPRODUCTION:** all five sourced pole equations pass but at
  least one of the 30 diagonal equations fails;
- **SOURCE_NORMALIZATION_MISMATCH:** at least one sourced pole equation fails.

Also report, without promoting it to full reproduction, whether all ten
phase-pair sums pass.  If sums pass while individual triples fail, label this
**STRUCTURAL RESTRICTED-ACTION CANCELLATION**.

The overall control advances only if at least the pole/source normalization
is reproduced.  A full pass would validate the published sandwich directly
on the complete one-slab action.  A pole-only result would show that the
published Hamiltonian/source balance is present but its full evolution uses
additional lapse/shift, adjacent-slab or equation-selection structure not
captured by calling the one-slab equality point stationary.

## 6. Scope boundary

This protocol does not:

- reproduce the paper's entire multi-tick Newton evolution;
- assume that a one-slab boundary-value variation is identical to its
  sequential initial-value algorithm;
- test a cosmological constant;
- search any 35- or 65-dimensional neighborhood;
- select `M`, `tau`, the equality ansatz or a physical clock rate.

After the result and the mandatory post-result literature check, decide
whether the correct next carrier is one slab, two adjacent slabs with a shared
Cauchy surface, or the explicit five-stage canonical update.  Do not choose
that architecture merely because it makes the published numbers pass.
