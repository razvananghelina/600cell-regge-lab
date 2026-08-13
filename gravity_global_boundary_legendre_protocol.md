# Preregistration: global Regge boundary-Legendre rank

Date: 2026-08-13

Upstream certified commits: `d9fe159`, `0466ade`

Status: **protocol frozen before evaluating any 65-variable action gradient or
rectangular Jacobian**.

No root is searched here.  No physical target, `a1=5`, coupling, speed or
Planck scale is compared.

## 1. Corrected dynamical framing

The equal-boundary search fixed both spatial metrics and found no stationary
internal solution from its twelve frozen starts.  That is stronger than an
ordinary initial-value question and must not be treated as the definition of
dynamics.

For each ordered phase parity, keep the old 600-cell boundary regular and
allow the new boundary metric to vary.  The complete invariant coordinates
are:

```text
30 staircase-diagonal squared lengths q_j > 0
 5 pole squared magnitudes             rho_k > 0
30 final-boundary squared lengths      f_l > 0
------------------------------------------------
65 variables.
```

Only the 35 internal-edge derivatives are constraints.  The 30 derivatives
with respect to `f_l` are post-momenta.  They must be computed and real, not
set to zero.

## 2. Exploratory counts already seen

Before this commit, a scratch orbit enumeration found for each parity:

```text
720 old boundary edges   = 30 orbits x 24
720 final boundary edges = 30 orbits x 24
2400 boundary triangles  = 100 orbits x 24.
```

For every one of the ten unordered phase pairs, the 72 final edges again
split into three orbits of 24.  The verifier must reconstruct these counts.

No action, momentum or Jacobian using variable final edges was evaluated
before this protocol.

## 3. Frozen action construction

Extend the certified 35-variable evaluator without changing its angle branch:

- old boundary edges remain squared length one;
- internal diagonals use `q_j`;
- poles use squared length `-rho_k`;
- final boundary edges use `f_l`;
- every triangle area and its derivatives are computed from its three actual
  edge values;
- internal hinges carry curvature constant `2*pi`, boundary hinges `pi`;
- the same 100 four-simplex orbits and order-24 multiplicities are used.

The resulting restricted gradient has 65 components:

```text
(g_int, p_post),       dim(g_int)=35, dim(p_post)=30.
```

The old pre-momenta are also reconstructed from a separate 30-component
variation of the old boundary lengths at the regular control, but old lengths
remain fixed in the rank problem.

## 4. Frozen evaluator controls

For both phase parities, compare the 100-orbit evaluator with a direct sum over
all 2400 four-simplices at:

```text
B0: q=1, rho=1/4, f=1
B1: q_j=1+(j+1)/1000,
    rho_k=1/4+(k+1)/1000,
    f_l=1+(l+1)/1500
B2: q_j=1-(j+1)/2000,
    rho_k=1/4+(5-k)/1500,
    f_l=1-(l+1)/2500.
```

Require every simplex to remain Lorentzian and the real branch residual below
`2e-7`.  Require relative agreement below:

```text
action                 2e-8
all 65 derivatives     2e-8
triangle curvatures    2e-9.
```

At `B0` and `B1`, centered differences of the complete action in all 65
directions at relative step `3e-5` must reproduce the analytic gradient to
relative `3e-5`.

At `B0`, the first 35 components must reproduce the previously certified
equal-boundary gradient to relative `2e-10`.

## 5. Frozen rectangular rank audit

At `B0`, form the centered Jacobian of the 35 internal per-edge residuals

```text
R = (1/24) g_int
```

with respect to logarithmic variables.  Use relative step `5e-4` and split

```text
J = [J_internal | J_final],
shape(J)          = 35 x 65,
shape(J_internal) = 35 x 35,
shape(J_final)    = 35 x 30.
```

Report singular spectra and ranks at relative thresholds `1e-7`, `1e-9`,
`1e-11`.  Report also:

- condition number of `J_internal`;
- rank and five-dimensional-or-larger left nullspace of `J_final`;
- projection of the regular residual onto `im(J_final)` and its orthogonal
  complement;
- minimum-norm linear steps solving `J delta=-R` and, separately when
  possible, `J_final delta_f=-R`;
- phase-parity differences in all singular spectra.

No rank value is anticipated as a PASS target except dimensional consistency
and stability across the three frozen thresholds.  A rank difference between
parities is a result, not a script failure.

## 6. Boundary Legendre controls

At `B0`, calculate 30 old pre-momentum orbits and 30 final post-momentum
orbits from the full action.  Verify:

- both are real;
- the old/new momentum multisets transform into each other under complete
  time reversal of the staircase slab with the correct action-derivative
  sign convention;
- no momentum is called zero unless its numerical value is below `2e-10` and
  confirmed by direct action differences;
- momentum-orbit equality is imposed only when the schedule stabilizer proves
  it.

This establishes the kinematic canonical relation carried by the slab; it
does not establish an on-shell evolution because `B0` is nonstationary.

## 7. Acceptance and next boundary

The boundary-Legendre route advances if:

1. the 65-variable evaluator passes every full-action control;
2. the rectangular ranks are stable;
3. pre/post momenta are real and obey time reversal;
4. the final-boundary response is nonzero.

Only after committing these results may a continuation/root protocol choose a
boundary perturbation.  That protocol must use the observed rank/nullspace
without selecting a direction because it happens to yield a root.

Failure of `J_final` to have rank 30 does not kill dynamics, but it further
restricts which internal residual directions boundary data can control.
Failure of the complete 65-variable evaluator or time-reversal control kills
this formulation.

Only the targeted verifier will be run.  The full suite remains outside this
mission.
