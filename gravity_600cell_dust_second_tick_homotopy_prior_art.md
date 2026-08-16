# Prior-art gate: canonical target homotopy for the second 600-cell dust tick

Date: 2026-08-16

Accepted first tick: `46a7361`.

Direct second-tick solver boundary: `6346ad0`.

Status: **completed before evaluating any intermediate homotopy target**.

This search establishes context, not external novelty.

## 1. Exact proposed construction

Use the accepted first-tick root as the exact start of a numerical homotopy.
Let `a1=log(L1/L0)`, let `t0` be the incoming momentum target matched by the
accepted first slab, and let `t1` be that slab's mapped post-momentum, which is
the desired incoming target for the second slab.

For `lambda in [0,1]` set

```text
lower_log(lambda) = lambda*a1,
target(lambda)    = (1-lambda)*t0 + lambda*t1.
```

The unknowns are absolute logs

```text
b = log(L_upper/L0),
r = log(rho/rho0),
```

with exact geometry

```text
q_old    = exp(2*lambda*a1)*L0^2,
q_new    = exp(2*b)*L0^2,
pole     = exp(r)*rho0,
diagonal = exp(lambda*a1+b)*L0^2-exp(r)*rho0.
```

At `lambda=0`, `(b,r)=(a1,r1)` is the already accepted first-tick root.  At
`lambda=1`, any accepted endpoint is exactly the desired second-tick problem.
The interpolation is only a numerical path between boundary-value problems;
its intermediate points are not asserted to be physical time slices.

## 2. KNOWN

Canonical simplicial evolution through consecutive moves and propagation of
pre/post data are standard:

- Dittrich and Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>.

The general consistent-discretization mechanism in which discrete equations
fix continuum multiplier-like data is also known:

- Gambini and Pullin, <https://arxiv.org/abs/gr-qc/0511096>.

Curved Regge backgrounds can replace constraints by background-dependent
pseudo-constraints, so tracking the weak singular value along the path is
essential:

- Bahr and Dittrich, <https://arxiv.org/abs/0905.1670>.

Regge FLRW studies explicitly warn that global and local variations can give
different viable equations and that finite carriers can terminate at null
struts; refinement changes the stopping behaviour:

- Liu and Williams, *Regge calculus models of the closed vacuum Lambda-FLRW
  universe*, <https://arxiv.org/abs/1501.07614>.

Predictor/corrector or homotopy continuation is a standard numerical method,
not a new physical principle.  The located Regge sources do not prescribe the
particular linear boundary-data homotopy above.

## 3. CONTROL

The calculation must reconstruct the exact accepted `lambda=0` root from the
committed artifact before advancing.  It must ignore the failed direct
trajectory as a seed and use only the preceding accepted homotopy point.

At every point it must retain:

- the same fixed mass and complete local Regge+dust action;
- the vertex-derived orbit map;
- the same two residual equations and all 65 substitution gates;
- the calibrated rank test on the two-variable Jacobian;
- the Lorentzian/complex-angle branch gate;
- independent even/odd schedule evaluation.

## 4. OPEN difference

Before calculation it is open whether:

1. the accepted first root lies on a connected rank-two branch all the way to
   the second boundary datum;
2. the weak singular scale shrinks, recovers or becomes unresolved;
3. the lapse tends toward zero before `lambda=1`;
4. the endpoint continues contraction, stalls or turns around;
5. both staircase schedules follow the same branch;
6. any endpoint survives refinement or anisotropic perturbations.

No located primary source gives this exact homotopy, endpoint or carrier.
External novelty remains **OPEN**.

## 5. Framing boundary

A passing homotopy would defeat only the direct Newton failure and establish a
connected local second-tick solution in the homogeneous subspace.  The chosen
homotopy is a numerical existence device, not the physical evolution law.

A failed fixed homotopy does not prove nonexistence of every disconnected
root.  Conversely, a passing endpoint does not derive absolute time: `tau0`
and the interpolation parameter remain supplied externally.
