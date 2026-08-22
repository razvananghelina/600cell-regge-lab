# Result: one homogeneous internal-constraint tangent survives

Date: 2026-08-22.

Status: **DERIVED COMPUTATIONAL, MECHANICALLY DIFFERENTLY REPLICATED;
PHYSICAL CLASSIFICATION OPEN.**

## Complete hypotheses

Fix the first positive-height homogeneous 600-cell dust slab generated from
the incoming state `v=3/2`, with zero cosmological constant, conserved total
dust mass `M=mu(v)`, positive proper height, positive endpoint scales and the
repository's fixed Lorentzian branch.  Fix either of the two certified H4
staircase parities and the exact geometry-selected rank-240 scale-plus-strut
carrier

```text
G_p : R^240 -> R^1560.
```

Hold the old spatial boundary fixed.  Let the `840` internal action equations
be the `720` oriented cross-diagonal equations and `120` same-vertex
pole/strut equations.  The tested map is

```text
R_p = H_p[internal,active] G_p : R^240 -> R^840.
```

It is the derivative of internal stationarity along this carrier.  It is not
the complete canonical boundary evolution equation.

## Reproducible verdict

The preregistered 180-digit identity-row orbit calculation reports

```text
FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_SELECTED_PRIMARY
25/25 PASS
```

and a mechanically different calculation using complete physical
`2280 x 2280` Hessians, direct global binary64 SVD/QR with an explicit
roundoff allowance, and nonlinear action-gradient secants reports

```text
FINITE_HEIGHT_INTERNAL_CARRIER_KERNEL_ADVERSARIALLY_REPLICATED
19/19 PASS
```

The first adversarial execution stopped after its even rank census because a
numeric prediction error had been formatted as text before applying its
already frozen threshold.  The failure is preserved.  Commit `c71750e`
changes only the value's operational type through classification.

Accepted commits:

```text
primary implementation        1472213
primary artifacts             2b9c124
adversarial protocol          b28c4a9
adversarial implementation    5065a7c
failure preservation/correction c71750e
adversarial artifacts         3ac04d5
```

Accepted artifact hashes:

```text
primary JSON
513fdea33f6b868efa6d6f2b2526bade7ce615ea949f955588916a8d0baee0c8

primary matrices
97f5b8318be2b3ccf843db87e678ac1ac6ce402db262023c6bbc63a7b647321b

adversarial JSON
ddd4704b7d1deb6360e752b2ebfe5cc0b66d03819c9f8df7b74e24373aa98fb5

adversarial matrices
45ee09642485dc0e18c6378b9454882414bb07d7e5ad9dbd0c3a8896fd8a7f74
```

Only the two targeted verifiers and static registry checks ran.  The full
suite did not run.

## Rank ledger

Both methods and both parities give

| Map | Shape | Rank | Nullity |
|---|---:|---:|---:|
| oriented diagonals only | `720 x 240` | `119` | `121` |
| diagonals plus poles | `840 x 240` | `239` | `1` |

The pole equations remove every diagonal-only survivor except one
homogeneous direction.  All nontrivial 2T sectors have full column rank in
the complete map.  The sole kernel occurs in the trivial sector.

The adversarial global classifier has

```text
conservative normalized error       5.29008467897e-8
normalized null singular value      about 2.7e-13
next normalized singular value      about 4.30e-3
```

Thus the decision is separated from both classifier bands by many orders of
magnitude.  Pivoted QR independently gives rank `239` and `119` under the
same interval.

## Kernel identity and direct-action control

In the common physical source order

```text
(sigma_0,...,sigma_119,c_0,...,c_119),
```

the kernel vector is homogeneous to binary64 reconstruction error:

```text
all sigma_v equal,
all c_v equal,
c/sigma = 0.4589898592210...
```

Its even/odd projector distance is

```text
4.09083479608e-12,
```

below the preregistered agreement gate

```text
9.32587340685e-10.
```

The primary/real-space projector distances are about `2.3e-11`, also below
that gate.

The direct nonlinear control reevaluates the complete action gradient on all
`2400` four-simplices at three centered secant steps.  Relative to the
smallest nonzero singular direction, the finest kernel secant is suppressed
by

```text
even  5.5225e-11,
odd   5.4810e-11.
```

The nonzero-direction secant agrees with the assembled Hessian response to
about `7.9e-11`.  The step differences decrease by the expected factor four.
This rules out an orbit-Fourier or Hessian-assembly cancellation at the
tested precision.

## What follows, and what does not

**DERIVED COMPUTATIONAL:** within the complete stated hypotheses, internal
stationarity eliminates every nonhomogeneous direction of the fixed
scale-plus-strut carrier at this finite-height slab and retains exactly one
schedule-independent homogeneous tangent.

**STRUCTURAL:** the survivor is a tangent to the internal constraint surface.
It is not yet a canonical solution tangent because no fixed incoming momentum
condition has been imposed in `R_p`.

**OPEN:** whether the observed homogeneous line is exactly the differential
of the closed homogeneous lapse constraint, whether it integrates to a
nearby internal-stationary curve with the required boundary data, and whether
the canonical momentum equation removes it.

This result does not contradict the older complete-carrier no-go: that test
used a different background and the complete canonical graph.  The present
test uses the new finite-height background and internal equations only.

It does not derive a graviton, local propagation, a wave equation, stability,
a continuum limit, a physical tick, `c`, `G`, Planck units, particle masses
or Standard-Model physics.

## Exact next gate

Restrict the already certified homogeneous one-slab equations to the common
kernel coordinates and compare the kernel with the exact differential of

```text
C(h,q)=0.
```

Then impose the fixed incoming canonical momentum equation

```text
P(h,q)=0.
```

The repository already contains the identity

```text
det partial(C,P)/partial(h,q)=8*pi^2*h*mu(q)^2>0
```

at every positive-height root.  A registered exact reconciliation must check
that its coordinate conversion is the present carrier's
`(sigma,c)` conversion and that the observed projector is the `C` tangent.
If so, `P` removes the last line and this carrier supplies no nonzero fixed-
incoming canonical tangent at the finite-height slab.  That would be a clean
bounded negative, not evidence for a tick.
