# Prior-art and framing gate: canonical meaning of the internal kernel

Date: 2026-08-22.

Status: **TARGET-DISCLOSED AFTER AN EXPLORATORY RATIO CHECK; FROZEN BEFORE A
REGISTERED SYMBOLIC/PROJECTOR CERTIFICATION.**

## Exact question

Under the complete hypotheses of the replicated finite-height internal-rank
result, the map

```text
R_p=H_p[internal,active]G_p
```

has a common one-dimensional homogeneous kernel.  Determine whether this line
is exactly the tangent to the homogeneous internal lapse constraint `C=0` in
the present carrier coordinates.  Then determine whether the fixed incoming
canonical momentum condition `P=0` annihilates or removes it.

This is a reconciliation calculation.  It makes no external novelty claim.
An unregistered high-precision scratch evaluation already found
`c/sigma=0.4589898592210...`, matching the numerical kernel, and a nonzero
momentum response.  Those values are disclosed rather than treated as blind
evidence.  The registered calculation must derive the identities exactly,
use both independently stored projectors and retain hostile coordinate
controls.

## Repository facts that must be composed, not rediscovered

### Replicated internal rank

The primary and full-real-space verifiers independently give

```text
rank R_diag =119, nullity R_diag =121,
rank R_full =239, nullity R_full =1,
```

for both staircase parities.  The sole complete survivor is homogeneous.

### Exact homogeneous equations

For normalized fixed incoming data `(m,pi)` the certified one-slab equations
are

```text
C(h,q)=8*pi[mu(q)-m]+4*pi*h*q*mu(q),
P(h,q)=p(q)-pi-2*pi*h*mu(q).
```

Here `C=0` is the internal lapse equation and `P=0` fixes the incoming
canonical momentum.  The exact state identity

```text
4*pi*mu'(q)+q*p'(q)=0
```

gives the already certified local Legendre determinant

```text
det partial(C,P)/partial(h,q)=8*pi^2*h*mu(q)^2>0
```

at every positive-height root.

### Present carrier coordinates

On a homogeneous carrier vector, `sigma_v=sigma` and `c_v=c`.  The carrier
protocol fixes

```text
sigma = delta lambda,
c     = delta log rho,
lambda=1+h*q,
rho   =h^2.
```

Therefore the coordinate conversion to be checked is

```text
delta h = h*c/2,
delta q = sigma/h-q*c/2.
```

Its determinant is `-1/2`, so it cannot create or remove a homogeneous
nullity.

## Framing attack

A one-dimensional tangent after imposing one homogeneous internal constraint
on two homogeneous variables is not evidence for a clock or a dynamical
degree of freedom.  It may be only the tangent to `C=0`.

Conversely, if the fixed-input equation `P=0` removes the tangent, this is not
evidence that the discrete theory has no evolution.  A locally regular
boundary-value relation should have an isolated output for fixed incoming
canonical data.  Physical perturbative evolution requires the derivative of
the output with respect to varying incoming canonical data, not a nonzero
fixed-input kernel.

Thus the possible bounded negative is narrow:

> the internal kernel cannot itself be interpreted as a free tick or physical
> mode at fixed incoming data.

It does not close an action-derived forced Jacobi/Legendre map with varying
incoming geometry and momentum.

## Primary literature already gated in the repository

- Dittrich and Höhn's action-generated canonical simplicial evolution
  distinguishes pre/post constraints and propagating data
  ([arXiv:1108.1974](https://arxiv.org/abs/1108.1974)).
- Their Regge analysis explains background-dependent pseudo-constraints
  ([arXiv:0912.1817](https://arxiv.org/abs/0912.1817)).
- Bahr and Dittrich show why curved Regge backgrounds generically break exact
  lapse gauge symmetry ([arXiv:0905.1670](https://arxiv.org/abs/0905.1670)).

These works make the constraint-versus-evolution distinction known.  They do
not compute the present 600-cell carrier ratio.

## Falsifiable outcomes

### `INTERNAL_KERNEL_IS_LAPSE_CONSTRAINT_TANGENT_FIXED_INPUT_REMOVES_IT`

Use only if an exact symbolic coordinate calculation derives the `C`-tangent
line, both independently reconstructed kernel projectors agree with it, and
the derivative of `P` on that line is certified nonzero.

This is **DERIVED EXACT/COMPUTATIONAL** and a **BOUNDED NEGATIVE** for reading
the line itself as fixed-input evolution.

### `INTERNAL_KERNEL_NOT_THE_HOMOGENEOUS_CONSTRAINT_TANGENT`

Use if all provenance and exact-coordinate controls pass but the replicated
projector is resolved different from the exact `C` tangent.  This refutes the
proposed reconciliation and reopens the carrier/Hessian interpretation.

### `CANONICAL_RECONCILIATION_OPEN` or `CONTROL_FAILED`

Use for a numerical comparison in its open band or any provenance, formula,
coordinate, branch or hostile-control failure.

No outcome derives a physical tick, `c`, `G`, Planck units, gravitons or
particle physics.
