# Disclosure: exploratory derivation of the full scale--strut edge response

Date: 2026-08-19

## Why this file exists

During the coordinate audit for the next gravity mission, an exact symbolic
calculation was performed before the mission-specific prior-art gate and
preregistration had been committed.  No action, Hessian, strong-equation
matrix, representation target, continuum label, or desired rank was loaded.
Nevertheless, the calculation exposed the candidate formula that a later
verifier must test.  It therefore cannot be presented as a blind first result.

The status of everything below is:

```text
DERIVED EXPLORATORY; TARGET-DISCLOSED;
NOT ACCEPTED UNTIL A PREREGISTERED INDEPENDENT VERIFIER PASSES,
FOLLOWED BY A MECHANICALLY DIFFERENT ADVERSARIAL REPLICATION.
```

This disclosure freezes the entire exploratory observation before any new
action calculation is authorized.

## Exact local object and conventions

Use a regular tetrahedron with spatial vertices `p_i`, normalized by

```text
|p_i-p_j|^2 = 8,     |p_i|^2 = 3,
```

and Lorentz metric `diag(1,1,1,-1)`.  Its homothetic upper copy is

```text
q_i = lambda p_i + tau n,     n^2=-1.
```

The four `sigma_i` are infinitesimal upper radial/scale data, normalized so
that an upper edge has raw squared-length response

```text
delta |q_i-q_j|^2 = 8 lambda (sigma_i+sigma_j).
```

The four `s_i` are raw signed squared-length variations of the four struts
`p_i q_i`.  For an oriented cross diagonal from lower vertex `i` to upper
vertex `j`, `i!=j`, the endpoint-supported ansatz was

```text
delta d_(i->j)^2 = A sigma_i + B sigma_j + C s_i + D s_j.
```

Complete hypotheses for the generic calculation are

```text
lambda != 1,
tau != 0,
(lambda-1)^2 - 3 tau^2 != 0,
the local 22 x 16 squared-length Jacobian has its generic rank 16,
the exact affine transition across the shared lateral face exists.
```

The final denominator is the non-null lateral-face factor
`lambda^2-2 lambda+1-3 tau^2`.  No static `lambda -> 1` limit is asserted.

## What one frustum determines

An exact `22 x 16` Jacobian was built from the six lower edges, six upper
edges, four struts and twelve oriented cross diagonals of one tetrahedral
frustum.  Substitution of the endpoint-supported ansatz into its exact left
kernel gave only

```text
A+B = 8,
C+D = 1.
```

**DERIVED EXPLORATORY.** A single frustum does not select the four endpoint
coefficients.  Any claim of local uniqueness would be false.

## What shared-face consistency appeared to determine

Two regular tetrahedra `(0,1,2,3)` and `(0,1,2,4)` were developed on opposite
sides of their common face `(0,1,2)`.  The exact affine Lorentz transition
between their frames had determinant `-1`.  Its lower-face Poincare
stabilizer had dimension one and was retained during elimination rather than
set to zero by hand.

After eliminating that connection column, exact shared-face consistency
reduced to

```text
B (lambda-1)^2 = 2 (lambda-1)^2 + 2 tau^2,
(lambda-1) D   = lambda.
```

Together with the one-frustum equations this gives the unique candidate

```text
A = 6 - 2 tau^2/(lambda-1)^2,
B = 2 + 2 tau^2/(lambda-1)^2,
C = -1/(lambda-1),
D =  lambda/(lambda-1).
```

At the two exact rational controls already present in the frozen universal
local-lift artifact, this gives

```text
(lambda,tau)=(2,5):
    (A,B,C,D)=(-44,52,-1,2),

(lambda,tau)=(3,11):
    (A,B,C,D)=(-109/2,125/2,-1/2,3/2).
```

These values agree with the physical responses reconstructed from the frozen
exact `6 x 8` local blocks.  Agreement with already-known controls is not an
independent proof of the generic formula.

## Relation to the accepted two-frustum theorem

The accepted two-frustum theorem says that two *local Poincare motions*
compatible on a shared face reduce to the diagonal common motion; no hidden
relative face mode survives.  The present candidate concerns a different
map: endpoint boundary data `(sigma,s)` to cross-diagonal squared-length
responses after compatible local representatives are glued.

Thus the statements are compatible.  The old theorem removes an arbitrary
relative Poincare choice; the exploratory calculation suggests that this
removal is precisely what fixes the two coefficients left free by one cell.
This reconciliation is **STRUCTURAL** until the new verifier is accepted.

## Conversion to the frozen curved slab

Let the physical old boundary edge satisfy `L0^2=8 alpha^2`, let `rho>0` be
the magnitude of the timelike pole squared length, and set

```text
tau^2 = 3 (lambda-1)^2 + 8 rho/L0^2,
q_diag = lambda L0^2-rho.
```

Use action coordinates

```text
sigma_v                 additive vertex-scale data,
c_v = delta log(rho_v)  logarithmic strut-magnitude data.
```

The candidate `1560 x 240` pre-Legendre carrier has rows consisting of 840
internal edges (720 oriented cross diagonals and 120 poles) followed by 720
new-boundary edges.  Its nonzero log-squared-length responses are:

```text
pole v:
    delta log(rho_v) = c_v,

new boundary edge {u,v}:
    delta log q_new = (sigma_u+sigma_v)/lambda,

cross diagonal u(lower) -> v(upper):
    delta log q_diag
      = (L0^2/8q_diag) (A sigma_u+B sigma_v)
        + rho/((lambda-1)q_diag) (c_u-lambda c_v),

A = -16 rho/(L0^2 (lambda-1)^2),
B = 8+16 rho/(L0^2 (lambda-1)^2).
```

The strut half is exactly the already accepted corrected pure-strut formula.
The sum of all scale columns gives the homogeneous derivative

```text
cross diagonal:  L0^2/q_diag,
new edge:        2/lambda,
pole:            0.
```

These are proposed controls, not yet a frozen result.  Near the accepted
background `lambda` is close to one, so individual scale coefficients are
large and the later numerical conditioning audit is load-bearing.

## Scientific firewall

Even if the carrier is accepted, it is kinematic.  It does not prove that
any of its columns solves the canonical equations, is gauge, propagates, or
defines time.  It does not select a tick, `c`, `G`, a Planck scale, a graviton
or a particle mass.  Those questions begin only after the geometric carrier
is frozen and pulled through the independently frozen action derivative.

