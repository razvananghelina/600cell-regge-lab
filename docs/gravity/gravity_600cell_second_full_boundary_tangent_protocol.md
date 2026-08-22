# Preregistration: second finite-height full-boundary tangent and two-step map

Date: 2026-08-22

Prior-art gate commit: `eae691a`.

Status: **FROZEN BEFORE THE FIRST SECOND-SLAB HESSIAN, PRE-LEGENDRE SINGULAR
VALUE, DETERMINANT BALL, TANGENT ENTRY OR TWO-STEP COMPOSITION.**

Only the mission-specific verifier and its already registered 43 geometry
controls may run.  The full suite is forbidden.

## 1. Frozen inputs

```text
docs/gravity/gravity_600cell_second_full_boundary_tangent_prior_art.md
  d3740e0b08b2f3ec6adf2c69c762e5e5dc0cdd87a571d6d27bc62e78518e70be

reproducible/gravity_600cell_finite_height_asymptotic_map.json
  a93837d2bbec340ddbac528c0be4da52aefe45c8f0d4310496eb1aef6a7b19b6

reproducible/verify_gravity_600cell_finite_height_asymptotic_map.py
  3aafdb326eb9299d9e69ef79c0726eeb09f214b9dee1dc848848e34e0920b208

reproducible/gravity_600cell_finite_height_full_boundary_tangent.json
  266638aeaa825b327b63a84eda36a499456dc4b4f9a86f964cee5f79d6d6e930

reproducible/gravity_600cell_finite_height_full_boundary_tangent.npz
  0c34f179821f9d0b74de4906051bbcb7149b4e79881410ea662241adc0aa19bf

reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent.py
  c4e60d6ef87131d87a93b64d5381d16d8de8d3990340efd5405ec983f64db94d

reproducible/gravity_600cell_finite_height_full_boundary_tangent_adversarial.json
  ee9491b2ae5fdf3f2a9d0d78c0e837c8c2692797d87ccd8e1757efeadd8060e7

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py
  834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf
```

Require the first tangent to retain `21/21`, all 14 sector maps regular and
canonical, and `SCHEDULE_ROBUST`.  Require its dense adversarial artifact to
retain `22/22` and the replicated schedule-robust outcome.  These are input
and formula controls; their tangent entries remain unopened until the new
second-slab labels are frozen.

The new verifier must be registered exactly once, with no duplicate registry
names.

## 2. Independent two-slab homogeneous reconstruction

At 180 decimal digits define

```text
epsilon(q)=2*pi-5*acos((q^2+2)/(2(q^2+3))),
mu(q)=180*epsilon(q)/(pi*sqrt(q^2+4)),
p(q)=180*q*epsilon(q)/sqrt(q^2+4)
     -600*sqrt(3)*asinh(q/sqrt(8*(q^2+3))).
```

For incoming normalized `(m,pi)` define

```text
E(q;m,pi)=4*pi*(mu(q)-m)+q*(p(q)-pi).
```

Set `v=3/2`, `m0=mu(v)`, `pi0=p(v)`.  Use deterministic bisection only on
the frozen brackets

```text
q1 in (9,10),
q2 in (31,32).
```

After each root set

```text
h=(p(q)-pi)/(2*pi*mu(q)),
r=1+h*q,
m_next=m/r,
pi_next=p(q)+2*pi*h*mu(q)/r.
```

Require both residuals below `1e-140`, bracket widths below `1e-150`, all
`h,r,m` positive, junction identities below `1e-140`, and agreement with the
committed branch-B history below `1e-65`.  No other second root or branch may
replace this state.

Record without rounding the derived `q1,h1,r1,m1,pi1,q2,h2,r2,m2,pi2`.

## 3. Exact scale-lift theorem and hostile family

In logarithmic signed-squared-edge coordinates, prove symbolically from

```text
S_phys(x;M0)=r1^2 S_norm(x-2*log(r1);m1)
```

that the Hessians and logarithmic canonical momenta scale by `c=r1^2`.
Define on a phase block with `n=30d` configuration coordinates

```text
Omega=[[0,I_n],[-I_n,0]],
D_c=diag(I_n,c I_n).
```

Certify exactly

```text
D_c^T Omega D_c=c Omega,
(D_c T D_c^-1)^T Omega (D_c T D_c^-1)=Omega
```

for a symbolic symplectic block matrix `T`.

Also certify that scale covariance does not uniquely select a canonical
frame: `diag(sI,s^-1 I)` and `[[I,0],[bI,I]]` are symplectic for arbitrary
nonzero `s` and arbitrary `b`.  This is a framing control, not a positive
physics result.

The wrong conjugacies using `diag(I,r1 I)` and the identity must be retained
as hostile alternatives.  At least one actual sector must separate each from
the directly assembled physical-unit map above 100 times its comparison
uncertainty.

## 4. Common carrier before any second-slab Hessian

For each staircase parity construct the complete carrier

```text
O=720, X=840, N=720,
2400 four-simplices, 6240 triangles,
95 orbit types = 30+35+30.
```

Construct the state-independent right-regular sectors separately at both
`r1` and `r2`.  Require literal equality, before any Hessian result, of:

- the order-24 action tables;
- the 720 physical final-to-old edge shifts;
- all old orbit seeds;
- sector dimensions `[1,1,1,2,2,2,3]`;
- splitters and central eigenvalues;
- every minimal-basis entry below `1e-140`.

This is required before composing with the first-slab archive.

## 5. Second normalized and direct physical Hessians

Use 180 decimal digits and the frozen centered logarithmic derivative steps

```text
1e-25, 5e-26, 2.5e-26, 1.25e-26.
```

For the normalized second slab use

```text
old=1,
internal=r2-h2^2,
pole=-h2^2,
new=r2^2,
mass=m1.
```

Separately rerun every local derivative and assemble a direct physical-unit
kernel with

```text
old=r1^2,
internal=r1^2*(r2-h2^2),
pole=-r1^2*h2^2,
new=r1^2*r2^2,
mass=m0.
```

The physical evaluator must not obtain its kernel by multiplying the
normalized one.  At all four raw levels and all three Richardson levels
classify

```text
K_phys = r1^2 K_norm
```

using the sum of both adjacent Richardson variations, both arithmetic
floors and direct Flint/binary conversion bounds.  Equality inside ten times
that frozen uncertainty is `SCALE_LIFT_CONFIRMED`; separation above 100 times
is `SCALE_LIFT_REFUTED`; the gap is `SCALE_LIFT_OPEN`.

Require the same Lorentzian inertia, positive leading-minor and
angle-argument controls as the first tangent.  Require raw projected
reciprocity before any symmetrization.

## 6. Second pre-Legendre and canonical classifiers

For every parity, sector and Richardson level form

```text
J=[[K_XX,K_XN],[-K_OX,-K_ON]].
```

Use exactly the scaled-SVD plus 140-decimal Flint determinant-ball
classifier from the first tangent protocol.  A sector is `REGULAR` only if
all three scaled SVD gaps exceed `100*e_J` and all determinant balls exclude
zero.  Otherwise it is `NUMERICALLY_OPEN`.

For every regular block construct

```text
R=[[-K_XO,0],[K_OO,I]],
Y=J^-1 R,
T=[[Y_N],[[K_NO,0]+K_NX Y_X+K_NN Y_N]].
```

Apply only the exact physical boundary relabelling.  Use the same complex
ball symplectic and determinant controls as the first protocol.  No fitted
phase, basis rotation or Schur coefficient is allowed.

For every level require the independently solved physical tangent to agree
with

```text
D_c T_norm D_c^-1
```

under a ball-radius, Richardson-variation and binary-roundoff uncertainty
frozen before comparison.

## 7. Second-slab schedule classification

Compare the even and odd direct physical blocks in the literally common
minimal bases.  For each sector use the same three-way classifier:

```text
ROBUST     if distance <= 10*uncertainty,
DEPENDENT  if distance > 100*uncertainty,
OPEN       otherwise.
```

The global second-slab label is dependent if any sector is dependent, open
if none is dependent and at least one is open, and robust otherwise.  Freeze
all scale-lift, rank, canonicality and schedule labels before opening a first
slab tangent entry.

## 8. Four physical two-step compositions

Only now open the first tangent archive.  For every minimal sector and each
schedule pair

```text
(p1,p2) in {even,odd} x {even,odd}
```

form at every Richardson level

```text
T20[p1,p2]=T2_phys[p2] T1_phys[p1].
```

Propagate entrywise complex-ball radii by

```text
R_BA <= |B| R_A + R_B |A| + R_B R_A.
```

Compare all six pairs of the four operational maps.  The uncertainty is the
sum of both maps' adjacent Richardson variations propagated through the
product, all product ball radii, and a conditioning-aware binary roundoff
term.  Assign `TWO_STEP_SCHEDULE_ROBUST`, `DEPENDENT`, or `OPEN` with the same
10/100 classifier.

Check the symplectic block identities for all four product midpoints inside
the propagated uncertainty.  Require a synthetic `1e-3` schedule corruption
to be detected.  Also require the identity-lift and `r1`-rather-than-`r1^2`
hostile products to differ from the direct physical product above 100 times
the applicable uncertainty in at least one sector.

No eigenvalue or tangent singular-value spectrum is computed.

## 9. Outcome hierarchy

Assign exactly one:

1. `SECOND_FULL_BOUNDARY_TANGENT_CONTROL_FAILED` for failed provenance,
   history, branch, carrier, basis, exact algebra or hostile controls;
2. `SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_REFUTED` if direct physical
   assembly refutes the degree-two lift;
3. `SECOND_FULL_BOUNDARY_TANGENT_SCALE_LIFT_OPEN` for its comparison gap;
4. `SECOND_FULL_BOUNDARY_TANGENT_RANK_OPEN` if any second block is not
   certified regular;
5. `SECOND_FULL_BOUNDARY_TANGENT_CANONICALITY_FAILED` if a regular tangent or
   its direct physical conjugacy fails;
6. `SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_DEPENDENT` if the second slab or
   two-step composition is schedule dependent;
7. `SECOND_FULL_BOUNDARY_TANGENT_SCHEDULE_OPEN` for either schedule gap;
8. `TWO_STEP_FULL_BOUNDARY_TANGENT_SCHEDULE_ROBUST_PRIMARY` only if every
   gate passes.

A primary positive remains **PRIMARY ONLY** until a separately preregistered,
mechanically different dense real-space replication succeeds.

## 10. Interpretation boundary

Outcome 8 establishes only a two-step, first-order, action-generated
canonical response along the representative homogeneous branch.  It does
not make a transfer-matrix spectrum physical, perform the curved-background
constraint reduction, identify gravitons, prove stability or convergence,
or derive a wave equation, limiting speed, `G`, Planck units or particles.

The homogeneous invariant-region theorem still proves a unique successor at
every later finite step only.  Convergence, infinite total proper duration
and completeness remain **OPEN**.

