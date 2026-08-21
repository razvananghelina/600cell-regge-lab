# Consolidated result: local Regge curvature selects an on-shell refined mass

Date: 2026-08-21

## Headline

> **DERIVED COMPUTATIONAL / STRUCTURAL, adversarially corroborated after a
> disclosed post-hoc discovery:** on the induced static product over
> `K0=P(sd K_600)`, every local vertical Regge equation is exactly the
> endpoint-half spatial curvature balance
>
> ```text
> dS_grav/dlog(rho_v)=tau0*K_v/2,
> m_v=K_v/(8*pi).
> ```
>
> The curvature-matched masses are positive, conserve the previously selected
> total mass and make the inherited product fill stationary on every one of
> the 24 staircase schedules.

This bypasses the deferred equal-boundary root census for a new, explicitly
different matter branch.  It does not rescue the earlier homogeneous `P1`
dust branch.

## Provenance and reproduction

| stage | commit |
|---|---|
| prior-art gate plus exploratory disclosure | `f12f56c` |
| primary protocol | `bb80f28` |
| primary verifier registered before execution | `e76f756` |
| first `13/15` control failure preserved | `ead1a1c` |
| narrow control correction preregistered | `2db5457` |
| corrected primary implementation | `1fa48a2` |
| primary result | `2881e4c` |
| adversarial protocol | `f4336e7` |
| adversarial verifier registered before execution | `46c2c0b` |

Primary verifier and artifact:

```text
reproducible/verify_gravity_600cell_refined_local_curvature_mass.py
reproducible/gravity_600cell_refined_local_curvature_mass.json
SHA-256 180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
15/15 PASS twice, byte-identical.
```

Adversarial verifier and artifact:

```text
reproducible/verify_gravity_600cell_refined_local_curvature_mass_adversarial.py
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
SHA-256 c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
16/16 PASS twice, byte-identical.
```

The static registry audit reports 372 entries, 372 distinct names, zero
duplicates, zero unregistered verifiers, zero stale registrations and two
reasoned deliberate exclusions.  No full suite and no nested root census was
run.

## Exact object and complete hypotheses

The carrier is the exact projected barycentric 600-cell boundary with

```text
f=(2640,17040,28800,14400).
```

Use its rank-derived chordal metric, the induced static Minkowski product at
the supplied `tau0=0.0102`, the corrected complex Lorentzian Regge action with
boundary terms and any of the 24 colour-ordered staircase triangulations.
The matter action consists of conserved vertex masses on the vertical
worldlines,

```text
S_dust=-8*pi sum_v m_v sqrt(rho_v).
```

For every spatial edge put `C_e=l_e epsilon_e` and define the endpoint-half
curvature

```text
K_v=(1/2) sum_(e incident on v) C_e.
```

This is a declared weak vertex localization of hinge curvature.  It is not a
claim that endpoint-half localization is the unique scalar-curvature density
for arbitrary irregular Regge meshes.

## Primary route

The primary reconstruction used only the six exact rank-pair hinge types to
form the four rank totals `K_r`.  Independently, it removed the analytically
known `P1` dust derivative from the already certified 100/140-digit action
residuals.  Across all `24*4=96` vertical equations it found

```text
max |G_grav,r-tau0*K_r/2| = 1.4016460e-76,
schedule spread                    = 0.
```

The matter response is exactly `-4*pi*tau0*I_4`, so it has rank four and
selects a unique mass vector.  The four mass fractions are

```text
(0.1287831657723389984083291990898173...,
 0.3657000761313201399462782272254918...,
 0.3759856918014127686434122686089311...,
 0.1295310662949280930019803050757598...).
```

Their sum is one and every entry is positive.

## Mechanically independent route

The adversarial verifier imported no primary function and did not subtract a
`P1` residual.  It reconstructed the actual 2,640 vertices, 17,040 edges and
14,400 chambers, computed every edge incidence and accumulated half of every
`l_e epsilon_e` at each actual endpoint.  It found:

```text
within-rank curvature spread                   0,
endpoint/edge curvature conservation error     1.54e-101 relative,
maximum disagreement with primary fractions    4.74e-76.
```

It then parsed only the old stationary-fill action definitions, set the dust
mass to zero and directly evaluated the complete Lorentzian action for a
schedule and its reverse.  The result was

```text
max |G_grav,r-tau0*K_r/2|                      4.63e-98,
max residual after adding m_r=K_r/(8*pi)       4.63e-98,
Richardson/analytic derivative disagreement     1.44e-52 relative.
```

The known unrefined 600-cell control returned exactly `1/120` curvature per
vertex; dropping one refined endpoint contribution failed conservation; and
the `P1` masses retained a direct residual `0.038204082...`.

## Logical consequence for stationarity

The earlier all-schedule census already certified that every cross-diagonal
equation vanishes at the induced product and that every actual edge of a
given state type has one unmixed incidence signature.  Dust changes only the
vertical equations.  The present actual-incidence audit proves that `K_v` is
constant within each rank orbit.  Therefore assigning

```text
m_v=K_v/(8*pi)
```

cancels every individual vertical equation, not just the sum of four orbit
variations.  Hence the inherited induced product is an internal stationary
fill for all 24 schedules on this curvature-matched matter branch.

This is the refined local version of the previously derived global static
balance `M=K/(8*pi)`.

## Framing attack and physical cost

The positive result does not mean that the refined geometry supports the same
homogeneous dust approximation previously assumed.  The affine-exact `P1`
dual volumes put one quarter of volume and mass in every rank, while the
curvature-matched density relative to those volumes is

```text
(0.5151326631, 1.4628003045, 1.5039427672, 0.5181242652),
```

with maximum/minimum ratio `2.919525...`.  This is a large finite-lattice
rank contrast, not a small perturbation.

Accordingly:

- **DERIVED NEGATIVE:** homogeneous `P1` dust is not locally on shell on this
  fixed projected carrier;
- **DERIVED / STRUCTURAL:** solving the local time-symmetric constraint
  selects positive curvature-matched dust;
- **OPEN:** whether this selected density converges toward uniform density on
  the next canonical refinement;
- **FORBIDDEN:** calling the four fractions particle multiplets or fitting a
  familiar pattern to them.  The nearby `(1,3,3,1)/8` vector is explicitly
  refuted as an exact identity.

The discovery provenance is explicitly post-hoc: the numerical equality was
seen before the formal protocol.  Its evidential force comes from the general
local balance and independent action/incidence derivations, not from pretending
the comparison was blind.

## Literature reconciliation

Hinge-supported Regge curvature is standard since Regge's original action
(DOI `10.1007/BF02733251`).  Vertex scalar-curvature localizations using dual
cells are developed by McDonald and Miller
([arXiv:0805.2411](https://arxiv.org/abs/0805.2411)).  Local versus global
variation and the time-symmetric initial-value constraint in closed Regge
cosmology are analyzed by Liu and Williams
([arXiv:1501.07614](https://arxiv.org/abs/1501.07614)).  Conserved simplicial
dust worldlines with action minus mass times proper time are used by Dittrich,
Gielen and Schander
([arXiv:2109.00875](https://arxiv.org/abs/2109.00875)).

A post-result search using the terms `static product strut variation`,
`endpoint half edge curvature`, `vertex mass` and `Hamiltonian constraint`
located no primary source proving this exact projected-600-cell four-rank
identity.  Search absence is not proof; external novelty remains **OPEN**.

## Next decisive gate

The stationary seed now supplies six refined boundary momentum orbits from
the action.  The next target-blind calculation is therefore small:

1. compute those six old/new boundary covectors for every schedule;
2. apply the already proved coarse/fine cotangent pullback;
3. compare the pulled-back homogeneous momentum with the accepted coarse
   static momentum at the same induced coarse geometry.

If they agree, the action has selected one concrete point in the previously
five-dimensional lift fiber and the refined canonical route advances.  If
they disagree, the curvature-matched seed is a different coarse canonical
state and does not repair transport.  Only after that gate may a refined
boundary Hessian or mode spectrum be considered.

No tick, `c`, `G`, Planck scale or particle mass is derived here.
