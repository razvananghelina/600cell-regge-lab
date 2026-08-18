# Preregistered protocol: Lorentzian tent pre/post Legendre map

Date: 2026-08-12

Status at registration:
**PROTOCOL ONLY -- THE MIXED-HESSIAN RANK HAS NOT BEEN COMPUTED**

## 1. Question

The asymmetric tent has a certified real stationary pole and a nonzero
pole-only Hessian. The next question is whether the complete local Regge
action generates a regular canonical relation between the twelve old and
twelve new cone-edge lengths.

This is the correct replacement for the overstrong demand that bare symmetry
choose one unique boundary state. A discrete action may evolve supplied
initial data while remaining fully invariant.

## 2. Complete geometric hypotheses

Use the same 4-ball

```text
T_v=[v,v']*L_v
```

and the same witness certified in commits `cc71574` and `2e92b4a`. Work in
units `a=1`. Introduce squared-length variables

```text
p_u = s(v,u)       for the 12 old cone edges,
q_u = s(v',u)      for the 12 new cone edges,
s(v,v')=-rho       for the internal pole,
s(u,w)=1           for every link edge.
```

At the frozen point,

```text
p_u=1,
rho=1/4,
q_u=(x*,3/2,4/5,3/2) on shells (1,5,5,1),
x*=0.44333089835748125745... .
```

The decimal only locates the already certified bracket root. It is not an
exact input constant.

All final tetrahedra are spacelike, every four-simplex is Lorentzian, every
internal hinge is timelike, and the volume/cosmological, matter and
higher-curvature terms vanish. The link lengths remain fixed in this local
canonical test.

## 3. Frozen complex-angle convention

Use the `+` complex-angle branch and principal square-root/logarithm convention
of Borissova and Dittrich,
[arXiv:2303.07367](https://arxiv.org/abs/2303.07367), appendix A. In their
notation,

```text
i S_Regge^+ = sum_h sqrt(+V_h) epsilon_h^+.
```

For a four-simplex, construct the extended Gram matrix of outward facet
normal covectors. For the hinge opposite facets `a,b`, compute

```text
C_ab=-N_ab/[sqrt(+N_aa)sqrt(+N_bb)],
theta_ab^+=-i log_-[C_ab-i sqrt(1-C_ab^2)].
```

The verifier must recover the already certified real interior angle at every
timelike internal hinge, with `theta^+=-theta_real`.

Use the hinge census

```text
12 internal triangles [v,v',u]:       k=2, incidence 5,
30 old boundary triangles [v,u,w]:    k=1, incidence 2,
30 new boundary triangles [v',u,w]:   k=1, incidence 2,
20 common-link corner triangles:      k=0, incidence 1.
```

Thus

```text
epsilon_h^+=pi*k+sum_(sigma contains h) theta_(sigma,h)^+.
```

The exploratory branch audit performed before registration found the first
three curvature types respectively real, purely imaginary and purely
imaginary as required for a real Lorentzian action. The common-link `k`
choice affects only an additive constant while all link lengths are fixed;
the verifier must confirm that changing `k=0` to `k=1` leaves every derivative
used below unchanged.

Define the real Lorentzian action

```text
S_L=-i (i S_Regge^+).
```

No analytic continuation of the old Euclidean code is admissible as sole
evidence.

## 4. First-derivative controls

Use the complex Schlaefli identity to compute

```text
partial_x S_L
 =-i sum_h [partial_x sqrt(+V_h)] epsilon_h^+.
```

The verifier must:

1. reconstruct all `92` hinges with incidences `(12x5,60x2,20x1)`;
2. show `S_L` and all 25 derivatives are real to `2e-10`;
3. reproduce the internal pole equation `S_rho=E=0` and the certified
   `S_rhorho=6.37870935921427...`;
4. compare every Schlaefli gradient component with centered differences of
   the full complex action at relative tolerance `2e-6`;
5. report the pre/post momenta in squared-length coordinates,

   ```text
   P_u^-=-partial S_L/partial p_u,
   P_u^+=+partial S_L/partial q_u.
   ```

They must be constant on the four stabilizer shells where the boundary data
have that symmetry.

## 5. On-shell mixed Hessian

Let `b=(p_1,...,p_12,q_1,...,q_12)`. First compute the complete `25x25`
Hessian in `(b,rho)` by centered differences of the Schlaefli gradient.
Require:

- symmetry residual below `2e-5` relative to its maximum entry;
- stable entries across step sizes `2e-5`, `1e-5`, `5e-6`, with maximum
  relative discrepancy below `2e-4`;
- agreement of its pole entry with the Arb result to relative `2e-6`.

Eliminate the regular internal pole by the exact Schur-complement formula

```text
H_eff = S_bb-S_b,rho (S_rho,rho)^(-1) S_rho,b.
```

The `12x12` old/new mixed block

```text
W = partial^2 S_on-shell/(partial p partial q)
```

is the local Lagrangian two-form in the star sector. Compute its singular
values for all three step sizes.

Preregistered numerical rank rule:

- a singular value counts nonzero only if it exceeds `1e-7` times the largest;
- the rank and every singular value above threshold must be stable to relative
  `5e-3` across the three steps;
- a smallest/largest ratio within a factor `10` of the threshold is
  **INCONCLUSIVE**, not rounded to either verdict.

The rank is invariant under switching from positive squared lengths to
lengths because this left/right coordinate change multiplies `W` by
invertible diagonal Jacobians.

## 6. Decision boundary

- **DERIVED REGULAR LOCAL LEGENDRE MAP:** `rank(W)=12` with the frozen gap and
  stability controls. Nearby old configuration plus admissible old momenta
  locally determine the new star configuration after eliminating the pole.
- **DERIVED CONSTRAINED/DEGENERATE STAR MAP:** stable rank below `12`; report
  all null singular directions without calling them gauge unless an actual
  symmetry/constraint proof is supplied.
- **INCONCLUSIVE:** threshold or step-stability rule fails.
- **REFUTED ACTION BRANCH:** the complex action/gradient is not real, the
  Schlaefli control fails, or the old internal equation is not reproduced.

Even full rank does not prove a global evolution. This protocol fixes all
link edges, covers one local move, and does not test overlap scheduling,
constraint matching, first-class closure, continuum Lorentz invariance or
Planck units.

Only the targeted verifier and static registry guard may run. No full suite
and no PDF build.
