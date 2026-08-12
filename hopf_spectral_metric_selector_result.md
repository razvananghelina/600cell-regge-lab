# What the present spectral data do and do not select

Date: 2026-08-12  
Protocol commit: `a946935`  
Registered verifier: `reproducible/verify_hopf_spectral_metric_selector.py`  
Machine-readable result: `reproducible/hopf_spectral_metric_selector.json`

## Headline

There is a genuine but conditional advance:

> **DERIVED CONDITIONAL SHAPE SELECTION.**  Among all smooth left-invariant
> metrics on `S3=SU(2)` with fixed volume, the ordinary full-de Rham heat
> coefficient `A2` has the round metric as its unique global minimum.

This strengthens the previous local Hessian result to a global theorem on the
complete homogeneous metric cone.  It says that the already derived smooth
curvature term drives every homogeneous anisotropy toward the round shape.

It is not yet a complete emergent-gravity action.  The finite theory still
does not evaluate one all-form spectral functional on both the smooth round
metric and the singular fixed-Regge metric.  Nor do the finite moments or the
number `31/11` select an absolute length, heat time, Newton constant or Planck
scale.

The targeted verifier passes `24/24`.  No full-suite run was requested or
performed.

## 1. Exact global round selection on the smooth branch

For a positive left-invariant metric Gram matrix `G`, the independently
derived Koszul formula is

```text
R(G)=2 ((Tr G)^2-2 Tr(G^2))/det(G).
```

Fix `det(G)=1`, and let its eigenvalues be `a,b,c>0`, so `abc=1`.  Put

```text
p=a+b+c,   q=ab+bc+ca.
```

Then

```text
R(G)/2=4q-p^2.
```

Schur's degree-three inequality and AM--GM give

```text
p^3+9 >= 4pq,   p>=3,
R(G)/2 <= 9/p <= 3.
```

Equality in both steps forces `a=b=c=1`.  Hence

```text
R(G)<=6,
```

with unique equality at the round metric.  The verifier supplies an exact
positivity certificate for Schur: after ordering `a>=b>=c` and writing

```text
a=z+y+x, b=z+y, c=z,   x,y>=0, z>0,
```

its deficit is

```text
x^3+2x^2y+x^2z+xyz+y^2z >= 0.
```

For the ordinary heat trace of the full exterior de Rham operator in three
dimensions,

```text
A2=-(2/3) integral R dVol.
```

At fixed volume this therefore obeys

```text
A2 >= -4 Vol,
```

with the round metric as the unique homogeneous minimizer.  The graded/index
supertrace cancels this coefficient and supplies no analogous selector.

This is an action-like shape-selection result, not merely an equality of
representations.  Its complete hypotheses are load-bearing: smooth metric,
left invariance, fixed volume, ordinary rather than graded heat trace, and use
of the `A2` coefficient rather than a fully specified finite-cutoff action.

## 2. Why `31/11` does not select the missing scale

The fixed 2640-state operator has exact finite moments

```text
c0=2640, c1=14880, c2=55920,
r=c1/(2c0)=31/11.
```

The number `r` is an exact spectral invariant of that normalized matrix, but
it has inverse-length-squared units.  Under a constant metric rescaling,

```text
g -> c^2 g,        D_g^2 -> c^-2 D_g^2,
r -> c^-2 r,
K_(c^2 g)(t)=K_g(t/c^2).
```

Thus `31/11` is not an absolute metric or Planck scale.  Its reciprocal can be
used to form a scale-covariant heat time,

```text
t=alpha/r,
```

but every positive dimensionless `alpha` is equally covariant.  The active
paper correctly says that the continuum interpretation of `r` remains open;
it does not derive `alpha=1` or identify `r` with a physical cutoff.

There is still a useful opening here: `r` can provide an internally normalized
spectral unit once an action supplies `alpha`.  Calling that unit the Planck
scale before deriving the action would be circular.

## 3. A heat trace without its time does not order geometries

The exact preregistered control spectra

```text
A={0,1,10},   B={0,2,3}
```

have heat-trace difference, with `x=exp(-t)`,

```text
K_A-K_B=x+x^10-x^2-x^3.
```

At `x=100/101` the difference is exactly negative, while at `x=1/2` it is
exactly `129/1024>0`.  The preferred spectrum reverses without changing
anything except heat time.  This is a general counterexample to treating
`Tr exp(-tD^2)` with unspecified `t` as a variational ordering.

Overall scale also runs to a boundary.  For positive `k`,

```text
Tr((D_(c^2g)^2)^k)=c^(-2k) Tr((D_g^2)^k)
```

is strictly decreasing with `c`, while a positive heat trace is strictly
increasing with `c`.  Neither has an interior scale extremum.  This does not
obstruct the fixed-volume shape theorem in Section 1; it separates shape
selection from scale selection.

## 4. Why round versus Regge is still not one spectral comparison

The active implementations cover different domains:

| Construction | Full exterior algebra | Variable metric | Metric domain | Function/scale selected |
|---|---:|---:|---|---:|
| 2640-state incidence `D` | yes | no | fixed Euclidean cochain adjoint | no |
| smooth homogeneous `D_g` coefficient | yes, at `A2` | yes | smooth left-invariant `SU(2)` | coefficient only |
| exact Whitney tower | yes | no in the certified spectrum | fixed Regge metric | no |
| projected barycentric diagnostic | no, scalar `P1` only | embedding-dependent | projected smooth mesh | no |

The endpoints also belong to different analytic categories.  The round
metric is smooth with distributed scalar curvature.  The Regge metric is flat
inside each tetrahedron and has curvature on a lower-dimensional skeleton.
Heat expansions on conical singularities generally acquire additional
surface/singularity contributions; this is the subject of Fursaev's
spectral-geometry analysis, the spin-dependent cone calculation of Fursaev
and Miele, and Vassilevich's heat-kernel review:

- <https://arxiv.org/abs/hep-th/9405143>
- <https://arxiv.org/abs/hep-th/9605153>
- <https://arxiv.org/abs/hep-th/0306138>

Therefore the smooth identity `A2=-(2/3) integral R` cannot simply be applied
to the Regge deficit sum.  The repository contains neither the required
singular full-form coefficient nor a verified common `D_n(u)` spectrum.
Reporting a numerical preferred `u` now would compare different operators,
not minimize one defined action.

## 5. Attack on the decision criterion

Demanding the same minimizer for every heat time is stronger than ordinary
spectral-action theory.  A physical theory is allowed to specify a particular
cutoff function and scale.  The problem here is narrower: this repository has
not specified them.  Consequently cutoff independence is an appropriate
robustness gate for a claim of selection **from the existing data**, but it is
not a theorem that every future spectral action must pass.

Likewise, restricting to homogeneous metrics is physically severe.  It proves
a global result for the five Hopf anisotropy directions and their nonlinear
completion, not for arbitrary local metric fields.  It supplies no spatial
gradient term, Hamiltonian constraint or diffeomorphism quotient.

## 6. Status ledger

| Claim | Status |
|---|---|
| `R(G)<=6` at fixed volume on smooth left-invariant `SU(2)` | **DERIVED** |
| Equality is unique at the round metric | **DERIVED** |
| Ordinary de Rham `A2` uniquely minimizes at round shape in that class | **DERIVED CONDITIONAL SHAPE SELECTION** |
| The old positive Hopf Hessian extends to a global homogeneous theorem | **DERIVED** |
| The index/graded trace selects the same metric | **REFUTED** |
| `31/11` is an invariant absolute cutoff or Planck scale | **REFUTED** |
| `31/11` can define an internal inverse-length unit | **STRUCTURAL POSITIVE** |
| That unit fixes the dimensionless heat parameter `alpha` | **REFUTED** |
| An unspecified heat trace orders all candidate metrics | **REFUTED** |
| Existing positive moments select an interior overall scale | **DERIVED NEGATIVE** |
| One complete all-form action compares `g_R` and `g_0` | **OPEN / absent** |
| Smooth `A2` applies unchanged to the Regge endpoint | **NOT ESTABLISHED; analytically unsafe** |
| The theory now derives Newton's constant or the Planck scale | **OPEN** |
| The theory now derives a propagating Lorentzian graviton | **OPEN** |
| Emergent gravity is impossible | **NOT CLAIMED** |

## 7. Documentation consistency found during the audit

The source audit exposed stale statements in `one_integer_paper_v5.tex` and a
few legacy passages in `one_integer_paper.tex`: they still called the
alternating `Box` nullity `-4` a derived spacetime dimension and then used that
reading to label the `3/4` Morrey exponent and the Planck hierarchy as derived.
Those statements contradicted the already binding dimension audit elsewhere
in the same files.

The contradictions have been removed rather than covered by a later caveat:

- `-4` is an alternating count of a non-complex operator hierarchy;
- the actual finite/cochain geometry is three-dimensional;
- finite summability does not select `p=4`;
- `Tr(phi^3)=4` and the weight `4 phi^2` remain exact algebraic facts;
- their spacetime, Morrey and Planck readings are **PATTERN/OPEN**.

No PDF was built.

## 8. Consequence and next gate

### Result of the singular-coefficient gate

The first option below has now been completed in
`regge_de_rham_cone_selector_result.md`.  Cheeger's skeleton-local expansion
and the exact scalar/vector Hodge cone coefficients give a common ordinary
all-form `A2` at the fixed Regge endpoint.  At equal volume the round endpoint
is lower by `0.0848366160...`; using only the smooth Regge-curvature limit
would reverse the sign and is therefore explicitly refuted.  The subsequently
preregistered full-path audit includes the mandatory de Rham transmittal face
term and finds a converged monotone preference for round at all 201 frozen
grid points.  That interior result is **PATTERN**, not a continuum theorem,
and the calculation remains a single-coefficient result rather than a
complete metric action.

The most accurate statement at this stage was:

> The ordinary geometric operator contains a genuine restoring functional for
> all homogeneous spin-two shape deformations, and its unique smooth
> fixed-volume vacuum is round.  The theory has not yet selected the analytic
> branch, absolute scale or Lorentzian dynamics.

The later interval certificate `round_regge_a2_interval_result.md` proves the
full continuous affine-path inequality. The later sign audit
`round_regge_spectral_action_sign_result.md` also proves that every standard
positive cutoff gives this `A2` the favorable asymptotic sign. The magnitude,
finite-cutoff dominance, transverse metric Hessian, Lorentzian time, local
propagation, universal stress-energy coupling and Newton/Planck normalization
remain independent gates.
