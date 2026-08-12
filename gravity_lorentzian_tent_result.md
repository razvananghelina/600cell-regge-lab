# The symmetric Lorentzian vacuum tent has no real stationary pole

Date: 2026-08-12

Preregistered protocol commit: `14dd9aa`

Registered verifier:
`reproducible/verify_gravity_lorentzian_tent.py`

Machine-readable result:
`reproducible/gravity_lorentzian_tent.json`

## Headline

> **DERIVED SCOPED LORENTZIAN VACUUM NO-GO.** On the canonical local
> 600-cell tent carrier, use ordinary real Lorentzian Regge calculus, a
> timelike tent pole, equilateral old/link edges, one arbitrary positive
> final-edge orbit `a'`, and zero volume/cosmological coefficient. Then the
> Regge equation has no nondegenerate real stationary tent-pole length.

This is stronger than failure of the static ansatz `a'=a`: it holds for every

```text
rho=tau^2/a^2>0,       q=a'^2/a^2>0,
```

including the physically stricter two-spacelike-boundary sector `q>1/3`.
The exact reason is an angular gap: every allowed internal dihedral angle is
at most `arccos(1/3)=70.5288... degrees`, while five congruent simplices need
`72 degrees` for zero deficit.

The targeted verifier passes `23/23`. No full suite was run.

This closes one candidate clock; it does not close Lorentzian dynamics in
general. A nonzero volume term can create a root, but only after supplying an
as-yet unselected dimensionless coefficient.

## 1. Complete scope

The carrier is the already certified join

```text
T_v=[v,v']*L_v,
```

where `L_v` is the icosahedral link of a 600-cell vertex. There are 12
internal triangles `[v,v',u]`, and five four-simplices meet at each one.

The result assumes:

1. old cone and link edges have squared length `a^2`;
2. new cone edges have squared length `q a^2`, with `q>0` fixed;
3. the pole has squared interval `-rho a^2`, with `rho>0` varied;
4. every four-simplex has signature `(-,+,+,+)`;
5. timelike hinges use their real area magnitude and the ordinary rotation
   angle in their Euclidean normal plane;
6. the Lorentzian Regge internal-edge equation follows from the Schlaefli
   identity;
7. the volume/cosmological coefficient, matter and higher-curvature terms
   vanish.

The Regge dynamics and causal-angle convention are external inputs, not
selected outputs of the theory. Lorentzian simplex realizability is discussed
by [Tate--Visser](https://arxiv.org/abs/1110.5694); the Regge action with
timelike triangles is treated in
[Liu--Han](https://arxiv.org/abs/1810.09042), and the canonical role of tent
moves in [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974).

## 2. Exact Lorentzian simplex

Use one old vertex as origin. Let `y=v'-v` and let `u_1,u_2,u_3` be a link
triangle. In units `a=1`, their scalar products are

```text
y.y=-rho,
u_i.u_i=1,
u_i.u_j=1/2,
y.u_i=(1-q-rho)/2.
```

The signed Gram determinant factorizes as

```text
det G = -[3(q+rho-1)^2+8 rho]/16.
```

It is strictly negative for all `q,rho>0`. The link `3 x 3` block has
eigenvalues `(2,1/2,1/2)`, so the Schur complement proves that the complete
Gram matrix has exactly three positive and one negative eigenvalue. Unlike
the Euclidean tent, this Lorentzian simplex has no upper proper-time bound.

There is a framing correction worth preserving. Pairwise spacelike new
edges (`q>0`) do not alone make the final tetrahedral star spacelike. Its
three-Gram eigenvalues are

```text
(1/2,1/2,3q-1),
```

so a genuine second spacelike boundary requires `q>1/3`. The no-go was proved
on the larger domain `q>0`, hence it remains valid after this restriction.

## 3. Timelike hinge and exact angle gap

The internal triangle is timelike. Its real area magnitude is

```text
A_t = (a^2/4) sqrt[4rho+(1-q-rho)^2].
```

At fixed boundary data,

```text
dA_t/drho = a^2(q+rho+1)
             / [4 sqrt(4rho+(1-q-rho)^2)] > 0.
```

The facets meeting along this hinge have spacelike normals. Inverting the
signed Gram matrix gives the ordinary dihedral rotation angle

```text
cos(theta)=
[q^2+2q rho-2q+rho^2+1]
/
[2(q^2+2q rho-2q+rho^2+rho+1)].
```

Put

```text
D=(q+rho-1)^2+3rho.
```

Then

```text
cos(theta)=1/2-rho/(2D),
1/3 <= cos(theta) < 1/2.
```

The comparison with the fivefold flat angle is exact, not numerical:

```text
cos(theta)-cos(2*pi/5)
 = (3-sqrt(5))[(q+rho-1)^2+rho/phi^2]/(4D)
 > 0.
```

Since cosine decreases on `(0,pi)`, this proves

```text
pi/3 < theta <= arccos(1/3) < 2*pi/5.
```

Therefore the common fivefold deficit obeys the uniform bound

```text
epsilon=2*pi-5theta
       >=2*pi-5 arccos(1/3)
       =0.128388220475713... rad > 0.
```

The lower bound is attained inside the algebraic family when
`q+rho=1`; part of that locus also obeys `q>1/3`. The gap is not a numerical
finite-domain artifact.

## 4. Independent Minkowski reconstruction

The symbolic angle was attacked without analytically continuing the
Euclidean formula. At five deterministic static and nonstatic points, the
verifier:

1. constructs explicit edge vectors in `R^(1,3)`;
2. reproduces every prescribed signed scalar product;
3. obtains the two facet normals as numerical null vectors;
4. checks both normals have squared norm `+1`;
5. computes the interior angle directly from their Minkowski product.

The maximum edge-Gram residual is `5.6e-16`; the maximum angle residual
against the exact formula is `6.7e-16`. Every sampled Gram matrix has exactly
one negative eigenvalue. This is an independent branch/sign control, not the
sole evidence for the all-parameter theorem.

## 5. The pole equation

The standard Lorentzian Schlaefli reduction gives

```text
dS_R/dtau = 12 epsilon dA_t/dtau.
```

All 12 area derivatives and deficits have the same strictly positive sign.
Consequently there is no cancellation and no stationary pole on
`q,rho>0` at zero volume coefficient.

This conclusion is insensitive to an overall nonzero action sign or the
conventional real factor multiplying timelike-hinge terms: such a factor can
reverse or rescale the equation but cannot create a zero. It would fail if a
different legitimate oriented-angle branch supplied cancellations; the
explicit convex Minkowski normals and the standard timelike-hinge rotation
branch found none.

## 6. What happened to the Euclidean golden root

On the static Lorentzian slice `q=1`, the exact formulae are

```text
det G = -rho(8+3rho)/16,
cos(theta)=(2+rho)/[2(3+rho)],
A_t=(a^2/4)sqrt[rho(4+rho)],
V_4=(a^4/96)sqrt[rho(8+3rho)].
```

At the Euclidean golden value `rho=phi^-2`, the Lorentzian deficit is

```text
epsilon_L=0.228577374383... rad,
```

not zero. Thus the conditional Euclidean result `t/a=phi^-1` is a
positive-definite reflection geometry, not a physical timelike solution of
the corresponding vacuum Regge equation.

This also sharpens the relation to `a1=5`: incidence five is genuinely
load-bearing, but here it creates the obstruction. Five wedges demand a
larger angle than this Lorentzian simplex family can provide. The arithmetic
symbol `a1=5` still does not supply time, `c` or a unit scale.

## 7. Volume-term hostile control

Add a generic term

```text
-lambda V_total,       V_total=20 V_4,
```

and define `ell=lambda a^2`. On the static branch, stationarity at a chosen
`rho` would require

```text
ell(rho)=(72/5) epsilon(rho)
         (rho+2)/(3rho+4)
         sqrt[(3rho+8)/(rho+4)].
```

For example,

```text
ell(0+)      = 1.307292211049...,
ell(1)       = 3.214668192665...,
ell(infinity)= 8sqrt(3)pi/5
             = 8.706236948324....
```

So a nonzero coefficient can manufacture a stationary pole, but the desired
pole fixes the coefficient rather than the geometry selecting the pole. The
repository currently has no derived transfer from its finite spectral
moments to this dimensionful Regge normalization. Comparing one value of
`ell` with a numerically convenient repository ratio would be fitting.

No monotonicity or uniqueness claim for the full coefficient curve is needed
or made here.

## 8. Status ledger

| Claim | Status |
|---|---|
| The symmetric signed Gram matrix is Lorentzian for every `q,rho>0` | **DERIVED** |
| A spacelike final tetrahedral boundary requires `q>1/3` | **DERIVED** |
| Every internal tent hinge is timelike | **DERIVED** |
| Its facet-normal plane is Euclidean and the displayed angle is real | **DERIVED** |
| The angle is always strictly below `2*pi/5` | **DERIVED** |
| The fivefold deficit is uniformly positive | **DERIVED** |
| A zero-volume symmetric Lorentzian vacuum pole exists | **REFUTED** |
| Releasing the symmetric final radius `a'` repairs the vacuum equation | **REFUTED** |
| The Euclidean golden pole is a Lorentzian vacuum tick | **REFUTED** |
| Ordinary Lorentzian Regge calculus is selected by the theory | **STRUCTURAL / OPEN** |
| A nonzero volume coefficient can create a chosen static root | **DERIVED CONDITIONAL** |
| The theory selects that coefficient or its normalization | **OPEN** |
| Nonsymmetric data, matter or higher curvature can never repair the route | **NOT CLAIMED** |
| A global non-overlapping causal tent schedule exists | **OPEN** |
| `c`, `G`, Planck time or Planck mass follows | **OPEN** |

## 9. Consequence and next gate

There is now a real dynamical negative, not merely a missing construction:

> The most canonical local candidate for turning one 600-cell frame into a
> Lorentzian next frame is not a vacuum solution, even after its symmetric
> final spatial radius is released.

The next step must not be to rename `a/phi` as a tick. The shortest honest
fork is:

1. determine whether an already derived term fixes the Regge volume
   coefficient and normalization without using a desired root;
2. if not, record that coefficient as absent rather than fit it;
3. then test nonsymmetric tent data or a different causal carrier, with its
   degrees of freedom preregistered before looking for a solution;
4. only after a real solution exists, compute its pre/post Legendre map,
   constraints, light cone and continuum/refinement behavior.

The geometry still supplies a three-dimensional spatial frame, a canonical
metric phase space and a canonical four-dimensional carrier. It does not yet
supply a selected physical time evolution.

The subsequent registered audit
`gravity_lorentzian_volume_selection_result.md` closes item 1 for the current
repository: the spatial spectral action fixes at most a favorable relative
sign, while its free cutoff spans every positive coefficient and has no
derived 3D-to-4D transfer. Thus it cannot select or evidence a tent root.

## 10. Reproduction history

The preliminary all-`q` no-go and formulae were disclosed in protocol commit
`14dd9aa`. The first implementation produced `20/22`: both failures were
SymPy representation artifacts (`sqrt(x^2)` retained as `Abs(x)` and a
structurally unequal exact radical expression). Equality after squaring plus
the independently proved positive branch, and exact subtraction instead of
structural equality, repaired the checks without changing any formula or
numerical result.

Hostile review then found that `q>0` does not by itself make the final
tetrahedral boundary spacelike; the correct condition is `q>1/3`. This was
added as a new check. Because the no-go holds on the larger `q>0` domain, its
physical conclusion was unchanged. The final targeted verifier passes
`23/23`.

No full suite and no PDF build were run.
