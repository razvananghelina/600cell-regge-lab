# The fivefold Regge tent equation has an exact golden static root

Date: 2026-08-12

Preregistered protocol commit: `749429e`

Registered verifier:
`reproducible/verify_gravity_tent_move_regge.py`

Machine-readable result:
`reproducible/gravity_tent_move_regge.json`

## Headline

> **DERIVED CONDITIONAL EUCLIDEAN GOLDEN TENT.** For the canonical local tent
> carrier over one 600-cell vertex, using the ordinary Euclidean
> zero-cosmological-constant Regge action and holding both spatial boundary
> stars equilateral with edge length `a`, the unique nondegenerate stationary
> tent-pole ratio is
>
> ```text
> t/a = phi^-1.
> ```

The factor five is not fitted: every internal tent triangle is incident to
exactly five four-simplices. Zero deficit therefore requires their common
dihedral angle to be `2*pi/5`, and the exact simplex Gram matrix then forces
the golden ratio.

> **DERIVED NON-SELECTION WARNING.** If the final spatial edge length `a'` is
> released, zero deficit is a continuous flat boundary-data family. The
> golden point is selected only by the additional static-return condition
> `a'=a`. A nonzero volume/cosmological coefficient shifts it immediately.

The targeted verifier passes `23/23`. No full suite was run.

The golden root was found exploratorily before preregistration and disclosed
in the protocol. It is a certified conditional derivation, not a blind
prediction and not a comparison with Planck units.

## 1. The exact local carrier

Choose any 600-cell vertex `v`. Its link is the icosahedral two-sphere

```text
f(L_v)=(12,30,20).
```

Introduce one new vertex `v'`. The tent region is the simplicial join

```text
T_v=[v,v']*L_v.
```

It consists of 20 four-simplices and has

```text
f(T_v)=(14,55,92,70,20),       chi(T_v)=1.
```

Its boundary is the initial cone star `v*L_v` and final cone star `v'*L_v`,
sharing the link at their corner. The internal hinges are the 12 triangles

```text
[v,v',u],       u a vertex of L_v.
```

Every link vertex occurs in five icosahedral triangles. Consequently every
internal hinge occurs in five congruent four-simplices. This is the geometric
origin of the integer five in the field equation.

The order-120 vertex stabilizer is transitive on all 12 link vertices, all 30
link edges and all 20 link triangles. Thus the symmetric new metric has two
length variables before specialization:

```text
a' = common length [v',u],
t  = length [v,v'].
```

The equality `a'=a` is not implied by the stabilizer; it is a static-boundary
condition.

## 2. Exact four-simplex geometry

Set

```text
r=t^2/a^2,       q=a'^2/a^2.
```

Under the static condition `q=1`, take the tent-pole vector and three link
edge vectors from `v` as a basis. Their Gram determinant is

```text
det G = a^8 r(8-3r)/16.
```

Hence the nondegenerate Euclidean domain is

```text
0<r<8/3.
```

The inverse-Gram normal formula gives the internal dihedral angle

```text
cos(theta(r))=(2-r)/(2(3-r)).
```

The internal hinge area and one four-simplex volume are

```text
A_int=(a^2/4)sqrt(r(4-r)),
V_4=(a^4/96)sqrt(r(8-3r)).
```

All formulae are exact algebraic identities, not numerical fits.

## 3. Why five produces the golden ratio

The common internal deficit is

```text
epsilon(r)=2*pi-5 theta(r).
```

The Regge equation for the internal tent pole follows from the Schlaefli
identity:

```text
dS_R/dt = 12 epsilon(r) dA_int/dt.
```

The verifier also constructs all internal and boundary dihedral angles of the
complete action directly. At three deterministic off-shell pole lengths, the
largest residual between the full finite-difference derivative and this
reduced expression is `2.6e-10`. Thus omitted boundary terms are not being
used to manufacture the root.

At a nondegenerate stationary point `dA_int/dt` is nonzero, so

```text
theta=2*pi/5,
cos(theta)=cos(72 degrees)=(sqrt(5)-1)/4=1/(2 phi).
```

Combining this with the Gram result gives

```text
(2-r)/(2(3-r))=1/(2 phi),
r=(3-sqrt(5))/2=phi^-2,
t/a=phi^-1.
```

The cosine is strictly decreasing throughout the Euclidean domain, so this
root is unique. The full action derivative at the root is `-1.8e-10` in the
independent numerical reconstruction.

More generally, if `n` congruent four-simplices meet at the hinge, the same
static equation would use `theta=2*pi/n`. What is specific here is the
600-cell incidence `n=5`; the golden value is not inserted separately.

## 4. Precise relation to the programme's `a1=5`

The arithmetic bootstrap's symbol `a1=5` and the local Regge incidence are not
the same theorem. The honest chain is

```text
arithmetic bootstrap a1=5
    -> chosen/structural 600-cell realization
    -> five tetrahedra around each spatial edge
    -> five 4-simplices around each corresponding tent hinge
    -> epsilon=2*pi-5 theta
    -> static Euclidean root t/a=phi^-1.
```

Therefore:

- **DERIVED within the 600-cell realization:** the same integer five enters a
  genuine variational gravity equation and produces `phi^-1`;
- **STRUCTURAL bridge:** identifying that incidence with the independently
  bootstrapped symbol `a1`;
- **NOT an independent derivation:** this calculation does not prove that
  `a1=5` uniquely forces the 600-cell;
- **NOT physical units:** it fixes only a ratio to the spatial edge length.

This is substantially stronger than merely observing another occurrence of
five, because removing the five changes the field equation. It is weaker than
a Planck-time prediction because the overall edge scale and Lorentzian clock
are absent.

## 5. The flat family: why this is not yet a selected tick

For unrestricted `q`, the exact Gram calculation gives

```text
cos(theta)=
(q^2-2qr-2q+r^2+1)
/
[2(q^2-2qr-2q+r^2-r+1)].
```

Let `y=t/a`. The zero-deficit equation has two algebraic branches

```text
q=1+y^2-y/phi,
q=1+y^2+y/phi.
```

The minus branch is the direct flat embedding obtained by moving `v'` along
the normal to the affine hyperplane of the 12 neighbours. It is a continuous
family of boundary metrics. Imposing the static condition `q=1` yields

```text
y(y-phi^-1)=0.
```

After excluding the degenerate `y=0`, the golden solution remains.

Thus the Regge equation alone does not select an elapsed time. It supplies a
flat vertex-displacement family; equality of the initial and final intrinsic
stars selects one nondegenerate return point on it. Interpreting the family as
a lapse/vertex-displacement gauge orbit is **STRUCTURAL**, consistent with
canonical flat Regge calculus, but a completed constraint algebra has not yet
been built here.

## 6. Independent geometric meaning of the root

At unit circumradius, choose the exact 600-cell vertex

```text
v=(-1,0,0,0).
```

Its 12 neighbours lie in the affine hyperplane

```text
x_0=-phi/2.
```

Reflection of `v` through that hyperplane gives

```text
v'=v/phi=(-phi^-1,0,0,0).
```

The old and new distances to all 12 neighbours are equal, and

```text
|v-v'| / |v-u| = phi^-1.
```

This independently reconstructs the same root without differentiating the
action. Geometrically, the tent region is the flat bipyramid between two
congruent tetrahedral stars on opposite sides of their common icosahedral
link.

The new point has radius `1/phi`, not one. It is not another vertex of the
original 600-cell on `S3`. Consequently this reflection does not define a
global rigid rotation from one 600-cell slice to another.

## 7. Hostile controls and physical cost

### Nonzero volume coefficient

The total volume of the 20 four-simplices has

```text
dV_total/dr = 0.367485838... a^4
```

at the golden point. Therefore adding any term

```text
-lambda V_total,       lambda != 0,
```

makes the action derivative nonzero there. The golden ratio depends on the
frozen `lambda=0` hypothesis. The repository has not derived a vanishing
cosmological/volume coefficient.

### Variational sign

For the positive action sign displayed in the protocol, the fixed-`q=1`
second derivative is negative:

```text
d^2 S_R/dr^2 = -1.583592135... .
```

The pole root is a maximum in that one Euclidean direction, not an energy
minimum. Reversing the conventional overall Euclidean action sign reverses
this label without changing the equation of motion.

### No synchronous global move

All 120 vertex tent carriers form one `H4` orbit, but each overlaps the 12
tent carriers at adjacent vertices. The complete orbit cannot be applied as a
simultaneous independent layer. A selected colouring/schedule or a sum over
moves remains necessary.

### Not Lorentzian

The derivation uses a positive-definite Gram matrix and ordinary Euclidean
angles. Replacing the pole by a timelike interval is not the substitution
`t^2 -> -t^2` inside this proof: Lorentzian simplex admissibility, boost angles
and action branches must be rebuilt. No value of `c` follows.

## 8. Status ledger

| Claim | Status |
|---|---|
| The local carrier is `[v,v']*Icos` with 20 four-simplices | **DERIVED** |
| Every internal hinge has incidence five | **DERIVED** |
| The full Regge tent equation reduces to `epsilon=0` | **DERIVED** |
| Under Euclidean, `Lambda=0`, `a'=a`, the unique ratio is `phi^-1` | **DERIVED CONDITIONAL** |
| The same root is the neighbour-hyperplane reflection | **DERIVED** |
| The integer five is load-bearing in the equation | **DERIVED** |
| This independently re-derives the arithmetic bootstrap `a1=5` | **REFUTED** |
| Zero deficit alone selects a unique elapsed time | **REFUTED** |
| The unrestricted zero-deficit boundary data form a flat family | **DERIVED** |
| The flat family is a physical gauge orbit of a completed theory | **STRUCTURAL / OPEN** |
| The golden root survives a nonzero volume term | **REFUTED** |
| The frozen-sign root is a minimum | **REFUTED; IT IS A MAXIMUM** |
| All 120 tent moves form a synchronous invariant tick | **REFUTED** |
| The pole is Lorentzian time or Planck time | **OPEN / NOT CLAIMED** |
| `c`, `G`, Planck mass or a graviton follows | **OPEN** |

## 9. What comes next

This result changes the next question. It is no longer “can five enter a
dynamical equation?” It can. The next falsifiable gates are:

1. derive the Lorentzian version of the same icosahedral tent carrier;
2. enumerate its admissible boost-angle branches before solving them;
3. include the theory's actual volume/higher-spectral terms rather than set
   their coefficients to a desired value;
4. compute the full pre/post Legendre map in `(a',t)` and classify the flat
   vertex-displacement direction;
5. find a target-independent orbit-complete schedule of non-overlapping tent
   moves;
6. only afterward test a scale-independent light cone or Planck units.

The most dangerous shortcut would be to call `a/phi` the Planck time after
setting `a` equal to the Planck length. That would assume both the missing
scale and the conversion speed. The certified result is the dimensionless
Euclidean ratio only.

## 10. Literature/novelty boundary

Tent moves and their canonical role are established in
[Dittrich--Hoehn](https://arxiv.org/abs/1108.1974). Targeted searches did not
locate a primary source stating this exact 600-cell/icosahedral
`t/a=phi^-1` calculation. Search absence is weak evidence, so no claim of
literature novelty is made here.

## 11. Reproduction history

The first targeted verifier passed `22/22`. Hostile review then released the
final spatial length and found the exact two-branch zero-deficit family. This
did not change the preregistered conditional root; it changed its physical
interpretation from a candidate selected tick to the static intersection of
a flat boundary-data family. The explicit control was added. The final
targeted verifier passes `23/23`.

No full suite and no PDF build were run.
