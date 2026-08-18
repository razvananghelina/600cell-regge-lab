# Preregistered protocol: Lorentzian 600-cell Regge tent gate

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- PRELIMINARY NO-GO DISCLOSED**

## 1. Question and complete hypotheses

Let `K` be the fixed equilateral 600-cell boundary with spatial edge length
`a>0`. Choose a vertex `v`, write `L_v` for its icosahedral link, introduce a
new vertex `v'`, and use the same canonical local tent carrier

```text
T_v = [v,v'] * L_v.
```

Freeze the following hypotheses before certification:

1. every old cone edge `[v,u]` and every link edge has spacelike squared
   length `a^2`;
2. every new cone edge `[v',u]` has spacelike squared length
   `a'^2=q a^2`, with `q>0` held fixed under the pole variation;
3. the tent pole has timelike squared interval
   `s[v,v']=-tau^2=-rho a^2`, with `rho>0`;
4. each four-simplex has Lorentzian signature `(-,+,+,+)`;
5. use ordinary real Lorentzian Regge calculus on the timelike internal
   triangular hinges: their normal planes are Euclidean, their deficit is

```text
epsilon = 2*pi - sum theta,
```

   and the internal-edge equation follows from the Lorentzian Schlaefli
   identity;
6. first set the volume/cosmological coefficient to zero;
7. vary only `tau`, holding all boundary squared lengths fixed;
8. include no matter, higher-curvature or spectral-cutoff term.

The same incidence already certified in `verify_gravity_tent_move_regge.py`
puts five congruent four-simplices around each of the 12 timelike internal
hinges `[v,v',u]`.

The Lorentzian Regge convention is external methodology, not a result of the
600-cell programme. The causal distinction used here is consistent with the
Lorentzian simplex analysis of Tate--Visser
([arXiv:1110.5694](https://arxiv.org/abs/1110.5694)); timelike triangles in a
Lorentzian Regge action are also treated in
[arXiv:1810.09042](https://arxiv.org/abs/1810.09042). Canonical tent moves are
reviewed in [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).

No claim is made that these hypotheses, `q`, or a volume coefficient are
selected by the theory.

## 2. Preliminary result and provenance

Before registration, exploratory symbolic algebra found a stronger result
than the static `q=1` test:

> For every real `rho>0` and `q>0`, the internal dihedral angle is strictly
> smaller than `2*pi/5`. Hence the fivefold deficit is strictly positive and
> the zero-volume-coefficient pole equation has no stationary point.

This negative is disclosed, so the verifier is a hostile certification, not
a blind search. The Euclidean golden root `tau/a=phi^-1` does not continue to
this real timelike branch.

## 3. Frozen exact checks

Use `v` as the origin of one four-simplex. Let `y=v'-v` and let `u_1,u_2,u_3`
be the three link vertices of its link triangle. In units `a=1`, derive the
signed Gram matrix from

```text
y.y = -rho,
u_i.u_i = 1,
u_i.u_j = 1/2,
y.u_i = (1-q-rho)/2.
```

The verifier must establish, rather than insert:

1. the exact determinant

```text
det G = -[3(q+rho-1)^2+8 rho]/16 < 0;
```

2. because the link `3 x 3` block is positive definite, `G` has inertia
   `(3 positive, 1 negative)` for every `q,rho>0`;
3. the internal triangle is timelike and has real magnitude

```text
A_t = (a^2/4) sqrt(4 rho + (1-q-rho)^2);
```

4. `dA_t/dtau>0` throughout the domain;
5. the two facets meeting at that timelike triangle have spacelike normals;
6. both the inverse-Gram normal formula and an independently constructed
   Minkowski-coordinate normal calculation give

```text
cos(theta) =
[q^2+2q rho-2q+rho^2+1]
/
[2(q^2+2q rho-2q+rho^2+rho+1)];
```

7. with

```text
D=(q+rho-1)^2+3 rho,
```

   this becomes `cos(theta)=1/2-rho/(2D)`, hence
   `1/3 <= cos(theta) < 1/2` and

```text
pi/3 < theta <= arccos(1/3) < 2*pi/5;
```

8. equivalently, subtraction of `cos(2*pi/5)` has a strictly positive
   numerator proportional to

```text
(q+rho-1)^2 + phi^(-2) rho;
```

9. therefore `epsilon=2*pi-5 theta>0` for all `q,rho>0`;
10. the standard Lorentzian Regge pole equation

```text
dS_R/dtau = 12 epsilon dA_t/dtau
```

   has no real nondegenerate stationary solution.

The Minkowski-coordinate control must be evaluated at deterministic static
and nonstatic points. It may not obtain the angle solely by substituting
`t -> i tau` in the Euclidean formula.

## 4. Static branch and volume-term control

For `q=1`, independently reduce the formulae to

```text
det G = -rho(8+3rho)/16,
cos(theta)=(2+rho)/[2(3+rho)],
A_t=(a^2/4)sqrt[rho(4+rho)],
V_4=(a^4/96)sqrt[rho(8+3rho)].
```

If a generic real term `-lambda V_total` is added, call
`ell=lambda a^2`. The verifier must derive the coefficient required to make
a chosen `rho` stationary:

```text
ell(rho) = (72/5) epsilon(rho)
           (rho+2)/(3rho+4)
           sqrt[(3rho+8)/(rho+4)].
```

This is not a prediction of `lambda`: it demonstrates that a nonzero volume
coefficient can manufacture a root and that the root then depends on an
unselected dimensionless input. Check the exact endpoint limits

```text
ell(0+) = (36 sqrt(2)/5)[2*pi-5 arccos(1/3)],
ell(infinity) = 8 sqrt(3) pi/5.
```

Do not identify `lambda` with a measured cosmological constant without a
separate normalization derivation.

## 5. Decision boundary and labels

- **DERIVED SCOPED LORENTZIAN VACUUM NO-GO:** all exact signature, causal
  angle, area-variation and fivefold-deficit checks pass under hypotheses
  1--8. Then neither the static ansatz nor any real symmetric `q>0` boundary
  data admits a zero-volume-coefficient Lorentzian tent solution.
- **REFUTED:** a valid Lorentzian simplex/angle branch reaches
  `theta=2*pi/5`, or the pole equation contains a cancellation omitted here.
- **OPEN/INCOMPLETE:** causal angle signs or the Lorentzian Schlaefli
  reduction cannot be reconciled.

Even a certified no-go is local and symmetry-reduced. It rules out this
specific vacuum tent clock; it does not rule out nonsymmetric tent data,
matter, a derived volume term, a different causal carrier, or a dynamical
fourth direction in general.

Only the targeted verifier and the static registry guard may run. No full
suite and no PDF build.
