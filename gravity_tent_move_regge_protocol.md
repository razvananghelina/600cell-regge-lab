# Preregistered protocol: the symmetric 600-cell Regge tent move

Date: 2026-08-12

Status at registration: **PROTOCOL ONLY -- PRELIMINARY GOLDEN ROOT DISCLOSED**

## 1. Complete hypotheses

Let `K` be the fixed equilateral 600-cell boundary with spatial edge length
`a>0`. Choose one vertex `v`; all choices are equivalent under `H4`. Let
`L_v` be its simplicial link and introduce one new vertex `v'`. The local
four-dimensional tent carrier is

```text
T_v = [v,v'] * L_v.
```

Freeze the following metric/action hypotheses before certification:

1. the calculation is **Euclidean**;
2. every old spatial edge has length `a`;
3. every new final edge `[v',u]`, for `u` in `L_v`, also has length `a`
   (a static equilateral boundary ansatz);
4. the tent pole `[v,v']` has variable positive length `t`;
5. use the ordinary zero-cosmological-constant Regge action with its standard
   boundary exterior-angle term,

```text
S_R = sum_internal_h A_h (2*pi-sum theta_h)
    + sum_boundary_h A_h (pi-sum theta_h);
```

6. vary only the internal tent-pole length while holding all boundary lengths
   fixed;
7. no matter term, higher-curvature term or spectral cutoff term is included.

These hypotheses are not claimed to have been selected by the theory. In
particular, `Lambda=0`, equality of old/new spatial lengths and Euclidean
signature are load-bearing.

The external canonical-Regge benchmark is
[Dittrich--Hoehn](https://arxiv.org/abs/1108.1974). The use of Regge action as
Hamilton's principal function is external methodology, not a prior result of
the 600-cell programme.

## 2. Preliminary result and provenance

Before this registration, an exploratory calculation found the candidate

```text
r=t^2/a^2=phi^(-2),      t/a=phi^(-1).
```

It also noticed that this is the Euclidean reflection of `v` across the affine
hyperplane of its 12 neighbours, sending a unit-circumradius vertex to
`v/phi`. Consequently the golden value is **not blind evidence** and must not
be compared with Planck time, `c` or any measured target. The purpose of the
registered audit is to certify or refute the complete variational chain and
its limitations.

## 3. Frozen combinatorial checks

The verifier must independently rebuild the complex and establish:

1. `L_v` is an icosahedral 2-sphere with `f=(12,30,20)`;
2. every link vertex lies in five link triangles;
3. `T_v` has exact `f=(14,55,92,70,20)` and Euler characteristic one;
4. its boundary consists of the two cone stars `v*L_v` and `v'*L_v`, sharing
   `L_v`;
5. it has 12 internal triangular hinges `[v,v',u]`, each incident to exactly
   five congruent four-simplices;
6. the vertex stabilizer makes the 12 old edges one orbit, the 12 new edges
   one orbit, and the tent pole a singleton. Thus the symmetric new metric
   has two variables `(a',t)` before the frozen specialization `a'=a`.

## 4. Frozen exact geometry and action checks

Put `r=t^2/a^2`. The verifier must derive, not assume:

1. the Gram determinant and nondegenerate Euclidean domain

```text
det G = a^8 r(8-3r)/16,       0<r<8/3;
```

2. the internal-hinge dihedral cosine

```text
cos(theta(r))=(2-r)/(2(3-r));
```

3. the internal area and one-simplex four-volume

```text
A_int=(a^2/4)sqrt(r(4-r)),
V_4=(a^4/96)sqrt(r(8-3r));
```

4. the common deficit

```text
epsilon(r)=2*pi-5*theta(r);
```

5. the Schlaefli-reduced tent equation

```text
dS_R/dt = 12 epsilon(r) dA_int/dt;
```

with a direct independent reconstruction of all internal and boundary angles
at several deterministic points;

6. `dA_int/dt` is nonzero at every claimed root;
7. the root is unique in the full Euclidean domain;
8. the exact root is `r=(3-sqrt(5))/2=phi^-2`;
9. with the displayed overall action sign, classify the second variation in
   the tent-pole direction;
10. independently reflect one exact 600-cell vertex through its neighbour
    hyperplane and compare the resulting length ratio.

## 5. Mandatory hostile controls

The machine certificate must also show:

1. adding a volume term `-lambda sum V_4` with any nonzero `lambda` gives a
   nonzero derivative at the golden zero-deficit point, because
   `d(sum V_4)/dt != 0` there;
2. allowing `a' != a` restores a genuine two-parameter symmetric metric
   family, so the golden root is not a general evolution law;
3. analytic continuation to a timelike tent pole is not covered by the
   positive-definite Gram/ordinary-angle calculation;
4. the local tent moves at all 120 vertices overlap, so the golden local move
   is not automatically a global synchronous tick.

## 6. Decision boundary

- **DERIVED CONDITIONAL EUCLIDEAN GOLDEN TENT:** all combinatorial, Gram,
  angle, full-action and uniqueness checks pass under hypotheses 1--7.
- **REFUTED GOLDEN TENT:** the full boundary Regge action does not reduce to
  the claimed tent equation or the unique stationary ratio differs.
- **OPEN/INCOMPLETE:** numerical angle conventions, boundary terms or exact
  identities cannot be reconciled.

Even a positive result is not physical time. It is a stationary Euclidean
flat-completion ratio for one local symmetric move. Lorentzian signature,
the cosmological/volume coefficient, an orbit-complete schedule, the
Hamiltonian constraint class, `c`, `G`, Planck time and Planck mass remain
open.

Only the targeted verifier and a static registry check may run. No full suite
and no PDF build.
