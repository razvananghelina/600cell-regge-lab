# Prior-art gate: conserved-mass homothetic 600-cell slab

Date: 2026-08-16

Status: **completed before evaluating the non-static candidate in the
repository action**.

Upstream exact static-family result: `f6b450a`.

This is a targeted primary-source map, not proof of external novelty.

## 1. Exact object and complete hypotheses

Use each of the two already derived order-24 staircase schedules on the
product of two regular 600-cell boundaries.  Retain the complete Lorentzian
Regge curvature action, De Felice--Fabri point-dust action, angle branch,
triangle constants, orbit multiplicities and logarithmic derivative
conventions certified upstream.

Let the lower and upper regular edges be `L_- > 0` and `L_+ > 0`.  Realize
the unit 600-cell vertices as vectors `u_i` in `R^4`, with

```text
u_i . u_i = 1,
u_i . u_j = phi/2 for adjacent vertices,
phi = (1+sqrt(5))/2.
```

A regular 600-cell of edge `L` therefore has embedding circumradius
`R_circ=phi L`.  Put its two homothetic copies at

```text
X_i^- = (0, phi L_- u_i),
X_i^+ = (Delta v, phi L_+ u_i)
```

in ambient signature `(-++++)`.  Define `tau>0` to be the proper length of
each same-vertex strut, so

```text
Delta v^2 = tau^2 + phi^2 (L_+-L_-)^2.
```

It follows without an adjustable cross-edge parameter that the intrinsic
squared lengths on the staircase carrier are

```text
q_old       = L_-^2,
q_new       = L_+^2,
pole square = -tau^2,
diagonal    = L_- L_+ - tau^2.
```

The calculation will independently rederive and mechanically check these
relations from the exact dot products.  It will require `tau^2<L_- L_+` and
will separately certify the Lorentzian inertia and complex-angle branch of
every simplex orbit; positivity of the displayed diagonal alone is not
assumed to be a sufficient branch proof.

Fix the total dust mass throughout the slab at the upstream static value

```text
epsilon3 = 2*pi-5*acos(1/3),
M0       = (90/pi)*epsilon3*L0,
```

and take `L_-=L0`.  The dust action remains the sum of the 120 equal particle
worldline lengths and does not replace `M0` by a function of `L_+`.
There is no cosmological constant.

Three inequivalent equations must not be conflated:

1. the **global homothetic lapse equation**, obtained after substituting the
   geometric ansatz and varying `tau`;
2. the **complete local internal system**, the 35 independent staircase and
   pole orbit derivatives before the homothetic identifications;
3. the **canonical junction equation** equating the new slab's pre-momentum
   to the preceding slab's post-momentum.

The first can hold while the second fails.  Neither one by itself supplies
the third.

## 2. KNOWN

Collins--Williams-type Regge cosmology, regular-polytopal spatial slices and
homogeneous struts are old.  The later simplicial Sorkin construction was
applied explicitly to a dust-filled 600-cell:

- Barrett et al., *A Parallelizable Implicit Evolution Scheme for Regge
  Calculus*, <https://arxiv.org/abs/gr-qc/9411008>;
- De Felice and Fabri, *The Friedmann universe of dust by Regge Calculus:
  study of its ending point*, <https://arxiv.org/abs/gr-qc/0009093>;
- De Felice and Fabri, *Singularities of the closed RW metric in Regge
  Calculus: a generalized evolution of the 600-cell*,
  <https://arxiv.org/abs/gr-qc/0106077>.

De Felice and Fabri use the same five-dimensional Lorentzian embedding of
homothetic spherical sections.  They keep the total mass fixed, choose the
timelike strut length as a lapse condition, omit the associated equation as
one of the Bianchi-related equations, and solve the remaining edge equations
for the next geometry.  In one published evolution they choose
`tau_i=rho*l_i`, i.e. an approximately constant conformal-time step.  Thus

```text
conserved dust mass  =>  a uniquely selected lapse
```

is **not** a known consequence and is already false as a general statement
about the standard Regge evolution framing.  A selected lapse in the present
fully varied finite carrier would have to be classified first as a possible
pseudo-constraint, not immediately as an emergent physical clock.

The distinction between global and local variation is also known to be
load-bearing.  Liu and Williams compare the two in closed Regge--FLRW models
and report that local variation does not generally yield the same viable
model:

- Liu and Williams, *Regge calculus models of the closed vacuum
  Lambda-FLRW universe*, <https://arxiv.org/abs/1501.07614>.

The canonical action/principal-function interpretation and the possibility
that curvature breaks exact lapse gauge into background-dependent
pseudo-constraints are known:

- Dittrich and Hoehn, *Canonical simplicial gravity*,
  <https://arxiv.org/abs/1108.1974>;
- Bahr and Dittrich, *(Broken) Gauge Symmetries and Constraints in Regge
  Calculus*, <https://arxiv.org/abs/0905.1670>.

## 3. CONTROL

At `L_+=L_-=L0`, the new implementation must recover the exact upstream
family for every admitted positive `tau`:

```text
diagonal = L0^2-tau^2,
all 35 internal equations = 0,
S_grav + S_dust = 0,
p_pre,e  = -epsilon3*L0*tau/4,
p_post,e = +epsilon3*L0*tau/4.
```

It must also reproduce the exact unit-600-cell adjacent dot product
`phi/2`, edge chord `1/phi`, and the homothetic diagonal formula before the
Regge action is evaluated.

## 4. OPEN difference

The following questions are open for the repository's complete order-24
staircase action:

1. Does fixed `M0` make the global homothetic lapse equation nontrivial when
   `L_+ != L_-`?
2. If so, does it determine only the invariant ratio `(L_+-L_-)/tau`, as a
   discrete Friedmann/Hamiltonian constraint normally would, or does it
   determine an absolute lapse?
3. Does any resulting global root satisfy all 35 unaveraged internal orbit
   equations, or is it only a global-variation minisuperspace solution?
4. Is there a locally connected non-static homothetic branch through the
   published sandwich?
5. If the internal equations hold, does the canonical pre-momentum match the
   already derived forward datum without an extra chosen impulse?

No located primary source prints this complete comparison for the exact
order-24 quotient and parity pair.  External novelty remains **OPEN**.

## 5. Framing verdict before calculation

The proposed dichotomy

```text
either mass conservation selects tau, or L is forced constant
```

is incomplete.  A third and standard possibility is that the equation fixes
only a velocity-like ratio such as `(L_+-L_-)/tau`, leaving the absolute
parameterization of the same geometric history free.  A fourth possibility
is a global-variation root that fails the 35 local equations because the
staircase triangulation breaks the homogeneous reduction.

Accordingly, a nonzero lapse equation is not an acceptance result.  The route
advances to a physical tick only if the same geometrically fixed slab passes
all local internal equations and the independent canonical junction.  A
lapse selected only by the finite triangulation is **STRUCTURAL / candidate
pseudo-constraint** until it survives refinement or an exact symmetry
argument.
