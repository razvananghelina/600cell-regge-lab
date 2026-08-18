# Prior-art gate: full vertex-lapse Schur sector of the 600-cell dust tick

Date: 2026-08-17  
This note precedes the new precision protocol and calculation.

## Exact proposed object and hypotheses

Keep both already-derived five-stage schedule parities, the accepted
non-static homothetic dust slab, its complete 2,400 Lorentzian 4-simplices,
the uniform 120-pole dust action, and all 2,280 logarithmic signed-squared
edge variables.  Let

```text
J : (840 internal + 720 new boundary variations)
      -> (840 internal equations + 720 old momenta)
```

be the complete `1560 x 1560` pre-Legendre Jacobian certified in commit
`c874ef9`.

The 120 timelike pole edges form exactly five free regular orbits of the
binary tetrahedral schedule stabilizer `2T`; no orbit mixes a pole with a
spacelike diagonal.  Partition the rows and columns of `J` into the 120 pole
coordinates/equations and the remaining 1,440 directions.  The proposed new
object is the pole Schur operator

```text
S_pole = D - C A^(-1) B,

J = [ A  B ],
    [ C  D ]
```

resolved separately in the seven irreducible `2T` sectors.  Its minimal block
sizes are `5d`, while the eliminated strong blocks have sizes `60d`, for
`d=1,1,1,2,2,2,3`.

Also define, before inspecting `S_pole`, 120 geometric vertex-lapse
directions.  For a new vertex `v+120`, vary its pole logarithm by one, every
spacelike cross-edge ending at that vertex by

```text
-rho / (exp(s) L0^2-rho),
```

and no new-boundary spatial edge at first order.  Their sum is exactly the
already frozen collective internal-lapse direction.  This definition follows
from differentiating `q=spatial_square-rho` under a temporal displacement of
one new vertex; it is not fitted to a weak singular vector.

## Primary prior art

### KNOWN

- Dittrich and Hoehn derive a canonical formalism from the discrete action in
  which pre/post constraints and the Legendre maps encode simplicial
  evolution.  See [From covariant to canonical formulations of discrete
  gravity](https://arxiv.org/abs/0912.1817) and [Canonical simplicial
  gravity](https://arxiv.org/abs/1108.1974).
- Bahr and Dittrich show that curvature generically breaks exact vertex
  displacement symmetry in Regge calculus, replacing continuum-like
  constraints with background-dependent pseudo-constraints:
  [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670).
- Hoehn identifies four lapse/shift variables and their vertex-displacement
  generators around flat backgrounds, and distinguishes them from propagating
  lattice curvature degrees of freedom: [Canonical linearized Regge
  Calculus](https://arxiv.org/abs/1411.5672).
- De Felice and Fabri study only the symmetry-reduced dust 600-cell evolution,
  including its causal stopping point: [The Friedmann universe of dust by
  Regge Calculus](https://arxiv.org/abs/gr-qc/0009093).

### CONTROL already in the repository

- The old static, orbit-reduced `35 x 35` internal Hessian exposed a five-pole
  Schur sector with one collective null and four relative stiffnesses.  That
  calculation explicitly did not cover all 840 internal edges.
- The accepted non-static tick solves all individual 840 internal equations
  and reproduces all 1,440 boundary momenta.
- The complete binary64 census gives a full `1560/1560` rank but a sharply
  separated 120-dimensional weak cluster near `4.2445e-9`.  Its minimal
  multiplicities are `5d`, exactly the representation count of five regular
  `2T` modules.  The weakest block clears the preregistered global threshold
  by only `1.092`, motivating a precision correction rather than immediate
  inversion.

### OPEN

- Whether the full 120-dimensional pole Schur operator is genuinely
  invertible at the non-static tick.
- Whether its weak subspace is the canonical lift of the 120 geometric
  vertex-lapse directions just defined.
- Whether the small stiffnesses are broken gauge directions,
  pseudo-constraints, dust-clock dynamics, or a mixture.
- Whether any complementary directions have the two-polarization and
  dispersion properties of gravitons.

## Search result and proposed difference

The pre-calculation search found the general canonical and
pseudo-constraint mechanism, but no primary source computing this exact
complete dust-600-cell pole Schur operator, its binary-tetrahedral
decomposition, or the proposed 120 vertex-lapse comparison.  This search does
not prove novelty; external novelty remains **OPEN**.

The calculation is worthwhile even if negative.  A certified zero closes the
current ungauged inversion route and demands an explicit quotient.  A
certified nonzero Schur operator supports local invertibility but still does
not turn the lifted lapse variables into physical gravitons.
