# Whitney--Kähler--Dirac inductive dynamics

Date: 2026-08-10

## Decision

There is now a genuine, nontrivial inductive spectral dynamics on the full
oriented cochain carrier of the barycentrically refined 600-cell boundary.
The construction is not the earlier ultrametric path model and contains no
free level-weight sequence.

At every level `n`:

- `K_n` is the full piecewise-flat barycentric subdivision;
- `V_n^k` is the lowest-order Whitney `k`-form space for `k=0,1,2,3`;
- the Hilbert metric is the exact geometric `L2` integral;
- `d_n` is the simplicial exterior derivative;
- `D_n=d_n+d_n*` is the metric Kähler--Dirac operator;
- `H_n=c D_n+gamma_n M` gives a nontrivial metric-unitary evolution.

The Whitney inclusions are exactly isometric, commute with `d`, and compress
the weak Kähler--Dirac form exactly.  This is an **ACCEPTED DERIVED
GALERKIN-INDUCTIVE SPECTRAL DYNAMICS**.  It satisfies the construction goal in
the mathematical sense.

It is not yet an accepted fundamental spacetime dynamics.  Two independent
obstructions remain:

1. on the calibrated circle, the consistent Whitney mass permits cutoff-scale
   group velocity `sqrt(2)c`; row-sum lumping restores the unit speed bound but
   breaks exact inductive isometry;
2. iterated tetrahedral barycentric subdivision is not shape-regular, so the
   standard FEEC convergence theorems cannot simply be applied to this exact
   tower.

Thus the construction advances the programme from **OPEN: no dynamics** to
**DERIVED: exact Galerkin dynamics / OPEN: causal continuum or fundamental
interpretation**.  It does not derive particle masses, physical units, a
fourth dimension, or Lorentzian causality.

Targeted certificates:

- `reproducible/verify_whitney_kahler_induction.py` -- 21/21;
- `reproducible/verify_whitney_circle_calibration.py` -- 10/10;
- `reproducible/verify_barycentric_shape_regular_gate.py` -- 12/12.

## 1. Complete hypotheses

Every positive claim in this note assumes all of the following.

1. The spatial object is the piecewise-Euclidean boundary complex of the
   regular 600-cell, topologically `S^3`.
2. Refinement is full affine barycentric subdivision inside every tetrahedral
   facet, with compatible face gluings.
3. All cochain degrees `0,1,2,3` are retained; this is not a top-cell-only
   model.
4. The finite spaces are the lowest-order Whitney spaces, with one degree of
   freedom per oriented `k`-simplex.
5. Their scalar products are exact `L2` integrals in the inherited
   piecewise-flat metric, not counting metrics or diagonal approximations.
6. Induction means isometric nesting and exact Galerkin compression of
   quadratic/weak operator forms.  It does **not** mean the stronger identity
   `D_(n+1) I_n=I_n D_n`.
7. The time parameter in `exp(-itH_n)` is an external real parameter.  No
   conversion to seconds is claimed.
8. In `H_n=cD_n+gamma_n M`, `c` and the internal self-adjoint `M` are supplied
   parameters.  The mass-shell identity is derived conditional on them; their
   values are not selected here.
9. A continuum interpretation requires a separate convergence theorem for
   this particular, geometrically degenerating mesh family.

Changing Hypotheses 2, 4 or 5 changes the dynamics.  In particular, mass
lumping is not silently substituted for the exact Whitney metric.

## 2. Exact all-degree construction

Let `lambda_i` be barycentric coordinates on a tetrahedron.  For an oriented
`k`-simplex `I=(i_0,...,i_k)`, its Whitney form is

`omega_I = k! sum_r (-1)^r lambda_(i_r)
                 d lambda_(i_0) wedge ... omit r ... wedge d lambda_(i_k)`.

On the full barycentric subdivision, a fine simplex is a chain of nonempty
faces.  The inclusion coefficient of a coarse Whitney form on a fine simplex
is its canonical degree of freedom, equivalently the determinant of the
corresponding barycentric-coordinate minor.  Denote the resulting inclusion
by `P_k`.

The exact rational reference-tetrahedron calculation gives, simultaneously
for `k=0,1,2,3`,

`d_f P_k = P_(k+1) d_c`,

`P_k^T M_(f,k) P_k = M_(c,k)`.

The first identity is cochain compatibility.  The second is exact `L2`
isometry.  Both are local affine identities; summing them over tetrahedra
gives the global identities on every barycentric level of the 600-cell.

The local f-vectors used in the certificate are

`(4,6,4,1) -> (15,50,60,24)`.

All eight local mass matrices are symmetric.  Their smallest numerically
audited eigenvalue is positive (`0.02449686...`); positivity itself follows
structurally because they are Gram matrices of linearly independent Whitney
bases.

## 3. The metric Kähler--Dirac form

For the block metric `M_n=diag(M_(n,0),...,M_(n,3))`, define the
codifferential by

`delta_(n,k)=M_(n,k)^(-1) d_(n,k)^T M_(n,k+1)`

and

`D_n=d_n+delta_n`.

It is preferable to store the weak operator

`A_n=M_n D_n`.

`A_n` is exactly symmetric and odd under form parity `gamma_n`.  With
`P=diag(P_0,...,P_3)`, the verifier proves

`P^T M_f P=M_c`,

`P^T A_f P=A_c`,

`gamma_f P=P gamma_c`.

Therefore `D_n` is an exact Galerkin compression of `D_(n+1)`.  This is the
operator content missing from the former arbitrary level-weight construction.

The fine weak stencil has 820 directed kinetic entries on the reference
subdivision.  Every one connects simplices contained in a common top
tetrahedron; there are zero locality violations.  **DERIVED:** the weak
variational generator is simplex-star local.

The strong coefficient matrix `M_n^(-1) A_n` need not have the same sparse
support because the inverse consistent mass is generally dense.  A strict
finite-level light cone is therefore not inferred from weak locality.

## 4. Compression is not strong intertwining

Exact compression does not make the inherited coarse space reducing.  The
weak residual equivalent to codifferential intertwining is

`R_k=d_f^T M_(f,k+1)P_(k+1)-M_(f,k)P_k delta_(c,k)`.

Its exact ranks are

`(rank R_0, rank R_1, rank R_2)=(3,3,0)`.

Thus:

- **DERIVED NEGATIVE:** lower degrees leak into new vertical modes;
- **DERIVED:** the top `2 <-> 3` codifferential does intertwine exactly;
- **DERIVED:** every leakage is orthogonal to the inherited sector and is
  invisible under coarse compression.

This corrects the earlier expectation that all three degrees would fail.
Strong intertwining is not part of the accepted construction.

## 5. Nontrivial dynamics and the mass shell

The finite generalized eigenproblem is

`A_n v=lambda M_n v`.

On the reference tetrahedron and its subdivision:

- the Kähler--Dirac kernels both have dimension one, the Betti sum of a
  contractible tetrahedron;
- the spectral radii are `3.872983...` and `11.936870...`;
- nonzero eigenvalues have exact grading-forced `+/-` pairing.

The dynamics is therefore not a static inclusion or a phase attached to the
level.

Because `{D_n,gamma_n}=0`, the Hamiltonian

`H_n=cD_n+gamma_n M`

satisfies

`H_n^2=c^2D_n^2+M^2`

when `M` commutes with the spatial factor.  A spatial eigenvalue `p` and an
internal eigenvalue `mu` give

`E^2=c^2p^2+mu^2`.

The same relation supplies rest mass, inertial curvature and continuum
limiting speed.  The finite-level numerical witness has mass-shell residual
`1.42e-14`; `exp(-itH_f)` preserves the exact Whitney metric to `8.89e-15`
and is nontrivial.  **DERIVED CONDITIONAL:** the algebraic unification works
for supplied `c,M`.  **OPEN:** selecting `c,M` remains a different contract.

## 6. Known-answer circle calibration

Before interpreting the tetrahedral tower, the identical Whitney construction
was calibrated on the unit circle.  The exact continuum first positive Dirac
eigenvalue is `2 pi`.

For `N=(8,16,32,64,128)` vertices, the errors are

`(1.625e-1, 4.045e-2, 1.010e-2, 2.524e-3, 6.308e-4)`.

Successive asymptotic ratios are

`(4.005,4.001,4.000)`,

certifying second-order convergence.  Fixed physical momenta and the massive
group velocity converge monotonically as well.  This agrees with the FEEC
Hodge--Dirac stability/convergence theory of
[Leopardi and Stern](https://arxiv.org/abs/1401.1576) and the underlying FEEC
framework of
[Arnold, Falk and Winther](https://arxiv.org/abs/0906.4325).  Energy-preserving
mixed FEEC time evolution for the Hodge wave equation is also established by
[Wu and Bai](https://arxiv.org/abs/2009.02844).

These references validate the mathematical class of construction.  They do
not prove convergence on the particular iterated barycentric 600-cell tower
without their mesh hypotheses.

## 7. Exact causality/induction tradeoff on the calibration

For the consistent Whitney mass on a uniform circle mesh of spacing `h`, the
positive momentum dispersion is

`p_h(k)=2 sin(kh/2)/(h sqrt(2/3+cos(kh)/3))`.

At fixed physical `k`, this converges to `k`.  But at cutoff-scale phase
`q=kh=2pi/3`,

`dp_h/dk=sqrt(2)`.

Thus a Hamiltonian with coefficient `c` permits group velocity `sqrt(2)c` at
that lattice scale.  **DERIVED NEGATIVE:** the exact consistent Whitney
metric does not enforce the strict finite-level bound `|v|<=c`.

Canonical row-sum mass lumping changes the dispersion to

`p_h^lump(k)=2 sin(kh/2)/h`,

whose derivative is `cos(kh/2)` and is bounded by one.  However, its
refinement residual has exact rank three already for `4 -> 8` circle
vertices.  **DERIVED TRADEOFF:** among these two canonical metrics,
consistent mass preserves exact induction and lumped mass preserves the
finite lattice speed bound.  This is not a theorem that no third metric can
do both.

Dispersion analysis for mass-lumped wave elements on tetrahedral meshes is a
standard issue, not a new physical mechanism; see
[Geevers, Mulder and van der Vegt](https://arxiv.org/abs/1802.10333).

## 8. Shape-regularity kill gate

The standard FEEC convergence results assume controlled mesh geometry.  The
unmodified barycentric tower fails that gate by an explicit deterministic
path.

Repeat the flag

`(v0)<(v0,v1)<(v0,v1,v2)<(v0,v1,v2,v3)`.

In parent edge coordinates its affine transform is

```text
T = [1/2  1/3  1/4]
    [  0  1/3  1/4]
    [  0    0  1/4].
```

Hence

- `det(T)=1/24`;
- the eigenvalues are `(1/2,1/3,1/4)`;
- `det(T^n)=24^(-n)`;
- `cond_2(T^n)>=2^n`.

Direct singular-value calculation gives `cond(T^10)=14980.55...`; in the
physical regular-tetrahedron coordinates it is `18642.81...`.  Therefore the
family is not shape-regular.  The fact that iterated barycentric simplices can
flatten is also documented probabilistically in dimension two by
[Diaconis and Miclo](https://arxiv.org/abs/1007.3385), but the three-dimensional
negative here is supplied by the exact matrix witness and does not depend on
that paper.

**DERIVED NEGATIVE:** citing standard shape-regular FEEC convergence as proof
for this tower would be invalid.

## 9. Acceptance and next boundary

### Accepted

- **DERIVED:** a construction at every barycentric level, with no fitted
  spectral coefficient sequence;
- **DERIVED:** exact all-degree cochain nesting and `L2` isometry;
- **DERIVED:** exact Galerkin compression of a self-adjoint, graded
  Kähler--Dirac form;
- **DERIVED:** simplex-star-local weak stencil;
- **DERIVED:** nontrivial spectra and metric-unitary evolution;
- **DERIVED:** second-order known-answer convergence on the shape-regular
  circle control;
- **DERIVED CONDITIONAL:** one mass shell from the graded product Hamiltonian.

### Not accepted

- strong operator intertwining;
- strict finite-level propagation bound for the consistent mass;
- convergence of the unmodified barycentric `600-cell` tower;
- a selected internal mass operator or value of `c`;
- Lorentzian reconstruction, four dimensions or matter physics.

### Next falsifiable decision

There are now exactly two honest routes.

1. **Emergent-continuum route:** prove Hodge--Dirac spectral convergence on
   this specific degenerating barycentric family despite loss of standard
   shape regularity, and show that the cutoff `sqrt(2)c` sector decouples.
2. **Fundamental-lattice route:** select, from geometry rather than fitting, a
   metric/refinement repair which is simultaneously inductive and obeys a
   finite propagation bound.

If neither can be done, this Whitney family remains a sound numerical
dynamics but the refinement-spacetime interpretation is closed.

## Status ledger

- **DERIVED:** full all-degree Whitney inductive system.
- **DERIVED:** exact weak Kähler--Dirac compression.
- **DERIVED:** nontrivial unitary spectral evolution.
- **DERIVED:** local weak finite-element stencil.
- **DERIVED NEGATIVE:** lower-degree strong leakage ranks `(3,3)`.
- **DERIVED:** top-degree strong leakage rank `0`.
- **DERIVED:** second-order circle spectral calibration.
- **DERIVED NEGATIVE:** consistent cutoff velocity reaches `sqrt(2)c`.
- **DERIVED NEGATIVE:** lumped mass breaks exact induction.
- **DERIVED NEGATIVE:** barycentric tower is not shape-regular.
- **OPEN:** convergence on that degenerating tower or a canonical causal
  repair.
- **OPEN:** selection of `M`, `c` and Lorentzian time.
- **NOT CLAIMED:** a complete spacetime or particle theory.
