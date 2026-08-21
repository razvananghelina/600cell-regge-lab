# Prior-art gate: global scale covariance and absolute tick selection

Date: 2026-08-21

Status: completed before the scale-covariance protocol and before evaluating
the complete action at any rescaled state.

## Exact question and complete hypotheses

Consider the repository's classical Lorentzian flat-simplex Regge action on a
fixed one-slab 600-cell staircase carrier, including its boundary terms and
pressureless point-dust worldline action.  Write every signed squared edge
length as `q_e`, every positive dust proper-length square as `rho_v`, and every
geometrized point mass as `m_v`.  There is no cosmological term, quantum term,
fixed external length, or dimensionful coupling other than the units in which
`G=c=1`.

The proposed transformation is, for one positive real `alpha`,

```text
q_e   -> alpha^2 q_e       for every spatial and temporal squared length,
rho_v -> alpha^2 rho_v,
m_v   -> alpha m_v.
```

The narrow question is whether the complete action and its equations are
exactly covariant under this transformation.  If they are, can this classical
model select an absolute nonzero proper duration, or only dimensionless ratios
such as `tau/L` and relative lapse ratios?

The statement is conditional on scaling the matter masses.  Holding a
geometrized mass fixed supplies an external length scale and is a different
problem.  Likewise, a fixed cosmological constant or any other dimensionful
coefficient is outside the hypotheses.

## KNOWN from primary literature

- Ordinary Regge calculus replaces the Einstein--Hilbert curvature integral
  by hinge volume times deficit angle.  In four dimensions the hinge volume is
  an area; under a uniform rescaling of lengths it has degree two, while the
  dimensionless deficit angle is invariant.  This is standard Regge structure,
  not a new 600-cell mechanism.  See the review
  [Quantum Gravity and Regge Calculus](https://arxiv.org/abs/gr-qc/9701052).
- Consistent discretizations can determine continuum multiplier-like variables
  and implement discrete dynamics as a canonical transformation, but this does
  not by itself manufacture an absolute unit:
  [Consistent discretization and canonical classical and quantum Regge calculus](https://arxiv.org/abs/gr-qc/0511096).
- Dust can supply a relational reference frame and a proper-time variable.  It
  does not follow that the coupled classical equations select the normalization
  of a dimensionful unit from no scale:
  [Dust as a Standard of Space and Time in Canonical Quantum Gravity](https://arxiv.org/abs/gr-qc/9409001).
- A recent Lorentzian Regge cosmology similarly uses matter to deparametrize a
  relational Friedmann equation and treats the frustum height dynamically:
  [Cosmology in Lorentzian Regge calculus: causality violations, massless scalar field and discrete dynamics](https://arxiv.org/abs/2312.11639).

These sources establish the broad mechanisms.  They do not serve as evidence
that the repository implementation obeys the claimed scaling identity.

## CONTROL from the repository

- The complete 2400-simplex action reproduces the published dust sandwich.
- On the regular one-slab family, the collective duration is an exact null
  direction after the associated cross lengths change with it.
- With fixed incoming canonical data, the homothetic two-variable solve selects
  a locally unique *relative* next lapse,
  `tau_next/tau0 = 0.999998220375...`, while `tau0=0.0102` remains supplied.
- The current selected dust rules are geometrized curvature masses.  Spatial
  curvature has length degree one, so those rules are expected to transform as
  `m_v -> alpha m_v`; this expectation must be checked rather than assumed.

## Dimensional derivation to be audited

Let `r=alpha^2`.  A simplex Gram matrix is linear in squared lengths and hence
scales by `r`.  Its dihedral-angle cosine and sine are ratios with equal total
degree, hence angles and deficit angles are invariant.  A triangle area scales
by `r`.  Therefore

```text
S_Regge(alpha^2 q) = alpha^2 S_Regge(q).
```

The dust action is a sum of `-m_v sqrt(rho_v)` (up to the fixed normalization),
so simultaneous mass scaling gives

```text
S_dust(alpha^2 rho, alpha m) = alpha^2 S_dust(rho,m).
```

Consequently the complete action is homogeneous of length degree two.  Its
derivatives with respect to logarithmic squared lengths scale by `alpha^2`;
its derivatives with respect to raw squared lengths are invariant.  Zeros of
the internal equations therefore occur in continuous scale families, and
canonical momenta scale covariantly.

This dimensional argument is **STRUCTURAL** until it is checked against both
independent complete-action implementations and the actual selected mass rule.

## Framing attack

Even a successful audit will not prove that time is gauge in every discrete
model.  It proves only that the stated scale-free action cannot select an
absolute unit.  Finite Regge pseudo-constraints may still select a lapse *ratio*
or a height relative to boundary geometry.  Conversely, a fixed input mass may
select a length relative to that mass, but then the scale was supplied by the
mass and was not derived from the 600-cell.

The no-go is also classical.  A quantum anomaly, a cosmological constant,
boundary data with physical units, or another dimensionful matter parameter
could break it.  None may be introduced after seeing a preferred numerical
tick and then called a derivation.

## OPEN before calculation

- exact covariance of the complete orbit action and all 95 logarithmic
  derivatives at an off-shell nonhomogeneous state;
- independent covariance of the direct 2400-simplex action and raw
  squared-length derivatives;
- exact length degree of the selected coarse and refined curvature-mass rules;
- whether the current action can select any absolute nonzero `tau` without an
  external dimensional anchor;
- which independently motivated scale-breaking ingredient, if any, belongs in
  the theory.

External novelty of the repository-specific audit is **OPEN**.  Search absence
cannot establish novelty.

