# No natural `A5`-equivariant Dirac can connect the 936-state carrier

Date: 2026-08-11

Protocol commit: `c2c32df`.

Registered verifier:
`reproducible/verify_hopf_six_equivariant_dirac_universal.py`.
Targeted exact result: `17/17`.

No Hessian, particle representation, mass, coupling or Standard-Model target
was used.  The full suite was not run, by explicit user instruction.

## Headline

The earlier normalized-rook calculation tested 32 particular operators.  The
new calculation removes that normalization restriction completely.

On the fixed 936-state off-diagonal carrier, for any of the eight derived
gradings, **no odd first-order Dirac operator natural under the full derived
`A5` action can satisfy connectedness**.

The result is independent of:

- equal versus unequal link magnitudes;
- zero versus nonzero legal blocks;
- signs and phases;
- any linear combination inside an equivariant Hom space.

The exact maximal-constraint census is

```text
spectral readings                                      8,
distinct legal central-link sets                       4,
readings with scalar maximal commutant               0/8,
maximal commutant dimensions             109 (x2), 141 (x2), 174 (x4),
maximal commutator ranks                  251 (x2), 219 (x2), 186 (x4).
```

This is a **DERIVED UNIVERSAL EQUIVARIANT NO-GO**, with the hypotheses stated
below.  The missing datum is no longer accurately described as merely a new
"noncentral incidence tensor".  It must either break the residual `A5`
symmetry in an independently selected way, or come with a different
carrier/algebra.

## Complete hypotheses

The theorem applies when all of the following hold:

1. the algebra is
   `B_R=M6(R)+M6(R)+M12(R)+M12(R)`;
2. the Hilbert carrier is
   `H_off=direct_sum_(i!=j) V_i tensor V_j*`, of real dimension 936;
3. the grading is one of the eight signed lexicographic readings selected
   from the exact joint spectrum `(u_edge,v_ref)`;
4. `D` is odd and satisfies first order for the full matrix algebra;
5. `D` is a scalar natural operator under every automorphism in the derived
   orientation-preserving `A5` action, equivalently it commutes with the
   induced diagonal `A5` action.

Self-adjointness and `JD=DJ` may restrict the candidate space further, but
are not needed for the obstruction.  The verifier deliberately gives the
candidate more freedom than a real self-adjoint Dirac has.

## Short analytic obstruction

Exact Frobenius reciprocity gives

```text
V0 = 1+5,
V1 = 3+3',
V2 = 3+4+5,
V3 = 3'+4+5.
```

All five `A5` irreducibles have Frobenius--Schur indicator `+1`, so these are
real multiplicity-one decompositions.  The complete real equivariant Hom
Gram matrix is

```text
G = [[2,0,1,1],
     [0,2,1,1],
     [1,1,3,2],
     [1,1,2,3]].
```

First order permits eight odd cell positions and three central node links
for each grading.  Across all eight readings:

- the link `V2--V3`, with its two-dimensional Hom space, never occurs;
- a link from `V0` to `V2` or `V3` sees only their common `5`;
- `Hom_A5(V0,V1)=0`.

Therefore the projector

```text
p = projector onto the trivial 1 inside V0
```

is an element of the `M6(R)` summand of the algebra and commutes with every
legal equivariant first-order block.  It is neither zero nor the global
identity.  Hence

```text
{a in B_R : [D,pi(a)]=0} != R 1
```

for every `D` satisfying the hypotheses.  This explicit projector is the
same exact non-scalar witness in all eight readings.

The obstruction is stronger than a failed numerical rank: one irreducible
sector is invisible to every admissible equivariant arrow.

## Why the maximal-span calculation proves a universal result

For each legal link `i--j`, the verifier places a complete exact basis of
`Hom_A5(V_i,V_j)` and imposes, separately for every basis tensor,

```text
A_j T = T A_i,
A_i T^t = T^t A_j.
```

This simultaneous system is stricter than the commutation system of any one
Dirac operator.  If an algebra element commutes with every basis tensor, it
commutes with every linear combination, independently of coefficients.
Thus the maximal common commutant is contained in the commutant of every
individual equivariant `D`.

Even this deliberately overconstrained system has dimensions 109, 141 or
174.  A particular `D` can only have an equal or larger commutant.  The
previous 109/141 values are recovered exactly on the four edge-first
readings; the four reflection-first readings are still weaker because their
required `V0--V1` Hom space is zero.

## Canonicity audit

The inference from geometry to equivariance is valid only with a precise
naturality clause:

> A scalar operator constructed from the unlabelled incidence-decorated
> geometry, without a chosen vertex, chamber, fibre or external field, must
> be fixed by every automorphism of that geometry and hence by its derived
> `A5` subgroup.

This does not say that every physically meaningful Dirac operator must be
`A5`-invariant.  It says that the current symmetric geometry cannot select a
single non-invariant operator by itself.

An `A5`-covariant family of tensors is a legitimate new object, but choosing
one component of that family requires an additional vacuum or dynamical
selection law.  Calling one component "canonical" merely because the whole
family is canonical would hide exactly the symmetry-breaking choice that
must be explained.

## Relation to crossed-product spectral-triple constructions

Published crossed-product constructions start from a spectral triple on the
base and add group data such as a proper Dirac weight or length function.
They can preserve first order and orientability under hypotheses, but they do
not canonically supply the missing `A5`-breaking tensor here.  In the current
programme the six-label base already has an orientability/KO6 obstruction,
and on a finite group properness alone does not select a unique length.

Primary references inspected for this scope check:

- A. Rubin and L. Dąbrowski,
  [*Real Spectral Triples on Crossed Products*](https://arxiv.org/abs/2012.15698);
- P. Antonini, D. Guido, T. Isola and A. Rubin,
  [*A Note on Twisted Crossed Products and Spectral Triples*](https://arxiv.org/abs/2110.05345).

This literature therefore provides a framework for a future dynamical input,
not an automatic repair of the present no-go.

## Status ledger

- **DERIVED:** all five `A5` irreps used here are of real type.
- **DERIVED:** the exact induced modules and Hom Gram matrix above.
- **DERIVED:** every one of the eight gradings has six positive cells, eight
  legal odd positions and three central links.
- **DERIVED:** the complete maximal equivariant Hom spans have commutant
  dimensions 109, 141 or 174.
- **DERIVED:** the projector onto `1 subset V0` is an exact non-scalar
  commutant witness for all eight readings.
- **DERIVED UNIVERSAL EQUIVARIANT NO-GO:** connected hit fraction `0/8` for
  the complete hypothesis list above.
- **STRUCTURAL:** selecting a member of a non-invariant covariant tensor
  family without a vacuum law.
- **OPEN:** an independently derived `A5`-breaking order parameter and its
  dynamics, or a different carrier/algebra that removes the invisible
  trivial sector.
- **NO TARGET COMPARISON:** no desired matter or phenomenological character
  was examined.

## Programme boundary and next legitimate move

Searching more `A5`-invariant incidence formulas on this carrier is now
mathematically redundant: all of them lie inside the maximal Hom span already
killed above.

There are two honest continuations:

1. derive an `A5`-covariant field space from local 600-cell incidence and ask
   whether the theory's own action selects a non-invariant vacuum whose
   induced `D` is connected; or
2. replace/reduce the 936-state maximal carrier by a geometry-selected
   correspondence in which every algebra sector is seen by the legal arrows.

The first route is the direct place where inertia, masses and spontaneous
symmetry breaking could eventually enter.  It must be preregistered at the
level of the complete covariant field space before any favorable vacuum is
chosen.
