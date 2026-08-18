# Hopf selector versus the three-dimensional continuum

Date: 2026-08-10

## Question and hypotheses

This note asks a narrower question than whether the exact sixth-order Hopf
selector exists.  That existence and its extrema were already certified.
The question here is whether its anisotropic part survives as a relevant
perturbation of a three-dimensional continuum limit.

The conditional RG statement below has all of the following hypotheses:

1. the continuum limit is governed by the three-dimensional `O(3)`
   Wilson--Fisher fixed point;
2. its order parameter is a three-component real field `n`;
3. the six-axis polynomial enters the continuum action as a local potential
   perturbation, without a different lower-dimensional operator carrying the
   same microscopic information;
4. linearized RG about that fixed point controls the regime being discussed.

The repository has not derived hypotheses 1--3 from its certified spectral
action.  In particular, the current action stops at `Tr(D^4)` and the
all-gate action audit found no generated sixth-order selector or its sign.

## Exact bridge to the spin-six perturbation

Let the six exact fivefold-axis projectors be `P_i` and set

```text
S6(n) = sum_i (n^T P_i n)^3.
```

The targeted exact verifier now also certifies

```text
H6(n) = S6(n) - (6/7) (n.n)^3,
H6 != 0,
Laplacian(H6) = 0.
```

`H6` is therefore a nonzero homogeneous harmonic polynomial of degree six.
Such polynomials form the rank-six symmetric-traceless irreducible
representation of `O(3)`.  The sum over the complete six-axis orbit makes
this particular vector invariant under the rotational icosahedral group.
The radial subtraction changes none of the extrema on the unit sphere.

**DERIVED:** our selector is not merely analogous to an icosahedral
anisotropy.  Its nonradial part is exactly an `O(3)` spin-six perturbation.

This agrees with the invariant-theory classification in
[Platonic Field Theories](https://arxiv.org/abs/1902.05328): for the full
icosahedral symmetry the independent even invariants start at degrees
`2, 6, 10`.  The degree-two generator is radial, so degree six is the first
place at which the six directions can be distinguished.

## Direct numerical RG evidence

The closest external falsifier found is not a fitted polynomial model.  It is
the generalized icosahedral spin model of
[Hasenbusch, Phys. Rev. B 102, 024406 (2020)](https://arxiv.org/abs/2005.04448).
Its microscopic spin takes exactly the twelve signed icosahedron vertices
(plus an optional zero state), so its six unoriented preferred axes are the
same orbit selected by `-g S6` for `g>0`.  Simulations reached cubic lattices
of linear size 400.  The paper finds emergent `O(3)` criticality and measures

```text
y_ico = -2.19(2).
```

With the convention `y = d - Delta` in `d=3`, this is

```text
Delta_ico = 5.19(2).
```

It is far above the relevance threshold `Delta=3`.  The exponent was
extracted from the decay of an explicit icosahedral anisotropy observable,
including checks over lattice-size cutoffs and leading-correction-suppressed
models.  It was not obtained by assuming that the spin-six operator was
irrelevant.

**STRUCTURAL (external numerical result):** in this explicit 3D model, even
the hard restriction to the twelve preferred directions flows toward the
`O(3)` fixed point.  The icosahedral perturbation is strongly irrelevant at
criticality.

**DERIVED CONDITIONAL:** under hypotheses 1--4, the repository's `H6`
coupling has the same linearized symmetry channel and therefore flows to zero
at the critical point.

## Independent checks and a circularity caught

[Six-loop perturbation theory](https://arxiv.org/abs/2208.04612) gives
resummed spin-six dimensions `5.45(1)` for `O(2)` and `4.99(2)` for `O(4)`.
The published table omits `O(3)`, so interpolation is only a **PATTERN**, not
evidence by itself.  It is nevertheless consistent with `5.19(2)`.

[Rong and Su](https://arxiv.org/abs/2311.00933) report an `O(3)` extremal
functional estimate `Delta_6 = 5.2252`, also consistent with the Monte Carlo
number.  However, their bootstrap setup imposed the gap `Delta_6 > 5.0`.
Therefore that computation could not have falsified irrelevance and must not
be cited as independent evidence for `Delta_6>3`.  It is a consistency check
only.  This caveat is load-bearing.

The perturbative fixed-point search in *Platonic Field Theories* finds new
icosahedral fixed points only below three dimensions.  That supports, but
does not prove, the absence of a distinct physical 3D icosahedral fixed point.

## What survives and what dies

**DERIVED NEGATIVE, conditional on the O(3) continuum hypotheses:** the
six-axis term cannot supply a distinct icosahedral critical continuum in
three dimensions.  Critical long-distance physics restores `O(3)` rather
than magnifying the selector.

This does not erase the six vacua in an ordered phase.  Hasenbusch explicitly
identifies the icosahedral perturbation as likely *dangerously irrelevant*:
it disappears at the critical fixed point while the thermodynamic ordered
phase can still retain only the twelve signed preferred directions.  This is
exactly the qualitative role our selector would need.

**STRUCTURAL, not derived here:** the six-axis selector is therefore a viable
vacuum-orientation term after ordering, not a source of the continuum or of
the phase transition itself.

**OPEN:** the theory still lacks all of the following:

- a derived three-component order parameter to which `H6` couples;
- a generated nonzero coefficient `g`;
- the sign `g>0` selecting the six Hopf axes rather than the ten threefold
  axes selected by the opposite sign;
- a derivation of the ordered phase and its relation to time, inertia, mass,
  or the internal algebra.

## Verdict and next boundary

This is promising in one precise sense: an independently studied 3D model
with exactly the same twelve-direction orbit realizes emergent `O(3)` at
criticality while retaining the possibility of discrete vacuum selection
below it.  The agreement `Delta_ico = 5.19(2)` versus the non-independent
bootstrap estimate `5.2252` is quantitatively strong.

It is not yet a dynamical result of this theory.  The next acceptance boundary
is no longer another invariant search.  It is to derive a local effective
action or certified spectral fluctuation containing `g H6`, with `g` nonzero
and its sign fixed.  Failure of every admissible action extension to generate
that term closes this selector route; adding it by hand remains only a
consistent model.
