# Prior-art gate: orbit/flag incidence reduction of the 600-cell slab action

Date: 2026-08-21

Status: completed after the frozen direct/orbit disagreement and before any
flag-orbit incidence matrix was enumerated.

## Exact object and hypotheses

Let the already-derived order-24 schedule stabilizer act on the finite sets of
four-simplices, triangles, and incident flags `(triangle subset simplex)` of
one 600-cell staircase slab.  The squared-length state is constant on the
corresponding edge orbits.

The direct action sums one dihedral contribution per incident flag and one
hinge term per triangle.  The current reduced action instead chooses one
representative per simplex orbit, evaluates its ten local hinges, accumulates
those angles by triangle orbit, and finally multiplies each hinge term by the
triangle-orbit size.

The question is whether that representative shortcut uses the exact flag
incidence multiplicities.  This is a finite combinatorial question independent
of `M16`, cosmology targets, roots, Hessians, or numerical fitting.

## KNOWN structure

For an automorphism group, orbit partitions are equitable: incidence row sums
between two orbits are constant.  A lossless quotient therefore uses the
actual quotient/incidence coefficients, not an unproved unit coefficient for
each pair of chosen representatives.  See the original-results discussion of
orbit partitions and lossless symmetry compression in
[Exploiting symmetry in network analysis](https://www.nature.com/articles/s42005-020-0345-z).

For the bipartite triangle--simplex incidence graph, an exact coefficient can
be obtained directly.  If `F_a` is an orbit of incident flags and `T_i` its
triangle orbit, then the number of flags from `F_a` incident to one fixed
triangle in `T_i` is

```text
c_a = |F_a|/|T_i|.
```

This follows by double counting and orbit transitivity.  It is exact rational
arithmetic.  No literature result predicts the specific coefficients of this
order-24 carrier.

## CONTROL from the repository

- The shortcut and full direct binary64 actions agree on the regular controls
  used by the older boundary-Legendre verifier.
- On the frozen order-24-invariant off-shell perturbation, the 100-decimal
  shortcut differs from both direct implementations by `1.3e-6`--`1.6e-6`.
- The binary64 and independently reimplemented 80-decimal direct actions agree
  with one another to about `2e-9` and both obey exact scale covariance.
- The high-precision adjudication was preregistered to select the shortcut but
  selected the direct action instead; that failed outcome is retained at
  commit `af862ab`.

## Framing attack

The disagreement does not yet prove a flag-multiplicity bug.  Other live
possibilities are:

1. old/internal/final orbit-coordinate orderings differ between the two model
   augmentations;
2. the frozen index perturbation is not invariant under the group action
   actually used by one evaluator;
3. boundary/internal hinge constants differ;
4. the independent direct action contains a repeated or missing triangle;
5. the representative shortcut omits or misweights flag orbits.

Every possibility must be checked before changing code.  Regular-state
agreement cannot falsify an incidence error whose contributions cancel when
all same-type angles coincide.

## OPEN before enumeration

- equality of all carrier and coordinate lookup maps;
- edge-state invariance under every stabilizer element;
- the exact triangle, simplex and flag orbit-size distributions;
- the exact `|F_a|/|T_i|` coefficient multiset;
- equality or inequality between those coefficients and shortcut counts;
- whether a flag-corrected reduced action matches the direct action at both a
  regular control and the frozen off-shell discriminator;
- the scope of any correction across prior nonhomogeneous results.

External novelty is **OPEN**.  The general quotient principle is known; only
the repository-specific audit could be new.

