# Exact rational canonical-data lift

Date: 2026-08-19

## Complete hypotheses

On the complete variable-face flat-frustum 600-cell slab, use the exact
rational first-order gluing equations

```text
F f + E e + S s = 0,
```

with 3600 cell-flex variables, 720 upper spatial squared-edge data, and 120
strut squared-length data.  Restrict the edge data to the derived unsigned
vertex-scale image

```text
e_{uv} = 8 lambda (sigma_u + sigma_v)
```

and retain arbitrary strut data.  This is a kinematic gluing calculation, not
an action, Hessian, constraint algebra, or evolution equation.

## Frozen provenance

- prior-art/framing gate: commit `4de7e8f`;
- target-disclosed protocol: commit `1cbe2e1`;
- registered exact implementation before first execution: commit `c02ad2b`;
- first artifact: commit `16ae791`;
- artifact SHA-256:
  `1b6ac46a0ea4889f476cc71d51ac464c27caa6d4b6a9b2f2d74ff93da77b123f`.

The targeted verifier passed 12/12 checks.  No full-suite run was performed.

## Exact result

For both rational representatives, both local right-inverse graphs, reversed
face orientation, odd canonical relabelling, and reversed metric sign:

```text
exact pivots of F                 3600 / 3600
candidate consistency             true
nonzero rows after direct check       0 / 6000
corrupted incidence rejected       true
nonzero coefficients in the lift  28800
```

The solution was derived by sparse rational elimination on 3600 pivot
equations and then substituted independently into every original equation,
including all 2400 non-pivot rows.  Every residual coefficient vanished
exactly.  No tolerance, floating-point rank, or finite-field inference was
used.

The one-row-corrupted incidence map has the same rank 120 as the candidate
but leaves an exact obstruction.  For the two baseline representatives the
first obstruction is, respectively,

```text
(lambda,tau)=(2,5):  row 4200, data column 33, coefficient -13
(lambda,tau)=(3,11): row 4200, data column 33, coefficient -125/8.
```

Thus the calculation distinguishes the proposed image from a nearby image of
the same dimension.

**DERIVED.** Over `Q`, every one of the 120 vertex-scale and 120 strut basis
data has a unique cell-flex response in the declared coordinates, and the
resulting first-order geometry satisfies all complete face-gluing equations.

**DERIVED NEGATIVE, retained.** The earlier hand-chosen local lift is not this
solution and remains refuted.  The exact global solve repairs the lift, not
the old formula.

## Support census

For every construction, the exact lift has the following support counts:

```text
each of 3600 flex coordinates depends on exactly 8 data coordinates;
each of 120 vertex-scale data touches exactly 20 tetrahedral cells;
each of 120 strut data touches exactly 20 tetrahedral cells.
```

Every 600-cell vertex is incident on 20 tetrahedra, and a tetrahedron has four
vertices, supplying eight natural local data `(sigma, strut)`.

**PATTERN.** The counts are exactly those of a cell-local star-supported
formula.  The current artifact records counts, not the identity of every
support set, so exact equality with the vertex stars is not yet claimed.

## Physical meaning and strict limits

The result replaces an unidentified modular kernel by an exact rational
first-order data domain and response map:

```text
120 vertex conformal/scale data + 120 strut data -> unique glued cell flexes.
```

**STRUCTURAL.** This is a clean canonical-data carrier for the chosen
flat-frustum linearization.  It is compatible with known discrete-conformal
and flat-background vertex-displacement mechanisms; external novelty remains
**OPEN**.

It does not show that the 240 coordinates are physical degrees of freedom.
It does not select the arbitrary struts, a tick, or a time direction.  It
contains a conformal scalar sector, not tensor gravitons.  It derives neither
`c`, `G`, an absolute length, nor Planck units.

## Next falsifiable steps

1. Before comparing supports, preregister the exact claim that each datum is
   supported on precisely the 20-cell star of its vertex and each cell-flex
   row uses precisely the four scales and four struts at that cell's vertices.
2. If this survives, extract and compare the exact 6-by-8 local response
   blocks under cell relabellings, and reconcile them algebraically with the
   failed old local formula.
3. Only then restrict the Regge boundary action/Hessian to these 240 exact
   directions.  That calculation, not kinematic gluing, must determine the
   constraint/gauge split and whether any lapse/tick is selected.

