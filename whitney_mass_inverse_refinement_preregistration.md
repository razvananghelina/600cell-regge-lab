# Preregistration: does Whitney inverse complexity grow at first refinement?

Date: 2026-08-11

## Fixed question and framing boundary

Protocol/result commits `3323174` / `a87cc9d` proved exact inverse-polynomial
degrees

\[
(8,21,26,0)
\]

for the four assembled Whitney mass blocks on the unrefined 600-cell.  These
are finite and far below Cayley--Hamilton, but no fixed-complex result can
show bounded depth under refinement.

Construct the complete first barycentric subdivision of the 600-cell and ask
the falsifiable one-step question:

> Does an exact lower bound for each refined minimal-polynomial degree already
> exceed the corresponding coarse exact degree `(9,22,27,1)`?

One refinement cannot prove divergence or a continuum law.  An increase
refutes only the strongest constant-depth hypothesis: that the base degrees
already suffice unchanged at the first refined level.

## Frozen fine complex and masses

The vertices of the barycentric subdivision are all nonempty coarse cells.
A fine tetrahedron is a maximal flag

\[
v\subset e\subset f\subset t.
\]

Generate all 24 flags inside every one of the 600 parent tetrahedra, then
deduplicate their faces.  The required refined f-vector is

\[
(2640,17040,28800,14400).
\]

Use the exact rational Whitney masses of a barycentric flag child of the
regular reference tetrahedron.  Verify explicitly that all 24 child orderings
give the same local mass after using flag-rank vertex order.  Assemble the
four complete refined mass blocks, clear denominators and divide by the gcd
of all entries to obtain primitive symmetric integer matrices.

No diagonal lumping, fitted weights or numerical coordinates are permitted.

## Exact lower-bound method

Full refined polynomials and dense inverse support are deliberately not
required in this first gate.  For a matrix (B), fixed integer vector (v)
and prime (p), compute the exact modular sequence

\[
s_k=v^TB^kv\pmod p.
\]

Berlekamp--Massey applied to any finite prefix gives a linear complexity
(L).  If a degree-(r) polynomial annihilates (B) over the rationals, its
monic reduction modulo (p) annihilates this sequence, hence (r\geq L).
Therefore every reported (L) is a deterministic exact lower bound, even if
the probe misses other spectral sectors.

Freeze:

- probe seed `60020260811 + form_degree`;
- primes `1000003`, `1000033`, `1000037`;
- sequence length (4s_p+32), where
  (s_p=(9,22,27,1)) is the preregistered exact coarse degree;
- the same probe-generation rule and Berlekamp--Massey implementation at both
  levels.

## Calibration and gates

1. **Coarse calibration:** for every form degree and all three primes, the
   prefix complexity must reproduce the exact degree `(9,22,27,1)`.  Failure
   invalidates the estimator rather than the geometry.
2. **Fine exact lower bounds:** record all 12 complexities.  The certified
   degree lower bound is their maximum for each form degree.
3. **Top-form control:** every refined top form belongs to one fine
   tetrahedron, all children have equal volume, and the primitive top mass
   must remain the identity.  Its complexity must remain one.

## Decision labels

- **DERIVED DEGREE GROWTH AT LEVEL 1:** the fine lower bound exceeds the
  exact coarse degree for a form block.  Then the old inverse polynomial
  cannot work unchanged after one refinement.
- **DERIVED NEGATIVE FOR GROWTH DETECTION:** the lower bound does not exceed
  the coarse degree.  This does not prove equality because a finite probe can
  miss factors.
- **REFUTED PROTOCOL:** any coarse calibration fails, the fine f-vector is
  wrong, child masses are not congruent as assumed, or exact assembly fails.

The result must report the ratio of the certified fine lower bound to the
coarse exact degree, but must not fit a scaling exponent from only two levels.
It must not claim that infinity exists.  At most it can show that fixed-depth
locality has already failed at the first refinement.

Only the new targeted verifier will be run.  The full suite remains excluded
by user instruction.
