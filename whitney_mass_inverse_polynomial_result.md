# The Whitney inverse is finite and spectrally short, but globally supported

Date: 2026-08-11

Preregistration commit: `3323174`

Targeted verifier:
`reproducible/verify_whitney_mass_inverse_polynomial.py`

Targeted result: **10/10 PASS**.  The verifier is registered in
`reproducible/run_all.py`, but the full suite was not run by explicit user
request.

## Framing result

The statement “a finite inverse polynomial exists” is automatic for every
finite invertible matrix and has no evidential weight by itself.  The
nontrivial result is that the exact Whitney mass algebra closes at degrees
far below the carrier dimensions, while the inverse nevertheless reaches the
full graph diameter.

> **DERIVED:** infinity is not algebraically required to implement the exact
> Whitney inverse on the fixed 600-cell.  It is an exactly selected finite
> polynomial of the local upper Laplacian.

> **DERIVED NEGATIVE:** except on top forms, that inverse is not one-step
> local.  Its exact support reaches the diameter of every connected mass
> graph.

Both statements are needed.  “Finite” does not mean “microscopically local.”

## Exact mass identities

The local Whitney masses were independently regenerated from their defining
affine-form integrals and assembled on the full complex.  After clearing
denominators and dividing by the integer gcd, the primitive mass blocks obey

\[
\begin{aligned}
B_0&=20I-d_0^Td_0,\\
B_1&=50I-3d_1^Td_1,\\
B_2&=20I-d_2^Td_2,\\
B_3&=I.
\end{aligned}
\]

All four residuals are exact zero integer matrices.  Thus the metric inverse
does not introduce a new fitted operator: it belongs to the spectral algebra
of the geometry's own upper Hodge Laplacians.

The actual rational masses are

\[
M_0=\frac23B_0,\qquad
M_1=\frac1{60}B_1,\qquad
M_2=\frac1{60}B_2,\qquad
M_3=\frac38B_3.
\]

## Exact polynomial and support census

| degree | dimension | minimal-polynomial degree | inverse degree | graph diameter | farthest inverse support | inverse nonzeros |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 120 | 9 | 8 | 5 | 5 | 14,400 / 14,400 |
| 1 | 720 | 22 | 21 | 10 | 10 | 446,400 / 518,400 |
| 2 | 1,200 | 27 | 26 | 15 | 15 | 1,310,400 / 1,440,000 |
| 3 | 600 | 1 | 0 | 0 | 0 | 600 / 360,000 |

The inverse-degree/Cayley--Hamilton ratios are approximately

\[
0.0672,\qquad0.0292,\qquad0.0217,\qquad0.
\]

So the degrees 8, 21 and 26 are not the vacuous bounds 119, 719 and 1,199.
They are a genuine symmetry compression.  Because each primitive mass is
real symmetric, the minimal-polynomial degrees are also the exact numbers of
distinct eigenvalues.

The factor degrees over the integers are:

| form degree | irreducible factor degrees |
|---:|---|
| 0 | five linear, two quadratic |
| 1 | six linear, six quadratic, one quartic |
| 2 | five linear, eight quadratic, two cubic |
| 3 | one linear |

For example,

\[
\begin{split}
m_0(x)={}&(x-20)(x-11)(x-8)(x-6)(x-5)\\
&\times(x^2-22x+76)(x^2-20x+80).
\end{split}
\]

The complete exact coefficients and factorizations for all four blocks are
stored in the JSON certificate rather than copied into the note; the degree-1
and degree-2 coefficients are large.

## Why the certificate is exact

The verifier does not count numerically clustered eigenvalues.  For each
block it:

1. obtains an exact scalar-Krylov lower bound with Berlekamp--Massey over 12
   recorded primes;
2. reconstructs the monic integer polynomial by Chinese remaindering beyond
   an a priori coefficient bound;
3. evaluates the polynomial on the entire matrix modulo independent primes;
4. multiplies moduli beyond twice the integer residual bound, turning every
   modular zero into a certified integer zero;
5. applies the same bounded modular census to every entry of the inverse
   numerator.

For all four blocks the Krylov lower degree equals the whole-matrix
annihilator upper degree.  Minimality and the support census are therefore
exact; a probe that missed an eigenspace could not have passed the
whole-matrix test.

## What “finite” buys us

If

\[
m(x)=c_0+c_1x+\cdots+c_{s-1}x^{s-1}+x^s,
\]

then

\[
B^{-1}=-\frac{c_1I+c_2B+\cdots+B^{s-1}}{c_0}.
\]

Horner evaluation therefore gives exact upper bounds of 8, 21 and 26 local
mass applications for degrees 0, 1 and 2.  The nonzero support at graph
diameters 5, 10 and 15 gives corresponding lower bounds for any scheme whose
information advances by at most one mass-graph edge per step.

On this fixed complex the exact linear operation therefore has finite depth
windows

\[
5\leq T_0\leq8,qquad
10\leq T_1\leq21,qquad
15\leq T_2\leq26.
\]

This is the cleanest current answer to the “do we need infinity?” question:
**not on the fixed carrier, at the level of linear algebra.**

## What it does not buy us

The polynomial contains additions and large signed rational coefficients.  It
is not itself a unitary tick.  A reversible block encoding, reflection
product or ancilla construction is still required.  Moreover:

- the exact inverse is globally supported within each connected block;
- no refinement family has yet shown that degrees 8, 21 and 26 remain
  bounded as the carrier grows;
- if those degrees grow under refinement, the continuum inverse remains
  nonlocal even though every finite member has a finite polynomial;
- the construction selects no physical duration for one polynomial stage.

Thus the old infinite-penalty proposal is no longer the only algebraic route,
but the physical locality problem is not solved.

## Status ledger

- **DERIVED:** exact affine upper-Laplacian formulas for all four mass blocks.
- **DERIVED:** exact minimal-polynomial degrees `(9,22,27,1)`.
- **DERIVED:** exact inverse-polynomial degrees `(8,21,26,0)`.
- **DERIVED:** strong spectral compression relative to Cayley--Hamilton.
- **DERIVED NEGATIVE:** degrees 0, 1 and 2 inverses reach full graph diameter.
- **DERIVED:** no algebraic infinity is required on the fixed 600-cell.
- **STRUCTURAL:** interpreting polynomial evaluation as a multi-tick local
  process.
- **OPEN:** a coefficient-free reversible/unitary realization.
- **OPEN:** bounded polynomial degree and causal locality under refinement.
- **NOT CLAIMED:** physical time, mass, inertia, Lorentz invariance or the
  measured speed of light.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_mass_inverse_polynomial.py
```

Expected result: `10/10`.
