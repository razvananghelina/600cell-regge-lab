# Fixed finite trace stiffness fails the known topology of the circle

Date: 2026-08-11

Preregistration commit: `d499ca0`

Targeted verifier:
`reproducible/verify_whitney_finite_stiffness_circle.py`

Targeted result: **9/9 PASS**.  The verifier is registered exactly once.  The
full suite was not run by explicit user request.

The first execution stopped before producing data because SymPy represented
the same radical in a factored form and the installed NumPy requires lists,
not tuples, in `numpy.block`.  The correction changed only those two code
comparisons.  The preregistered operator, grids, stiffness laws, and gates were
unchanged.

## Headline

Uniform separation of the conforming and mismatch sectors is not sufficient
for a correct continuum limit.

On the known-answer unit circle:

> **DERIVED NEGATIVE:** the local-weak-plus-pure-jump pencil at every fixed
> finite dimensionless stiffness loses the harmonic Whitney 1-form.  Its
> spurious eigenvalue diverges as the mesh is refined.

Even scaling stiffness proportionally to `1/h` is insufficient:

> **DERIVED NEGATIVE:** `kappa h -> constant` leaves a nonzero spurious
> topological gap.

Recovery requires

\[
 \kappa h\longrightarrow\infty,
\]

or a different operator containing consistency/flux terms.

This reaches the kill boundary for a **fixed-finite-`kappa` pure penalty** as a
general conforming continuum mechanism.  It does not kill geometry-selected
DG flux completions.

## Exact circle construction

Divide a unit circle into `N` equal edges, `h=1/N`.  Each edge retains local
coefficients `(a_n,b_n,e_n)` for its left vertex, right vertex, and 1-form.
The exact Whitney metric is

\[
 M_{0,h}=\frac h6
 \begin{pmatrix}2&1\\1&2\end{pmatrix},
 \qquad
 M_{1,h}=\frac1h.
\]

The local exterior derivative is `(-1,1)`, and the jump at vertex `n` is

\[
 a_n-b_{n-1}.
\]

The tested pencil contains exactly the terms already proposed in three
dimensions:

\[
 (W_h+\kappa R_h^*R_h)v=zM_hv.
\]

No continuum target enters its construction.  The target used only for the
gate is topological: the conforming circle has

\[
 (b_0,b_1)=(1,1),
\]

so its Kähler--Dirac operator has two exact zero modes at every finite `N`.

## Exact Bloch calculation

In the translation-invariant `q=0` sector, write

\[
 c=\frac{a+b}{\sqrt2},
 \qquad
 m=\frac{a-b}{\sqrt2}.
\]

The scalar `c` remains an exact zero.  After mass orthonormalization, the
remaining mismatch/1-form block is

\[
 hA_0=
 \begin{pmatrix}
 12\kappa&-\sqrt{12}\\
 -\sqrt{12}&0
 \end{pmatrix}.
\]

Its exact eigenvalues are

\[
 hz_\pm=6\kappa\pm\sqrt{36\kappa^2+12}.
\]

The negative branch is the would-be harmonic 1-form.  Rationalizing gives its
absolute spectral gap:

\[
 |z_-|
 =\frac{12}
 {h\left(\sqrt{36\kappa^2+12}+6\kappa\right)}.
\]

For large `kappa`,

\[
 |z_-|\sim\frac1{\kappa h}.
\]

This is an exact asymptotic, not a fit.

## Three frozen scaling laws

The verifier evaluated `N=8,16,32,64,128`.

### Fixed `kappa=1`

| `N` | spurious harmonic gap |
|---:|---:|
| 8 | 7.42563 |
| 16 | 14.8513 |
| 32 | 29.7025 |
| 64 | 59.4050 |
| 128 | 118.8100 |

The value doubles exactly with refinement and its symbolic limit is infinity.
The missing harmonic mode does not merely fail to improve; it is driven to the
cutoff.

### `kappa=N=1/h`

The gaps are

\[
 0.998701, 0.999675, 0.999919, 0.999980, 0.999995,
\]

with exact limit

\[
 |z_-|\longrightarrow1.
\]

More generally, `kappa=c/h` leaves the nonzero limit `1/c`.  Merely matching
the inverse mesh scale is therefore not enough.

### `kappa=N^2=1/h^2`

The gaps are

\[
 0.124997, 0.0624999, 0.0312500, 0.0156250, 0.0078125,
\]

and tend to zero.  Here

\[
 \kappa h=N\longrightarrow\infty.
\]

This recovers the harmonic mode only through a singular double scaling.  The
positive stiff branch simultaneously grows as approximately `12 N^3` in
this family.  That divergence is a spectral fact; it is not called an energy
or speed because no physical time unit has been derived.

## Independent full-matrix control

For `N=8,16` and `kappa=1,2,4`, the verifier constructs every complete
`3N x 3N` generalized pencil without Fourier reduction.

- every pencil is Hermitian with positive metric and semipositive penalty;
- every finite-`kappa` local pencil has exactly **one** zero mode;
- the two nonzero `q=0` eigenvalues occur in the full spectrum;
- maximum relative disagreement with the exact formula is `2.84e-14`;
- independently assembled conforming Whitney matrices have exactly **two**
  zero modes.

Thus the failure is not an artefact of selecting a Bloch subblock.

## Why the harmonic mode is lifted

A globally constant 1-form is harmonic because the codifferential
contributions from adjacent edges cancel at every assembled vertex.  In the
duplicated local operator those contributions live in different vertex
copies and do not cancel.

The pure jump term penalizes the copy mismatch but supplies no consistency
flux that performs the missing cancellation.  At finite stiffness the
1-form mixes with the mismatch coordinate and moves away from zero.

This is exactly the distinction exposed by the previous role audit:

- the penalty separates mismatch directions;
- it has zero tangent action on conforming fields;
- separation alone does not make the local weak operator consistent on the
  conforming topology.

## Consequence for the three-dimensional route

The circle is a counterexample to the general proposition that a fixed finite
`kappa` pure penalty automatically converges to conforming Kähler--Dirac
dynamics.  Therefore that proposition is closed.

The result does **not** prove the numerical value of the lifted 3-sphere top
form.  But the 3-sphere also has two Kähler--Dirac harmonic modes,

\[
 b_0=b_3=1,
\]

and its top-form cancellation likewise occurs through assembly.  A
three-dimensional finite-`kappa` claim now carries the burden of explicitly
showing that the top harmonic survives; sector separation is no longer
acceptable evidence.

## What remains alive

Two routes remain mathematically distinct:

1. **Singular double scaling:** choose `kappa_h` with
   `kappa_h h -> infinity`.  This recovers topology but revives the divergent
   microscopic branch and still selects no law for `kappa_h`.
2. **Consistency/flux completion:** add a face-local term selected by
   adjointness and geometry so that assembled smooth/cohomological modes are
   already consistent at finite stiffness.  This is the mechanism used in
   genuine DG discretizations, but its precise Kähler--Dirac version is not
   yet derived here.

The second route is the only one that could preserve a finite dimensionless
stiffness without the circle obstruction.  An arbitrary flux coefficient
would be fitting.

## Status ledger

- **DERIVED:** exact `q=0` block and eigenvalues.
- **DERIVED:** complete spectra reproduce the Bloch result.
- **DERIVED NEGATIVE:** fixed finite `kappa` destroys one Betti zero mode.
- **DERIVED NEGATIVE:** `kappa proportional to 1/h` leaves a nonzero gap.
- **DERIVED:** recovery requires `kappa h -> infinity` for the pure penalty.
- **STRUCTURAL NEGATIVE:** uniform mismatch-sector separation does not imply
  correct low-energy topology.
- **OPEN:** a canonical Kähler--Dirac consistency flux.
- **OPEN:** the corresponding three-dimensional harmonic and continuum
  spectra.
- **NOT CLAIMED:** particle mass, inertia, time, causal speed, or Planck
  scale.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_finite_stiffness_circle.py
```

Expected result: `9/9`.

