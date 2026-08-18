# Local Kähler--Dirac tick: refinement verdict

Date: 2026-08-11  
Preregistered protocol: commit `59ab096`

## Result

The unweighted signed Grover--Szegedy tick has an **exact dyadic refinement
law on the circle**, but it is **not the full metric Whitney dynamics on a
barycentric tetrahedron**.

The honest final label is:

> **DERIVED KINEMATIC ONLY:** the construction gives a valid local unitary
> clock and cone for the abstract incidence graph, but not yet for the full
> piecewise-Euclidean spectral geometry used by the theory.

The targeted verifier passes `17/17` in under one second.  No full suite was
run.

## Exact circle control

For an oriented circle with \(N\) edges, the positive singular value of the
normalized signed incidence at Fourier mode \(m\) is

\[
\lambda_N(m)=\sin\frac{\pi m}{N}.
\]

Below and including Nyquist, the walk quasienergy relative to its central
phase is therefore

\[
\varepsilon_N(m)
=\arcsin\lambda_N(m)
=\frac{\pi m}{N}.
\]

Hence, exactly and without a fitted slope,

\[
\varepsilon_N(m)=2\varepsilon_{2N}(m).
\]

Two fine micro-ticks reproduce the phase interval of one coarse micro-tick
for every resolved mode tested at \(N=8,16,32,64\).  This is a genuine
refinement result.

If the dimensionless tick duration is set equal to the Hasse-edge length
\(1/(2N)\), then

\[
\frac{\varepsilon_N(m)}{1/(2N)}=2\pi m.
\]

This fixes unit speed in lattice units.  It does **not** determine seconds or
metres; the overall conversion factor remains free.

## Tetrahedral metric test

One barycentric flag child of the regular reference tetrahedron has squared
face areas

\[
\left(\frac13,\frac12,\frac29,\frac1{18}\right),
\]

so it is strongly anisotropic.  Its volume is exactly \(1/24\) of the
parent.

The preregistered test examined the metric codifferential from 3-forms to
2-forms,

\[
\delta_2=M_2^{-1}d_2^T M_3.
\]

Contrary to the initial suspicion, this direction remains exactly
proportional to signed incidence:

\[
\delta_2=\frac{540}{29},d_2^T.
\]

Therefore the literal preregistered gate passed.  This is a small positive
structural result and is not discarded because it was unexpected.

## Why that pass was insufficient

The preregistration tested only the top-to-face block.  A full Kähler--Dirac
operator also contains the \(0\leftrightarrow1\) and
\(1\leftrightarrow2\) blocks.  The post-protocol hostile audit checked the
necessary proportionality in all degrees.

On the regular parent tetrahedron, symmetry gives the three exact ratios

\[
\left(\frac54,\frac52,\frac{15}{4}\right).
\]

On the anisotropic flag child, the lower two proportionalities fail:

\[
(c_0,c_1,c_2)=\left(\text{none},\text{none},\frac{540}{29}\right).
\]

More strongly, the exact metric adjoints create matrix elements outside the
corresponding incidence supports:

\[
N_{\rm off-incidence}=(10,12,0)
\]

for degrees \(0,1,2\).  The consistent Whitney 2-form mass and its inverse
are both non-diagonal.

Thus the unweighted one-incidence-per-tick walk cannot equal the accepted
metric Whitney Kähler--Dirac operator in the lower degrees.  The exact cone
belongs to the abstract Hasse graph metric.

## Methodological verdict

- **DERIVED:** exact dyadic quasienergy scaling on the circle.
- **DERIVED:** exact unit limiting speed in lattice units, conditional on the
  convention `tick duration = Hasse-edge length`.
- **DERIVED:** the top \(3\to2\) Whitney direction remains incidence-uniform
  on the barycentric child.
- **DERIVED NEGATIVE:** the lower metric adjoints are not
  incidence-proportional and contain 10 and 12 off-incidence entries.
- **DERIVED NEGATIVE:** the unweighted local tick is not the complete metric
  Whitney dynamics.
- **OPEN:** a canonical metric-aware unitary dilation.
- **NOT DERIVED:** Lorentzian time, SI units, \(c\), \(G\), \(\hbar\), Planck
  time or Planck mass.

The protocol's original positive label was too broad.  Its stated test did
pass, but one passing block did not cover the full operator.  The result is
therefore recorded separately as a protocol pass and a framing correction,
not silently promoted to a spacetime result.

## What remains possible

This closes only the **unweighted** walk as the complete metric dynamics.  It
does not prove that every metric-aware local unitary dilation is impossible.
The two obvious choices already exhibit the known tradeoff:

- consistent Whitney mass: exact Galerkin induction, but a nonlocal strong
  adjoint and no strict one-step cone;
- diagonal mass lumping: local finite-speed dispersion, but loss of exact
  refinement isometry.

A continuation must either prove a no-go theorem for satisfying both, or
construct a coefficient-free unitary dilation of the local weak pair
\((A,M)\) without hiding the dense inverse in nonlocal transitions.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_kahler_dirac_tick_refinement.py
```

Expected result: `17/17`.

