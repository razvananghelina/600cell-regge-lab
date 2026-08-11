# Weighted incidence coins cannot repair the Whitney metric

Date: 2026-08-11  
Preregistration commit: `ab6c0de`

## Result

Allowing arbitrary complex amplitudes and reversal phases on the existing
Hasse-incidence arcs does not repair the metric mismatch of the local
Kähler--Dirac tick.

> **DERIVED NO-GO FOR WEIGHTED INCIDENCE COINS:** every rank-one weighted
> Grover--Szegedy walk on the original directed-incidence carrier has zero
> discriminant entries away from Hasse incidence, whereas the exact Whitney
> metric requires 22 nonzero off-incidence entries.

The targeted verifier passes `10/10` in about 0.1 seconds.  No full suite was
run.

## Complete class covered

For every simplex \(x\), permit

\[
|s_x\rangle=\sum_{y\sim x}a_{xy}|x,y\rangle,
\qquad
\sum_{y\sim x}|a_{xy}|^2=1,
\]

with arbitrary unequal, zero or complex amplitudes.  Also permit arbitrary
involutive reversal phases

\[
S|x,y\rangle=\omega_{xy}|y,x\rangle.
\]

The coin and discriminant are

\[
C=2\sum_x|s_x\rangle\langle s_x|-I,
\qquad
T=A^*SA.
\]

This covers every rank-one weighted Szegedy coin on exactly one directed
state per codimension-one incidence.  It does not cover added arcs, ancillas,
nonorthogonal Whitney encodings or multi-step effective operators.

## Universal support obstruction

The exact matrix element is

\[
T_{xy}
=\overline{a_{xy}}\,\omega_{yx}\,a_{yx}
\]

when \(x\sim y\), and

\[
T_{xy}=0
\]

when no such incidence arc exists.  Changing amplitudes, phases, diagonal
basis gauges or an overall normalization can alter values but never create a
missing matrix position.

For the 15 simplices of one tetrahedron, the carrier has 28 undirected or 56
directed incidences.  The symbolic support of the complete weighted class is
exactly those 56 positions.

## Exact Whitney contradiction

On the barycentric flag child, direct rational integration gives

\[
\delta_k=M_k^{-1}d_k^TM_{k+1}.
\]

The numbers of its nonzero entries lying outside \(d_k^T\) incidence are

\[
(N_0,N_1,N_2)=(10,12,0).
\]

Therefore 22 entries required by the metric operator vanish identically for
every walk in the preregistered weighted class.  This is a support theorem,
not a failed parameter search.

Every one of those 22 simplex pairs is at exact Hasse distance three.  Hence
any rescue that still propagates only across Hasse incidences requires at
least three local substeps before it can even have the necessary matrix
support.  Depth three is necessary, not sufficient.

## Physical meaning

- **DERIVED:** simple reweighting cannot turn the exact graph tick into the
  metric Whitney dynamics.
- **DERIVED:** one effective metric step requires incidence depth at least
  three on the local child, unless the carrier is enlarged.
- **OPEN:** a canonical three-or-more-substep factorization.
- **OPEN:** a local dilation with auxiliary element states or a
  nonorthogonal Whitney encoding.
- **NOT CLAIMED:** impossibility of all local quantum dynamics.

The finite-speed idea survives, but the metric cannot be inserted merely by
choosing clever probabilities.  Some additional mechanism must store or
mediate the non-diagonal mass information.  That mechanism is now the place
where any genuine new physics would have to enter.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_weighted_szegedy_metric_nogo.py
```

Expected result: `10/10`.

