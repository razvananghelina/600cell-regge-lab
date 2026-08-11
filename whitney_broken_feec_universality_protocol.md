# Protocol: do the two broken-FEEC projections share a strong-limit spectrum?

Date: 2026-08-11

This protocol is committed before evaluating any `k=4` projection weights or
any strong-penalty constrained spectrum.  No phenomenological or smooth
manifold eigenvalue is used as a target.

## Correction that fixes the question

The weak control `alpha=1` in commit `8fce4c9` has the correct exact kernel,
but it is not in the spectral-convergence regime of the published
broken-FEEC theorem.  Campos Pinto--Güçlü require

\[
 \alpha_h\geq C h^{-s}
\]

for an `s`-regular domain, and their numerical section exhibits spurious low
modes at `alpha=1`.  Therefore this protocol will not track the lowest weak
penalty eigenvalue.

Instead it studies the coefficient-free strong-penalty limit at each fixed
mesh.  If

\[
 K^X_{h,\alpha}=K^X_{h,0}
 +\alpha(I-P^X)^TM(I-P^X),
\]

then every eigenbranch that remains finite as `alpha -> infinity` lies in
`im(J)=ker(I-P^X)`.  Its exact constrained pencil is

\[
 \widehat K^X_{h,p}=J_p^T K^X_{h,0,p}J_p,
 \qquad
 \widehat M_{h,p}=J_p^TM_{h,p}J_p.
\]

No stabilization coefficient remains to be chosen.

## Frozen carrier and candidates

Use

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

and the two already preregistered projections:

1. equal counting `P^C=J L^C`;
2. diagonal-Whitney/mass-lumped `P^D=J L^D`.

The `k=1,2` calculations calibrate the new strong-limit implementation.
Only `k=4` is new evidence.

## Exact microscopic audit

For `k=2,4`, degrees `p=0,1,2,3`, record:

- rows and coefficients on which `L^C_p` and `L^D_p` differ;
- the exact maximum absolute weight difference;
- the fraction of global rows that differ.

If a nonzero exact maximum at `k=2` repeats without decrease at `k=4`, report
a **DERIVED NEGATIVE FOR UNIFORM MICROSCOPIC CONVERGENCE**: the two projection
matrices do not converge coefficientwise in supremum norm.  This statement is
about the full broken carrier, not about their low-energy restrictions.

## Frozen strong-limit spectral audit

For every level, candidate and form degree:

1. construct `d_h=D_pw P` and its broken-metric adjoint;
2. construct `J^T K_{h,0} J` directly, without a large finite penalty;
3. compute the lowest eight generalized eigenvalues at `k=1,2` and the lowest
   six at `k=4`;
4. retain the exact Betti zero multiplicities `(1,0,0,1)` and compare the
   first four positive ordered eigenvalues;
5. require every reported Ritz pair to have relative residual below `1e-7`;
6. use deterministic starting blocks and no target-dependent branch choice.

At `k=1`, the two pencils must agree because the projections already agree;
this is the calibration gate.  At `k=2,4`, record the maximum relative
difference of the four ordered positive values in each degree.

The `k=4` solver may use sparse shift-invert or block LOBPCG.  A solver that
does not meet the frozen residual is a recorded **OPEN COMPUTATIONAL
FAILURE**, not evidence for agreement or disagreement.

## Frozen interpretation

- If the maximum relative difference decreases from `k=2` to `k=4` in every
  degree `p=1,2,3`, label this only **PATTERN TOWARD A COMMON LOW-ENERGY
  CLASS**.  Two comparisons cannot prove a continuum theorem.
- If it fails to decrease in any degree, label **PATTERN NEGATIVE FOR COMMON
  LOW-ENERGY FLOW**.  One failed step does not prove distinct continuum
  limits.
- Exact equality at `k=4` would be **DERIVED ON THE CONTROL**, not an all-level
  theorem.
- The continuum question remains **OPEN** unless the analytic approximation
  hypotheses for both `P` and `P*` are proved.  Numerical flow cannot replace
  that proof.

## Attack on the framing

The strong-penalty constrained spectrum is mathematically clean but is not a
derived microscopic tick: it is a limit that removes mismatch modes.  A
common constrained spectrum would establish robustness of a low-energy
effective operator, not uniqueness of the finite dynamics.

Conversely, persistent coefficient differences on the full broken carrier do
not by themselves refute common low-energy universality.  The protocol keeps
these two statements separate.

## Scope exclusions

- no fit of an exponent or limiting value;
- no use of the weak `alpha=1` lowest branch as physical evidence;
- no selection between the two candidates;
- no claim about a round `S^3` spectrum—the fixed limit geometry is the
  piecewise-flat boundary-of-simplex control;
- no Lorentzian time, causal speed, inertia, mass or Planck units;
- no full suite run, by explicit user request.

