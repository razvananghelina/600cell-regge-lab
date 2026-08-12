# Preregistration: can a strict-local recovery be a cochain map?

Date: 2026-08-12

This protocol is committed before evaluating the new exhaustive support
census or either candidate's commutator.  No spectrum, continuum target or
phenomenological number is used.

## Complete question

On the closed controls

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}\partial\Delta^4),
 \qquad k=1,2,4,
\]

let `V_p` be the direct sum of local simplicial `p`-cochains over all
tetrahedra, `W_p` the global conforming cochains, `J_p:W_p -> V_p` the signed
occurrence injection, `D_pw,p` the blockwise simplicial coboundary, and `d_p`
the global simplicial coboundary.

Consider every strict-occurrence-local recovery

\[
 (L_px)_s=\sum_{T\supset s}\epsilon(T,s)w_{T,s}x_{T,s},
 \qquad \sum_{T\supset s}w_{T,s}=1.
\]

The signs use the same increasing-simplex orientation convention as `J`.
The sum condition is precisely `L_pJ_p=I`; weights are otherwise arbitrary
real numbers.  In particular, positivity, metric weighting, symmetry and the
two previously tested formulas are *not* assumed.

Ask whether the recoveries can form a cochain retraction:

\[
 L_{p+1}D_{{\rm pw},p}=d_pL_p,
 \qquad p=0,1,2.
\]

Equivalently, because `D_pw J = J d` and `J` is injective, the associated
projection `P=JL` commutes with the raw piecewise coboundary.

This condition is independently meaningful: recovery followed by the global
exterior derivative would equal the recovery of the elementwise exterior
derivative.  It is the standard cochain-map requirement, not a condition
chosen from a preferred spectrum.

## Frozen exact support obstruction

Fix one local occurrence `(T,s)` in degree `p`.  The corresponding column of
`L_{p+1}D_pw,p` can have a nonzero row only at a global `(p+1)`-simplex `u`
with

\[
 s\subset u\subset T.
\]

The same column of `d_pL_p` has coefficient

\[
 \pm w_{T,s}
\]

at *every* global coface `u` of `s`.  Hence, if there is a global coface
`u` not contained in `T`, the commuting equation at that row forces
`w_{T,s}=0` exactly.

The verifier will exhaust all occurrences at all three levels and record,
for each degree:

- total occurrence count;
- number and fraction having an external global coface;
- minimum and maximum number of external cofaces;
- whether the resulting forced-zero weights hit every occurrence and hence
  contradict every left-inverse sum.

No linear solver or numerical rank tolerance is required for this theorem.
The certificate is an explicit combinatorial witness for every unknown.

## Candidate controls

Independently of the all-weight theorem, rebuild the equal-counting and
diagonal-Whitney recoveries at `k=1,2` and compute

\[
 C^X_p=L^X_{p+1}D_{{\rm pw},p}-d_pL^X_p,
 \qquad X\in\{C,D\}.
\]

Record exact support counts and maximum absolute coefficients.  These
candidate commutators are controls only; their failure cannot establish the
general theorem without the support argument above.

## Frozen decisions

- If every occurrence in any degree has an external coface, all its weights
  are forced to zero, contradicting `LJ=I`.  Report a **DERIVED NO-GO IN THAT
  DEGREE** for every strict-occurrence-local cochain retraction.
- If this holds in all degrees `p=0,1,2`, report a **DERIVED NO-GO FOR A
  STRICT-OCCURRENCE-LOCAL COMMUTING RECOVERY** on the frozen controls.
- A repeated all-level combinatorial pattern is not by itself an arbitrary
  triangulation theorem.  A general closed-manifold corollary may be stated
  only with the explicit hypothesis that every local occurrence has an
  external global coface.
- If an occurrence lacks an external coface, the support proof does not
  force its weight to zero; solve the remaining exact linear system rather
  than inferring a verdict.
- No candidate will be selected because its commutator is smaller.

## Attack on the framing

This requirement may be too strong for a discontinuous carrier.  Broken FEEC
defines the usable differential as

\[
 d_h=D_{\rm pw}P,
\]

precisely instead of assuming that a recovery commutes with the raw
piecewise derivative.  Therefore a no-go would close one proposed local
selection axiom, not refute broken FEEC or all commuting projections.

In particular, standard smoothing or macroelement projections may spread a
single broken coefficient over a bounded neighbourhood.  Such a recovery is
not strict-occurrence-local and is outside this protocol.  Enlarging the
support is the logically indicated next route if the frozen test fails, but
its radius and construction must be fixed before any spectrum is inspected.

## Outputs and exclusions

The registered verifier will write
`reproducible/whitney_commuting_recovery.json`.  The result note will preserve
the complete hypotheses and distinguish exact no-go, structural scope and
open macro-star alternatives.

Excluded:

- no spectrum or fitted tolerance;
- no positivity assumption on the general weights;
- no Hodge-star or dual mesh;
- no claim about time, causality, mass or Planck units;
- no full-suite run, by explicit user request.
