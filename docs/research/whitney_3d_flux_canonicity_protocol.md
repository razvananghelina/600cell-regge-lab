# Protocol: canonicity of a three-dimensional Whitney flux lift

Date: 2026-08-11

This protocol is committed before evaluating any of the new three-dimensional
averaging comparisons.  It uses no phenomenological target and performs no
low-spectrum calculation.

## Question and complete hypotheses

Let `K_k` be the closed three-complex obtained from the boundary of a
4-simplex by one barycentric subdivision followed by the canonical
rank-edgewise subdivision `Esd_k`, for `k=1,2`.  On every tetrahedron use the
exact Euclidean Whitney masses already fixed by the certified refinement
construction.

For degree `p`, let `C^p_loc` contain one copy of every global oriented
`p`-simplex for every incident tetrahedron, and let

\[
 J_p:C^p(K_k)\longrightarrow C^p_{\rm loc}(K_k)
\]

be the exact conforming copy injection.  The question is:

> Do combinatorial naturality, the supplied piecewise-Euclidean Whitney
> metric, the left-inverse condition `L_p J_p=I`, positivity of averaging
> weights, and strict occurrence locality select a unique `L_p` for
> `p=0,1,2`?

Here **strict occurrence locality** means that row `s` of `L_p` may depend
only on the local copies `(T,s)` of the same global simplex `s`.  This is the
bounded-star analogue actually used by the successful circle construction;
it is stronger than merely having finite support on each finite carrier.

The distinguished local basis is the oriented simplex basis.  Therefore a
diagonal entry of a Whitney mass is an intrinsic squared norm of that basis
cochain, not a freely chosen coordinate after the protocol.

## Frozen candidate class

For every global `p`-simplex `s`, let `O(s)` be its incident-tetrahedron
occurrences.  Three maps are evaluated.

### C: counting recovery

\[
 (L^C_px)_s={1\over |O(s)|}\sum_{(T,s)\in O(s)}x_{T,s}.
\]

### D: diagonal-Whitney recovery

If `q_{T,s}>0` is the diagonal entry of the exact local Whitney mass belonging
to `(T,s)`, define

\[
 (L^D_px)_s=
 {\sum_{(T,s)\in O(s)}q_{T,s}x_{T,s}
  \over\sum_{(T,s)\in O(s)}q_{T,s}}.
\]

Both rules are positive, strict-occurrence-local, invariant under relabelling
that preserves the simplicial metric, and reproduce constant copy data.  No
coefficient is fitted.  If they differ, then their normalized positive
powers already give further natural rules; the two frozen endpoints suffice
to refute uniqueness under the stated hypotheses.

### A: exact Hilbert-adjoint recovery

Let `M_loc,p` be the block-diagonal exact local Whitney metric and

\[
 M_{\rm conf,p}=J_p^*M_{\rm loc,p}J_p.
\]

Metric adjointness uniquely gives

\[
 L^A_p=M_{\rm conf,p}^{-1}J_p^*M_{\rm loc,p}.
\]

It is automatically a left inverse.  It is included to test the competing
possibility that exact adjointness, rather than strict locality, selects the
lift.  Its support is measured; locality is not assumed.

## Frozen computations

For every `k=1,2` and `p=0,1,2`, the verifier will:

1. reconstruct the exact carrier, copy injection and local Whitney blocks;
2. check the exact f-vector and positive local mass diagonals;
3. verify exactly `L^C_pJ_p=I` and `L^D_pJ_p=I`;
4. record the number of rows and coefficients on which `L^C_p` and `L^D_p`
   differ, their exact maximum difference, and an explicit first witness;
5. construct `L^A_p` numerically from the exact matrices, require a relative
   solve/left-inverse residual below `1e-10`, and record coefficients outside
   strict occurrence support using the frozen threshold `1e-11` times the
   maximum absolute coefficient;
6. record the maximum graph distance reached by those coefficients, where
   global `p`-simplices are adjacent when they occur in a common tetrahedron.

The numerical support statement for `L^A` is labelled **DERIVED NUMERICAL**,
not exact.  The exact comparison `L^C` versus `L^D` carries the logical
canonicity verdict.

Degree 3 is a positive control: every top simplex has one occurrence, so all
three recoveries must reduce to the identity.  It is not needed to lift a
forward block and is not included in the uniqueness decision.

## Frozen decision rules

- **DERIVED NEGATIVE FOR UNIQUENESS:** if `L^C_p != L^D_p` exactly for any
  frozen propagating degree and level.  Existing geometry and locality then
  admit at least two coefficient-free canonical constructions, so no 3D flux
  is selected and no spectral test will follow merely by choosing one.
- **DERIVED NUMERICAL LOCALITY CONFLICT:** if the certified `L^A_p` has any
  coefficient outside strict occurrence support.  Exact metric adjointness
  then selects a different, non-strictly-local operation on that control.
- **OPEN, NOT UNIQUE:** if `L^C=L^D` on both finite controls.  Agreement of
  two constructions is not a uniqueness theorem; the candidate class must
  then be enlarged before any positive selection claim.
- A **positive uniqueness result is forbidden** unless an algebraic theorem
  establishes uniqueness in the entire stated class.  Finite agreement alone
  cannot trigger it.

## Attack on the framing

The word "canonical" can mean either "defined without arbitrary labels" or
"uniquely selected by the axioms".  Many different natural transformations
can satisfy the first meaning.  Physics needs the second.  This protocol
therefore treats two distinct natural local maps as a negative, even if one
of them later gives an attractive spectrum.

Conversely, nonlocality of the exact Hilbert adjoint would not prove that no
other local flux exists.  It would prove only that **exact metric
adjointness plus strict occurrence locality** cannot both characterize that
particular recovery on the tested carrier.

## Scope exclusions

- no continuum eigenvalues, Betti modes, speed, mass, time or Planck units;
- no selection by spectral performance;
- no fitted mixture of the three maps;
- no claim about every conceivable higher-order star-local recovery;
- no full verifier-suite run, by explicit user request.

