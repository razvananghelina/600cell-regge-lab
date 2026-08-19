# Result: pure corrected struts have zero canonical intersection

Date: 2026-08-19

## 1. Provenance ledger

| stage | commit | status |
|---|---:|---|
| prior-art gate | `27b8e85` | exact graph-intersection object frozen |
| primary protocol | `786cae7` | no actual nullity inspected |
| primary verifier | `08b5237` | registered before execution |
| primary artifact | `b64fd83` | 14/14, all 14 nullities zero |
| adversarial protocol | `2afb0c7` | zero target explicitly disclosed |
| adversarial verifier | `b5739c8` | QR/Frobenius route registered before execution |
| first adversarial execution failure | `1438df3` | JSON NumPy-boolean failure preserved; no artifact |
| serialization-only repair | `e6fcb66` | no scientific matrix or threshold changed |
| adversarial artifact | `6236103` | 9/9, zero intersection corroborated |

Frozen artifact hashes:

```text
primary intersection census
422d8d8cb0fc0d72d842e3bf79609d4d985da6237c58e7c699b5f9cc21b65cec

adversarial intersection audit
c186260ee9520eac59658e3290fb1f4502fd9a7d92f533e8774506cd30e9d03b
```

## 2. Mathematical reduction

In a minimal `2T` sector of irrep dimension `d`, the corrected geometric
strut graph and canonical strong-equation graph are

```text
G,C : C^(5d) -> C^(65d).
```

Both have the same literal identity on the same five pole positions, in the
same geometry-derived coefficient order.  Therefore

```text
im(G) intersect im(C)  is isomorphic to  ker(G-C).
```

This reduction was checked before reporting a nullity.  The independent
stacked-image formula

```text
dim(im G intersect im C) = rank(G)+rank(C)-rank([G,C])
```

was also required in every sector.

## 3. Primary result

All primary controls passed 14/14.  The calibrated nullities, separately for
both staircase parities, were

```text
even: [0,0,0,0,0,0,0]
odd:  [0,0,0,0,0,0,0].
```

The four derivative variants agreed in every case.  Synthetic `G=C` controls
returned nullity `5d`; embedded full-rank differences returned nullity zero;
the common coefficient-basis transform preserved every nullity; and the
stacked-image formula reproduced all 14 results.

The smallest operational singular values were:

```text
six non-homogeneous irrep sectors  6.569739702 .. 15.033173993
homogeneous irrep sector           1.173190442e-5
```

The most delicate singular value was still `2.792e6` times its complete
matrix-uncertainty bound.  It is small relative to the other sectors, but not
numerically compatible with zero.

## 4. Adversarial calculation

The adversarial verifier did not call the primary intersection verifier or
use an SVD/eigendecomposition to decide rank.  For a column-pivoted QR

```text
D P = Q R
```

it used

```text
sigma_min(D) >= 1/||R^-1||_F
```

from a triangular solve, then subtracted the independently reconstructed
matrix uncertainty.  All four variants in all 14 sectors retained a robust
positive lower bound.  The smallest adversarial lower-bound-to-uncertainty
ratio was `8.536e5`.

All 9 controls passed, including zero/full-rank synthetic blocks, coefficient
mixing, row reversal with alternating phase, complex conjugation, QR residuals
and the frozen source/target-role corruption.  The corruption changed a raw
QR lower bound by as much as `1.831e-2`.

Accepted outcome:

```text
CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_CORROBORATED
```

## 5. Global statement and limits

The five pole orbits carry five copies of the regular `2T` representation.
Each irrep block of dimension `d` repeats with its regular multiplicity `d`.
Thus full column rank in every minimal block implies full global rank 120 for
each staircase parity.

**DERIVED COMPUTATIONAL, adversarially corroborated under the frozen slab and
action hypotheses:**

```text
im(G_corrected pure strut) intersect im(C_canonical) = {0}.
```

Equivalently, no nonzero pure-strut coefficient vector has a corrected
geometric response that also satisfies the frozen canonical strong equations.

This statement means:

- pure struts are kinematically admissible face-gluing data;
- they are not independent canonical lapse freedom on this curved slab;
- any nonzero canonical direction inside the complete accepted kinematic
  carrier must include non-strut data, hence a scale/form component.

It does **not** mean that struts are absent from evolution.  It means they
cannot evolve alone.  It also does not prove gauge freedom, a graviton,
stability, a continuum dispersion relation or a physical clock.

## 6. Post-result literature check

The post-result search added the explicit distinction between kinematical and
dynamical constraints to the terms already used for Regge pre/post and
pseudo-constraints.

- Khatsymovsky,
  [*On kinematical constraints in Regge calculus*,
  arXiv:gr-qc/9311005](https://arxiv.org/abs/gr-qc/9311005), studies how
  kinematical relations among Regge variables interact with equations of
  motion in a continuous-time Hamiltonian formulation.
- Khatsymovsky,
  [*Regge calculus in the canonical form*,
  arXiv:gr-qc/9310004](https://arxiv.org/abs/gr-qc/9310004), classifies
  constraints in a different connection/area-variable canonical model.
- Bahr and Dittrich,
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670), show that exact gauge
  symmetries are generically broken on curved Regge solutions and replaced by
  pseudo-constraints.
- Dittrich and Hoehn,
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974), derive simplicial
  canonical evolution and pre/post constraints from the action.

These sources support the general warning that kinematical admissibility is
not canonical freedom.  None of the located primary sources contains this
600-cell strut graph, its action graph or the zero-intersection calculation.
External novelty remains **OPEN**; a search is not a novelty proof.

## 7. Consequence for the programme

This is a useful negative and a genuine narrowing of the dynamics.  The old
idea that one could interpret 120 arbitrary struts directly as 120 lapse
directions is closed on the fixed curved slab.  The action forces coupling to
the other 120 directions of the exact kinematic carrier.

The next calculation is therefore not another pure-strut comparison.  It is
the target-blind assembly of the complete 240-column dynamic
scale-plus-strut carrier, followed by its pullback through the frozen action
Hessian/strong-equation map.  The first outputs must be rank, nullity and
sector multiplicities, before comparison with continuum scalar/vector/tensor
labels.

No full suite was run.  The scoped registry after the final registration had
352 distinct verifier names, 354 files including two documented exclusions,
zero duplicates, zero unregistered files and zero missing files.

