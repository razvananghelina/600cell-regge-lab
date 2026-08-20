# Prior-art gate: exact homogeneous scale--strut/canonical line

Date: 2026-08-20  
Status: **prior-art gate; no new result accepted here**

## Exact object and hypotheses

Work on the fixed nonstatic Lorentzian 600-cell slab already frozen by the
complete-carrier mission. The lower spatial length is `L0`, the upper/lower
ratio is real `lambda != 1`, the timelike pole magnitude is `rho > 0`, and
`q_diag=lambda*L0^2-rho > 0`. The full scale--strut carrier and the action
canonical graph retain their frozen source/target, parity, logarithmic-strut and
weak-pole conventions.

The only unresolved representation sector is the unique constant-overlap
dimension-one sector. Its reduced matrices are

```text
D = [G_scale, G_strut-C] : C^10 -> C^65,
K = [G_scale, G_strut,-C] : C^15 -> C^65.
```

The primary interval calculation certifies ranks at least 9 and 14, respectively,
but its P100-frozen candidate did not pass the P160 no-refit zero gate. The new
question is whether the homogeneous canonical graph itself supplies an exact
nonzero kernel vector, not whether a high-precision singular value looks small.

## Primary literature

- Dittrich and Hoehn formulate action-generated canonical simplicial evolution
  and weak/free data in Pachner evolution:
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- Hoehn derives canonical linearized Regge constraints and lattice-graviton
  counting on flat backgrounds:
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672).
- Bahr and Dittrich show that curved Regge backgrounds generically replace exact
  gauge constraints by pseudo-constraints:
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670).
- De Felice and Fabri study generalized 600-cell evolution:
  [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).
- Liu and Williams derive homogeneous and perturbed closed Regge lattice-universe
  evolution: [arXiv:1502.03000](https://arxiv.org/abs/1502.03000).

## KNOWN / CONTROL / OPEN

- **KNOWN:** an action serves as a generating function for canonical simplicial
  evolution; weak variables and pseudo-constraints are standard features.
- **KNOWN:** homogeneous 600-cell/closed-lattice Regge cosmologies and
  less-symmetric perturbations are prior art.
- **CONTROL:** the repository already has an independently derived closed
  homogeneous frustum action and exact flat-subdivision identity.
- **CONTROL:** the complete primary resolver certifies rank at least 9 for D and
  14 for K, so one constructed null line would fix their nullities exactly.
- **OPEN:** whether the present full carrier and weak-pole canonical graph share
  exactly one homogeneous line.
- **OPEN:** whether such a line also satisfies the omitted weak/pole equation,
  represents gauge, or gives physical evolution.
- **OPEN external novelty:** no searched primary source contains this exact
  `1560 x 240` carrier, weak-pole graph and binary-symmetry sector reduction.
  Search absence is not proof of novelty.

## Framing attack

An exact line in the weak-pole canonical graph would not by itself be a physical
mode. Because the pole equation is deliberately omitted when the lapse is treated
as weak/free data, the line may be an off-shell lapse response. Physical
interpretation is forbidden until that omitted equation is tested separately.
