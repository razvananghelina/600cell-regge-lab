# Primary multiprecision carrier/action intersection result

Date: 2026-08-20  
Status: **DERIVED primary result; adversarial replication pending**

## Frozen execution

The verifier source was committed in `7fdbd95` after the precision-aware repair
protocol was preregistered in `c4a6d74`.

- verifier SHA-256:
  `b836cad394fe8a54644d514b6f31cb899ff5c3697b6c2a5f4edfc2b5f0ac5d62`;
- artifact SHA-256:
  `75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706`;
- targeted execution: `17/17`;
- outcome: `FULL_SCALE_STRUT_CANONICAL_HOMOGENEOUS_OPEN`.

No full-suite execution was performed, by explicit project instruction.

## Precision-aware geometry control

Both P200G parity audits pass every preregistered condition. They use the same
finite-difference steps as P160 and do not rebuild or refit the scientific
intersection classifier.

| parity | P160 imaginary | P200G imaginary | ratio P200G/P160 | max real-entry change |
|---|---:|---:|---:|---:|
| even | `3.5680071e-119` | `2.3179983e-159` | `6.4966191e-41` | `4.6876921e-116` |
| odd | `7.4154131e-119` | `8.0283531e-159` | `1.0826576e-40` | `6.4032880e-116` |

Thus the earlier P160 imaginary-residue failure is **DERIVED** to be numerical
contamination under the preregistered convergence test, rather than a geometry
failure.

## Primary sector classification

### Nonhomogeneous sectors

- **DERIVED primary:** every nonhomogeneous D and K interval Gram determinant
  excludes zero at P100 and P160.
- The binary classifier does not disagree.
- The minimum singular-value / interval-ball-radius ratio is
  `2.8204006836e132`.
- The maximum P100/P160 relative singular-value change is
  `8.3424508240e-35`.

This is strong primary evidence that the canonical carrier/action intersection is
zero in every nonhomogeneous sector. It is not consolidated until a mechanically
different adversarial verifier agrees.

### Homogeneous sector

- **OPEN:** D has a certified rank-at-least-9 minor and K a certified
  rank-at-least-14 minor in each parity.
- Exactly one P160 midpoint singular value lies below `1e-50` in D and K.
- The next singular values are approximately `9.8234441621e-8` for D and
  `8.0184489625e-8` for K.
- The P100 candidate, frozen to 70 digits and evaluated at P160 without refitting,
  has residual approximately `2.4898227e-42`, failing the preregistered `<1e-50`
  exact-zero gate.

Therefore neither a zero-dimensional nor a one-dimensional homogeneous
intersection is derived. The honest classification is `OPEN`.

## Physical scope

This calculation does not identify a gauge mode, physical propagation, a unique
evolution tick, `c`, `G`, a Planck scale, or particle masses. Even after an
adversarial confirmation of the nonhomogeneous negative, the remaining
homogeneous near-null direction must be resolved before interpreting the carrier
as dynamics.
