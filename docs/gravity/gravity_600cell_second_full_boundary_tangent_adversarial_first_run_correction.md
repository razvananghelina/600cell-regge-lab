# Correction after the first dense second-tangent adversarial run

Date: 2026-08-22

Status: **FROZEN BEFORE MODIFYING OR RE-EXECUTING THE VERIFIER.**

Failure-preservation commit: `94d6b3b`  
Failure note SHA-256:
`56bf30a4046c8b66162c2cd7aaf77acecc92edd6d88c655123b81b9a10ddeab5`  
Failed artifact SHA-256:
`2338fd58ab50ff309b50bc7a9beeac5596bc362cd6d07732e23c3a125acc9`

Only two mechanical corrections are licensed. The dense Hessian assemblies,
rank classifiers, solves, canonicality tests, scale tests, schedule tests and
the frozen 10/100 classifier are unchanged.

## 1. Delayed basis precision

The primary archive was constructed at 180 decimal digits. The adversarial
dense decision correctly used its preregistered 120 digits, but its delayed
closure incorrectly reconstructed the primary coordinate basis at 120 digits
and compared entries with a 180-digit-basis archive.

A target-free diagnostic on the 24-element regular carrier compared the same
deterministic eigenspace construction at 120 and 180 digits. For sector
dimensions `3,2,2,2,1,1,1`, the Frobenius distances between bases were

```text
2.5227793, 2.3102993, 1.9552270, 1.4117469,
4.8359617e-17, 4.8359617e-17, 2.1159317e-18.
```

For the four non-scalar sectors, the 120-to-180 overlaps were unitary to
`9.8e-17` but differed from the identity by `1.41` to `2.52`. This explains
the observed split without using a second-tangent target: every scalar-sector
closure agreed and every non-scalar closure was expressed in a different
unitary frame.

The correction is to run only the post-firewall deterministic basis
construction inside `mp.workdps(180)`. This cannot alter a dense scientific
label. Entrywise comparison is meaningful only after matching the disclosed
primary coordinate convention. No conjugation variant is searched or fitted.

## 2. Binary error propagation for a product

The first implementation used `kappa_1*kappa_2` as the binary amplification
for an already computed product `P=T_2 T_1`. That is the condition estimate
for a different composed solve, not the first-order multiplication error.
For independently solved maps,

```text
delta P = delta T_2 T_1 + T_2 delta T_1 + delta T_2 delta T_1.
```

Hence the first-order relative binary term is proportional to
`u*(kappa_1+kappa_2)`, with the quadratic `u^2*kappa_1*kappa_2` term added
explicitly. The observed middle-level estimates were approximately

```text
kappa_1 = 114.28, kappa_2 = 117.26,
kappa_1+kappa_2 = 231.54, kappa_1*kappa_2 = 13401.01.
```

The corrected product conditioning term is therefore

```text
kappa_product = kappa_1 + kappa_2
                + machine_epsilon*kappa_1*kappa_2.
```

This formula is frozen before re-evaluating the `1e-3` hostile control. The
synthetic perturbation amplitude, schedule comparisons and 10/100 thresholds
remain unchanged.

If either correction does not resolve its corresponding control, the next
execution remains a failure and must be preserved separately. No full suite
is licensed.
