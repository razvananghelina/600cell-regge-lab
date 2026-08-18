# Homothetic forward-root test: result

Date: 2026-08-16

## 1. Provenance

- prior-art gate: `2ce58bd`;
- frozen protocol: `29653b9`;
- registered implementation before evaluation: `bbfcf04`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_homothetic_forward_root.py`;
- artifact:
  `reproducible/gravity_600cell_dust_homothetic_forward_root.json`;
- artifact SHA-256:
  `7a83d800156444bff2b549ad3eb34bd0ba506c96353fa452e33360c0d38d879a`.

Only the targeted verifier was run.  It returns **6/6**.  The full suite was
not run.

## 2. Mechanical verdict

Both parities return

```text
HOMOTHETIC_STATIONARY_NOT_CANONICAL.
```

This is an honest negative at the fixed published lapse `tau=0.0102`.  The
predicted non-static slab exists and satisfies every internal equation, but
it does not glue canonically to the preceding static sandwich within the
committed uncertainty.

## 3. The non-static stationary slab is real

The frozen bracket passes independently for both schedules:

```text
E(-2.5e-5)    = +4.78286735158e-6,
E(-1.5625e-6) = -2.12226529489e-8.
```

Eighty bisections give the identical root

```text
s_root = log(L_+/L_-)
       = -3.1160706675973036032799885512889e-6,

L_+/L_- = 0.99999688393418734585634357203427,
final bracket width = 1.9387045607e-29.
```

Thus the upper slice contracts by approximately 3.116 parts per million.
At the root,

```text
abs(E_lapse)             = 2.159e-31,
max abs(all 35 equations)= 4.318e-32,
max abs(diagonal eqs)    < 1.2e-96.
```

Every evaluated simplex remains Lorentzian and away from the angle branch
cut.  Even and odd roots, scale ratios and momenta agree far inside all
frozen parity gates.

**DERIVED COMPUTATIONAL.**  At the chosen lapse, fixed conserved mass and
the exact homothetic geometry select a non-static action-stationary
contracting slab.  This upgrades the earlier Taylor pattern to a computed
root and confirms that the 30 diagonal equations vanish on this state, not
only their global sum.

## 4. Continuum control

The time-symmetric closed-dust continuum estimate frozen before the root was

```text
s_cont = -2.88796961781375994885e-6.
```

The finite 600-cell result gives

```text
s_root/s_cont = 1.07898318887309...
```

so its contraction is about 7.9% larger in magnitude.  This is a
**STRUCTURAL CONTROL**, not a fit: the root was obtained solely from the
Regge equation.  The difference is compatible in scale with a coarse finite
discretization but no convergence statement follows from one resolution.

## 5. Canonical junction fails decisively

The root's 30 pre-momenta are uniform:

```text
p_pre(root) = 0.0009081060650005251076... per component.
```

The independently committed forward target is

```text
P p_post(static) = 0.0009081044489065315762... per component.
```

Hence every component misses by

```text
1.61609399353141e-9,
```

and

```text
||p_pre(root)-P p_post(static)||_2 = 8.85171135306e-9,
allowed bound                           = 3.65136539620e-21.
```

The mismatch is twelve orders of magnitude beyond the frozen bound.  Its
within-component spread is only `1.14e-39`, so this is not an anisotropic or
orbit-matching accident.  It is a homogeneous scalar mismatch.

**DERIVED NEGATIVE at fixed lapse.**  The stationary root is not the
canonical next frame of the published sandwich.  Calling it the first tick
would be wrong even though it is geometrically plausible and extremely
close to the desired momentum.

## 6. What the near miss means

The failure is small physically but not numerically: the per-component
relative mismatch is about `1.78e-6`.  It cannot be excused by tolerance.

It also sharpens Claude's proposal.  This mission fixed `tau`.  There are now
two scalar equations for a genuinely forward homogeneous slab:

```text
E_internal(s,rho) = 0,
p_pre(s,rho)-P p_post(static) = 0.
```

Allowing the next slab's positive `rho=tau^2` to vary may select a pair
`(s,rho)` and remove the mismatch.  If so, the selection would come from
internal stationarity plus canonical consistency, not from mass conservation
alone.  In a finite curved discretization it must initially be classified as
a possible pseudo-constraint/consistent-discretization effect, not
immediately as an emergent physical clock.

## 7. Post-result primary-source audit

The general possibility is known.  Gambini and Pullin's consistent
discretization of Regge gravity replaces continuum constraints by a
well-defined discrete canonical transformation:
<https://arxiv.org/abs/gr-qc/0511096>.  Bahr and Dittrich explain how broken
discrete diffeomorphism symmetry yields background-dependent
pseudo-constraints: <https://arxiv.org/abs/0905.1670>.  Canonical
pre/post-momentum matching is standard in Dittrich and Hoehn:
<https://arxiv.org/abs/1108.1974>.

Therefore a discrete equation determining a lapse-like variable would not
be conceptually new.  No located source gives the present two-variable
root on this exact carrier.  External novelty remains **OPEN**.

## 8. Next falsification test

The next mission should vary only the next slab's lapse and scale while
keeping all prior canonical data fixed.  Before evaluation it must freeze:

1. the two equations above and their coordinates `(s, log(rho/rho0))`;
2. a deterministic local Jacobian calibration at the stationary root;
3. one Newton seed fixed to `(s_root,0)`;
4. full 35-equation, 30-momentum, parity and Lorentzian gates;
5. outcomes separating a selected positive lapse, rank loss, branch loss,
   multiple-root ambiguity and numerical openness.

A passing pair would establish the first canonically glued homogeneous slab
and a discretely selected next lapse on this carrier.  It would still not
establish refinement stability, a universal clock or multi-tick evolution.
