# Third canonical homothetic dust tick

Date: 2026-08-16

## Provenance

- prior-art gate: `7b9a676`;
- target-independent root enumeration: `3401137`;
- disclosed fixed-lapse comparison: `7cf4e27`;
- clean local-correction protocol: `1782b29`;
- registered implementation before any new corrected state: `f969712`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_dust_third_tick_local_correction.py`;
- artifact:
  `reproducible/gravity_600cell_dust_third_tick_local_correction.json`;
- artifact SHA-256:
  `ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0`.

Only this targeted verifier was run.  It returned **6/6**.  The full suite was
not run.

The earlier target-comparison diagnostic was disclosed as non-preregistered.
The present new-state Newton calculation was separately and cleanly frozen
before evaluation; it used only committed contracting root 0.

## Verdict

```text
THIRD_HOMOTHETIC_TICK_ACCEPTED
```

**DERIVED COMPUTATIONAL:** the fixed-mass contracting root admits a resolved
local correction satisfying the pole equation and the complete canonical
seam on both independently derived schedule parities.

The accepted absolute state is

```text
C3 = log(L3/L0)
   = -1.86963991043457648317744099662105411297713175985666389073915e-5,

R3 = log(rho3/rho0)
   = -3.20333682847923061806883161989269477486921081920153782268906e-5.
```

Relative to the accepted second output,

```text
u3 = log(L3/L2)
   = -9.34821204543993769543618697355300739949703565847672386319533e-6,

v3 = log(rho3/rho2)
   = -1.77963407327825031845582616746794662148542711258488403011932e-5,

L3/L2     = 0.9999906518316489581306198885769769461487,
rho3/rho2 = 0.9999822038176211498605569266287745415058,
tau3/tau2 = 0.9999911018692222092609450519281653044973.
```

The geometry therefore continues to contract while the proper tick duration
changes by only about 8.90 parts per million.

## Complete gates

For each parity:

```text
Newton iterations                  3,
all accepted dampings              1,
reduced residual infinity norm     3.120e-41,
canonical junction norm            1.709e-40,
inherited junction bound            3.651e-21,
maximum diagonal residual          2.45e-99 even / 1.60e-96 odd,
maximum pole residual              1.58e-45.
```

The endpoint Jacobian has

```text
singular values = (5.8283562156e2, 4.2445618126e-9),
epsilon         = 1.3300463377e-22.
```

Thus its weak direction is ill-conditioned but separated from the calibrated
error by about `3.19e13`.  Parity differences in `(C3,R3)` are below `1e-89`,
and all evaluated action, derivative and trial states retain the Lorentzian
branch.

## Preregistered sequence audit

The third-tick protocol recorded, without using them as acceptance targets,

```text
u3/A1 approximately 3,
v3/R1 approximately 5,
C3/A1 approximately 6,
R3/R1 approximately 9.
```

The result is

```text
u3/A1 = 3.0000106915017549,
v3/R1 = 5.0000210878312066,
C3/A1 = 6.0000133643718067,
R3/R1 = 9.0000253053809844.
```

Together, the three accepted ticks give

| tick `n` | increment `u_n/A1` | cumulative scale log `/A1` | increment `v_n/R1` | cumulative lapse log `/R1` |
|---:|---:|---:|---:|---:|
| 1 | `1` | `1` | `1` | `1` |
| 2 | `2.000002673` | `3.000002673` | `3.000004218` | `4.000004218` |
| 3 | `3.000010692` | `6.000013364` | `5.000021088` | `9.000025305` |

**PATTERN:** cumulative scale logs follow the triangular numbers
`n(n+1)/2`, while cumulative lapse logs follow `n^2`, to relative deviations
of a few parts in `1e-6`.  This is now a preregistered three-tick pattern, but
still not an exact identity or general recurrence theorem.  The deviations
are resolved and nonzero; the integers are asymptotic structure, not exact
equalities.

## Physical status

- **DERIVED COMPUTATIONAL:** three consecutive non-static homogeneous slabs
  glue canonically with one conserved dust mass.
- **DERIVED COMPUTATIONAL:** contraction accelerates in the discrete scale
  coordinate: the first three log-scale decrements are approximately in the
  ratio `1:2:3`.
- **STRUCTURAL:** this is a three-step isotropic Regge minisuperspace
  trajectory on one fixed carrier.
- **PATTERN:** the triangular/square law strongly suggests a low-order local
  recurrence or Taylor structure.
- **OPEN:** whether that structure follows analytically from the reduced
  action, survives tick four, or converges under refinement.
- **OPEN:** anisotropic stability and physical gravitational-wave modes.
- **OPEN:** absolute time.  `tau0=0.0102` is still an input, so no fundamental
  time unit, `c`, Planck time or Planck mass has been derived.

## Post-result primary-source audit

Successive 600-cell Regge evolution is already known broadly from the
implicit scheme of Barrett et al. and the dust evolution of De Felice--Fabri:
<https://arxiv.org/abs/gr-qc/9411008> and
<https://arxiv.org/abs/gr-qc/0009093>.  Lorentzian Regge shells with dust can
approximate continuum closed cosmology, and their accuracy depends on the
discretization: <https://arxiv.org/abs/2109.00875>.  Dust can serve as a
proper-time reference only in a fuller canonical dust construction such as
Brown--Kuchař: <https://arxiv.org/abs/gr-qc/9409001>.

No located primary source states the present triangular/square law for this
fixed staircase canonical map.  External novelty remains **OPEN**; search is
not proof.

## Next decisive step

Do not merely accumulate ticks.  First derive or refute the local recurrence
from the homothetic reduced action:

1. compute a target-independent high-precision Taylor jet of the canonical
   map around the static family;
2. determine whether triangular/square coefficients follow from the leading
   nonzero terms and the half-step placement of the initial static slab;
3. preregister the tick-four prediction from that derived jet;
4. only then calculate tick four as an out-of-sample test.

This separates an explainable Friedmann/constant-acceleration Taylor effect
from a genuinely new discrete recurrence.
