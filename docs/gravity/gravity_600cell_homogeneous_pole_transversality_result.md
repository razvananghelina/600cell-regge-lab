# Result: the full pole equation removes the weak homogeneous freedom

Date: 2026-08-20  
Status: **DERIVED COMPUTATIONAL from certified endpoint data**

## Mechanical outcome

The registered targeted verifier reports

```text
HOMOGENEOUS_WEAK_LINE_TRANSVERSE_TO_POLE_EQUATION
6/6 tests passed
```

Artifact SHA-256:

```text
d8fd2b0cd71d428d6cef5874b0cd6cf0496f174db13471bdb818a0803d182e0a
```

Verifier source SHA-256:

```text
1e75e13e5e07ef1f36f9273b30ab955f0340c2b5dd30ffd14ba44340400c3c91
```

Only this targeted verifier and the static registry audit ran.  The full suite
did not run.

## Certificate

At the accepted nonstatic endpoint, in both staircase parities,

```text
J = d(F_pole,F_momentum)/d(s,z),
det J = -2.473952264927475232698678784929269...e-6.
```

The calibrated determinant error bound is

```text
7.7528410931671789082714876386462...e-20,
```

so the determinant excludes zero by a factor `3.1910266638998e13`.  If

```text
v=(-F_momentum,z, F_momentum,s),
```

then algebraically

```text
dF_momentum(v)=0,
dF_pole(v)=-det J
             =2.473952264927475232698678784929269...e-6.
```

The projector of `v` agrees with the independently derived exact-action weak
line to `7.55e-42`; the two parity projectors agree exactly in the loaded
representation.

## Correct interpretation

- **DERIVED:** the unique weak homogeneous carrier/canonical line is not tangent
  to the complete fixed-input equations.  It is an off-shell lapse response.
- **DERIVED LOCAL:** the accepted nonstatic homothetic endpoint is isolated when
  incoming canonical data and dust mass are held fixed.
- **DERIVED NEGATIVE:** this line supplies neither a free tick nor a gauge
  direction.
- **PRESERVED:** the endpoint itself satisfies the full pole and seam equations;
  the result does not erase or refute that nonstatic first step.
- **OPEN:** existence, uniqueness and stability of the evolution map under
  iteration and perturbation.
- **NOT DERIVED:** an absolute tick unit, `c`, `G`, Planck units or particle
  masses.

The framing matters: a locally unique root has zero solution freedom at fixed
inputs precisely because the evolution output is selected.  The physical
candidate is the state-to-state map, not the kernel of its fixed-input endpoint
equations.

## Literature status

The pre-result primary-source gate confirms that this mechanism is known in
canonical simplicial gravity and consistent discretization: initially free
data can be fixed by later constraints, and curved Regge gauge constraints can
become pseudo-constraints.  This verifier is an internal classification, not an
external novelty claim.

