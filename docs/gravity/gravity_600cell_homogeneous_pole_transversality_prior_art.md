# Prior-art and framing gate: homogeneous pole transversality

Date: 2026-08-20  
Status: **prior-art/framing gate; no new calculation accepted here**

## Complete hypotheses

Work only in the frozen two-variable homothetic endpoint problem

```text
y=(s,z)=(log lambda, log(rho/rho0)).
```

At the already accepted nonstatic root define

```text
F0(y) = mean of the five complete pole/lapse equations,
F1(y) = mean of the old-momentum/canonical-seam mismatch.
```

The endpoint satisfies all 35 internal equations and all 30 seam components.
Its calibrated `2 x 2` Jacobian `J=d(F0,F1)/d(s,z)` is rank two in both derived
staircase parities.  The exact homogeneous weak-pole line is the nonzero
direction preserving `F1` to first order.

No claim is made outside the homogeneous sector, away from the accepted root,
or after changing the fixed incoming canonical datum, dust mass, carrier,
branch or derivative conventions.

## Primary literature and known mechanism

Action-generated canonical simplicial evolution and the later fixing of
initially free data are standard in Dittrich and Höhn,
<https://arxiv.org/abs/1108.1974>.  Curved Regge backgrounds generically replace
exact gauge constraints by background-dependent pseudo-constraints in Bahr and
Dittrich, <https://arxiv.org/abs/0905.1670>, and Dittrich and Höhn,
<https://arxiv.org/abs/0912.1817>.  Consistent discretization can determine
multiplier-like variables and leave a canonical transformation rather than a
continuum-like constraint family in Gambini and Pullin,
<https://arxiv.org/abs/gr-qc/0511096>.

Thus the general mechanism is **KNOWN**.  The calculation below is an internal
classification of the project's already derived line, not a novelty claim.

## Framing correction

Write

```text
J = [[F0_s,F0_z],
     [F1_s,F1_z]].
```

A generator of the weak line in `(delta s,delta z)` coordinates is

```text
v=(-F1_z,F1_s).
```

Then

```text
dF1(v)=0,
dF0(v)=-det(J).
```

Therefore calibrated full rank of `J` already predicts that the pole equation
is transverse to the weak line.  A positive transversality result must not be
described as “no evolution”.  It means:

- the weak line is not an additional solution freedom or free tick;
- at the fixed incoming canonical data, the accepted nonstatic endpoint is
  locally isolated;
- the nonstatic root itself remains a valid first canonical step.

The evolution object is the locally selected root/map, not a kernel obtained
while holding its inputs fixed.

## Status before verification

- **DERIVED upstream:** a nonstatic complete homothetic root exists.
- **DERIVED upstream:** the weak canonical graph has one exact homogeneous line.
- **PREDICTED from the two results:** the full pole equation eliminates that
  line as an infinitesimal freedom.
- **OPEN until the registered check:** convention alignment, determinant error
  certificate and agreement of the two independently stored line projectors.

