# Protocol: does the spectral action preserve the round--Regge selector?

Date: 2026-08-12

## Provenance

This is a post-recognition protocol.  Before it was written, the continuous
certificate at commit `c97349c` had already proved that the normalized
ordinary complete-de-Rham coefficient obeys

```text
A2_eq'(u)<0,  0<=u<1,
```

on the frozen round--Regge path.  It was also already known from the standard
spectral-action literature that the cutoff function is normally required to
be positive.  The purpose of this audit is not blind discovery.  It is to
prevent three distinct statements from being conflated after the favorable
path result:

1. the sign of the asymptotic `A2` weight;
2. the magnitude of that weight and the absolute cutoff scale;
3. minimization of the complete finite-cutoff action.

No measured gravitational constant, desired endpoint, fitted cutoff, or
chosen heat time may enter the audit.

## Complete hypotheses

1. `D_u` is the self-adjoint complete-exterior de Rham operator with the
   transmittal/conic domain already used in the path certificate, and
   `P_u=D_u^2`.

2. The heat convention is exactly

   ```text
   Tr exp(-t P_u)
     ~ (4*pi*t)^(-3/2) [A0(u)+t*A2(u)+higher terms].
   ```

3. Equal-volume normalization is imposed before comparison, so `A0(u)` is
   constant on the path.

4. The bosonic spectral action has the standard form

   ```text
   S_chi(Lambda,u)=Tr chi(D_u/Lambda),
   ```

   where `chi` is even, nonnegative, not identically zero, and has the decay
   and regularity needed for the asymptotic expansion.  Put
   `F(x)=chi(sqrt(x))`, so `S=Tr F(P/Lambda^2)`.

5. Positivity of `chi` is an explicit spectral-action axiom.  It is not to be
   relabelled as a result derived from the 600-cell.

6. A statement about the complete finite-cutoff action requires either its
   actual spectrum on the whole path or a uniform remainder bound including
   every earlier/same-order singular term.  The `A2` certificate alone may
   not be cited for that stronger statement.

## Frozen derivation

From Mellin/Laplace functional calculus, the first two terms in dimension
three must be written with all constants exposed:

```text
S_chi(Lambda,u)
 ~ (4*pi)^(-3/2) [
      Lambda^3 C0(chi) A0(u)
    + Lambda   C2(chi) A2(u)
    + lower powers / singular remainders],

C0(chi)=4/sqrt(pi) integral_0^infinity chi(v) v^2 dv,
C2(chi)=2/sqrt(pi) integral_0^infinity chi(v) dv.
```

The verifier must establish `C0>0` and `C2>0` under hypothesis 4.  Thus a
positive cutoff cannot reverse the sign of the `A2` contribution.

For the explicit positive family

```text
chi_a(v)=exp(-a v^2),  a>0,
```

derive exactly

```text
C0=a^(-3/2),  C2=a^(-1/2),
C2/C0=a.
```

This family is the frozen counterexample to any claim that positivity alone
fixes the relative weight or heat scale.

Reuse the already preregistered exact spectra

```text
X={0,1,10},  Y={0,2,3}
```

only as a general finite-cutoff control: their positive heat-trace ordering
reverses between `exp(-t)=100/101` and `1/2`.  It must not be represented as
a spectrum of the round--Regge path.

## Decision boundary

- **DERIVED CONDITIONAL SIGN:** `C2(chi)>0` for every standard positive
  cutoff, hence the asymptotic `A2` term preserves the certified round
  preference.  At each fixed `u<1`, the round endpoint wins for sufficiently
  large cutoff if the stated asymptotic remainder is valid there.

- **REFUTED SIGN:** an admissible nonnegative cutoff gives `C2<=0` or the
  coefficient convention reverses the path ordering.

- **DERIVED COMPLETE ACTION SELECTION:** in addition to the sign result, an
  already selected `chi`, `Lambda`, and controlled full remainder give the
  same unique round minimum uniformly for every `u`.

- **PARTIAL / OPEN PHYSICAL GATE:** the sign is fixed but the cutoff moment,
  scale, or higher terms remain free.  This outcome closes the sign objection
  only; it does not derive a gravitational action or Newton's constant.

## Mandatory framing

Even the strongest sign outcome is universal spectral-action structure, not
a special prediction of `a1=5` or the 600-cell.  What is specific to the
construction is the certified `A2(u)` geometry.  The note must report both
facts separately.

No full-suite run is authorized.  Register and run only the targeted
verifier.
