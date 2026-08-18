# Stationary momentum envelope: corrected base-bracket result

Date: 2026-08-16

## Provenance

- prior-art gate: `dedcbc6`;
- frozen protocol: `ed1cd6a`;
- initial implementation: `0bd3fc4`;
- initial reporting failure: `9322394`;
- reporting-only correction: `4f5ed12`;
- corrected artifact SHA-256:
  `c416d2f4956ac610fbcc38a8d499ccdc55864d0ed62ab4b6c44aa0794cd1067c`.

Only the targeted verifier was rerun.  It returned **5/5**.  The full suite
was not run.

## Mechanical verdict

```text
MOMENTUM_ENVELOPE_BASE_BRACKET_FAILED
```

The result is target independent: the verifier verified the accepted artifact
hash but did not parse any desired momentum.

## What the serialized attempts show

All 26 endpoints of the thirteen frozen symmetric brackets remain on the same
Lorentzian branch (`2400` simplices with one negative direction).  The failure
is not a branch-domain failure.

At the first two widths:

```text
b                         G(b,r1)
-9.3481787301e-6          -8.6128916602e-13
-3.1160595767e-6          -3.3956377441e-8

-12.4642383068e-6         +6.7909809150e-8
 0                        +4.2574745356e-31
```

Thus the symmetric endpoint product is positive at each width even though the
left endpoints at consecutive widths have opposite signs.  The root lies in
the annulus between `-12.4642e-6` and `-9.34818e-6`, not across either frozen
symmetric pair.

Moreover, `b=0` satisfies `|G|<1e-25`.  This is the structurally expected
time-reversed first slab: lower scale `L1`, upper scale `L0`, and the same
lapse.  It is a second stationary branch/root.

## Verdict on the framing

- **DERIVED COMPUTATIONAL:** the complete target-free bracket attempts and
  their Lorentzian branch data are now recorded.
- **DERIVED NEGATIVE:** the preregistered symmetric-bracket algorithm fails.
- **REFUTED:** the implicit assumption that one symmetric expanding bracket
  can isolate the base stationary root.  `G(b,r1)` is not monotone over those
  endpoints and contains multiple roots.
- **STRUCTURAL:** `b=0` is the time-reversal root.  It is not automatically a
  canonical forward continuation because time reversal changes the boundary
  momentum sign.
- **OPEN:** the exact count of stationary roots on a fixed `b` domain, which
  branch is forward-connected, and each branch's momentum envelope.

No momentum target has yet been compared.

## Required next correction

Do not choose the left root because it looks like the desired continuation.
Preregister a target-independent **root enumeration** over a fixed `b` domain
at `r=r1`, record every sign bracket and exact/near-zero node, and commit the
full root list before tracing any branch or comparing momentum targets.

Each enumerated root must then seed its own stationary curve.  This prevents
silently selecting the branch that best fits the desired canonical momentum.
