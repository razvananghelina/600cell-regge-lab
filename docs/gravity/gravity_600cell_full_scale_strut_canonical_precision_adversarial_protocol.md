# Preregistered adversarial replication of nonhomogeneous full rank

Date: 2026-08-20  
Status: **preregistered adversarial audit**

## Scope

The primary result frozen in commit `99efe9b` reports full column rank for every
nonhomogeneous D and K block and leaves the homogeneous block **OPEN**. This audit
tests only the material nonhomogeneous rank claim. It may not resolve or
reinterpret the homogeneous candidate.

## Independence of the decisive step

The primary verifier certifies full column rank by an interval determinant of the
full Gram matrix `M* M`. The adversarial verifier will not use a Gram determinant
or a singular-value threshold as its decisive certificate.

For each matrix M it will instead:

1. select a square row minor deterministically by column-pivoted QR of the
   floating midpoint transpose `M.T`;
2. wrap every selected multiprecision entry in a conservative complex ball whose
   real and imaginary component radii both equal the primary matrix's full
   Frobenius ball radius;
3. compute the determinant of that square minor directly in Arb/Acb;
4. certify full column rank only if that direct determinant ball excludes zero.

The full Frobenius radius bounds every component error; using it independently on
both real and imaginary components enlarges rather than understates the primary
uncertainty.

## Cross-precision anti-selection test

For each parity, each of the six nonhomogeneous sectors, and each of D and K:

- rows selected from P100 must certify the corresponding P160 direct minor;
- rows selected from P160 must certify the corresponding P100 direct minor.

Thus the row choice at a precision level cannot certify only the matrix on which
it was selected. The expected census is 48 certified direct determinant balls:
`2 parities x 6 sectors x 2 matrices x 2 cross-directions`.

## Controls and frozen conventions

- Known pass: an exact padded identity matrix must yield a direct minor excluding
  zero.
- Known fail: appending an exact duplicate column must yield a selected direct
  determinant containing zero.
- Every primary input, source and artifact is pinned by SHA-256; executing the
  primary source under capture must reproduce `17/17`, the frozen outcome, and
  byte-identical artifact content.
- The audit uses the already frozen source/target convention. The primary
  source/target reversal and pole-deletion falsifiers must still be recorded as
  passed in the frozen artifact.

No scientific target, physical interpretation, fit, or new tolerance is used to
select a row minor.

## Outcomes

- All 48 cross-precision determinants exclude zero and both controls behave as
  planted: `NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED`.
- Any nonhomogeneous determinant contains zero, a control misbehaves, or frozen
  provenance fails: `NONHOMOGENEOUS_DIRECT_MINOR_DISAGREEMENT`.

Only the first outcome permits consolidation of the nonhomogeneous zero
intersection as **DERIVED**. The homogeneous intersection remains **OPEN** in
either case.
