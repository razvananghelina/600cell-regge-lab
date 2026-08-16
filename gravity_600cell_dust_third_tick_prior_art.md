# Prior-art gate: third canonical homothetic dust tick

Date: 2026-08-16

## Exact object and hypotheses

Use the fixed 600-cell order-24 staircase carrier, both derived schedule
parities, the complete 100-decimal Lorentzian Regge plus five-orbit dust
action and the same canonical boundary permutation as the accepted first two
ticks.  The dust mass remains

```text
M=(90/pi)*(2*pi-5*acos(1/3))*L0
```

and is not recomputed from a later boundary scale.

The accepted second output is the next lower boundary:

```text
B2 = log(L2/L0)
   = -9.34818705890582713633822299265753373027428194008991504419612e-6,

R2 = log(rho2/rho0)
   = -1.42370275520098029961300545242474815338378370661665379256974e-5.
```

The first mission is target-independent: at fixed `(B2,R2)`, enumerate the
stationary roots in the absolute upper log `C=log(L3/L0)` without parsing or
comparing the accepted second tick's post-momentum.  Only after that multiset
is committed may a later mission load the canonical target.

## Primary prior art

- Barrett, Galassi, Miller, Sorkin, Tuckey and Williams describe an implicit
  Regge evolution scheme and illustrate it on a 600-cell Friedmann cosmology.
  Successive 600-cell evolution is therefore **KNOWN** in broad form:
  <https://arxiv.org/abs/gr-qc/9411008>.
- De Felice and Fabri evolve a dust-filled 600-cell and explain the finite
  causal endpoint encountered by the evolution:
  <https://arxiv.org/abs/gr-qc/0009093>.
- Dittrich and Höhn derive canonical discrete evolution from Hamilton's
  principal function, including pre/post momenta and later consistency
  conditions:
  <https://arxiv.org/abs/1108.1974>.
- Brown and Kuchař show how physical dust can supply proper-time reference
  variables in canonical continuum gravity:
  <https://arxiv.org/abs/gr-qc/9409001>.

The last result does not by itself turn the present fixed world-line mass term
into a Brown--Kuchař clock.  Our starting `tau0` remains an external datum.

## KNOWN / CONTROL / OPEN

- **KNOWN:** multiple-step 600-cell Regge evolution and action-generated
  canonical evolution exist in the literature.
- **CONTROL:** the first and second canonical ticks are committed and pass all
  internal, momentum, branch and parity gates.
- **CONTROL:** the second-tick artifact is hashed but must not be parsed by the
  root-enumeration verifier; `B2` and `R2` are frozen above as literal
  geometry values.
- **DERIVED:** mass is conserved across the first two ticks and unequal scales
  break the exact static all-lapse cancellation.
- **OPEN:** the stationary-root multiset at `(B2,R2)`.
- **OPEN:** whether a geometrically forward root exists and whether it can be
  corrected to the third canonical target.
- **OPEN:** the apparent triangular/square sequence, global recurrence,
  continuum convergence, stability and absolute clock selection.

## Framing attack

A successful third tick would establish a three-step homogeneous trajectory,
not a general evolution theorem.  Agreement with the preregistered integer
ratios would remain a **PATTERN** because three small weak-field increments can
share low-order Taylor structure.

The finite grid can exhaust sign-changing and grid-node roots only on its
frozen domain.  It cannot falsify even-multiplicity tangential roots between
nodes or roots outside that domain.  Such omissions must remain **OPEN** and
cannot be hidden behind the word "complete".

A literature search cannot prove external novelty.  Only targeted verifiers
will be run; the full suite will not be run.
