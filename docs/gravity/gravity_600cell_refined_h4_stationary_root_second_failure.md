# Second execution failure: an invalid Lorentzian reality criterion

Date: 2026-08-20

First correction commit: `69ca0f3`  
Corrected-anchor implementation: `7767343`

The second targeted execution again stopped at `10/11`, before all 120 root
attempts. Its artifact has SHA-256

```text
a57367a2766cd62e0b24fb84471f025cdb8e642c8606c608377673e9f75ebad9
```

The smaller preregistered anchors had nonzero angle arguments and angle
identities accurate to about `1e-78`. Direct diagnosis on the first anchor
gave

```text
|Im complete action|       = 7.39e-77
max |Im action gradient|   = 3.16e-77
max |Im hinge curvature|   = 3.24e-4
minimum angle argument     = 0.92796...
```

The implementation had included the final line in a generic
`maximum_imaginary` branch gate. That criterion is physically wrong for the
complex Lorentzian Regge convention used here. Individual spacelike
hinge/boundary curvatures may be imaginary while their areas and the overall
`-i` continuation give a real Lorentzian action and real equations. The
repository's independently accepted tent calculation already documents this
structure.

Thus the second failure is not evidence that the anchors left the Lorentzian
branch. It exposes a category error in the control itself. Shrinking anchors
again would hide rather than repair it.

The next correction must gate the minimum logarithm-argument modulus, angle
identity, complete-action reality and complete-gradient reality. It must
record the hinge-curvature imaginary component diagnostically but must not
require it to vanish. The same correction applies to trial points and final
root validation. No search bounds, seeds, solver settings, attempts or
look-elsewhere denominator may change.
