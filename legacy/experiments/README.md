# Historical exploratory scripts

This directory contains the old `exp*.py` research scripts formerly stored
at repository root.  They range from useful construction probes to obsolete
or refuted speculation.  Their presence is archival and does not certify
their claims.

The authoritative executable evidence is the registered collection in
[`../../reproducible/`](../../reproducible/), together with the corresponding
protocol and result notes.

Because many historical scripts import modules from repository root, invoke
one from the root with the project on `PYTHONPATH`, for example:

```bash
PYTHONPATH=. /home/razvan/science/.venv/bin/python \
  legacy/experiments/exp015_600cell_spectrum.py
```

Two scripts intentionally remain at root because a registered provenance
verifier reads their source text at the original path:

```text
exp261_spectral_action.py
exp513_cc_spectral_action.py
```
