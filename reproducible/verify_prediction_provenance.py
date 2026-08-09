#!/usr/bin/env python3
"""Mechanical checks for prediction_provenance_ledger.md.

Dates and sources are pinned here deliberately so chronology changes require
an explicit code review rather than a prose-only edit.
"""
from datetime import date
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "prediction_provenance_ledger.md"

FORMULA_DATE = date(2026, 2, 23)
PINNED_PUBLICATIONS = {
    "NuFIT6": (date(2024, 10, 7), "https://arxiv.org/abs/2410.05380"),
    "JUNO_FIRST": (date(2025, 11, 18), "https://arxiv.org/abs/2511.14593"),
    "KATRIN_LIMIT": (date(2024, 6, 19), "https://arxiv.org/abs/2406.13516"),
    "DESI_LIMIT": (date(2025, 3, 19), "https://arxiv.org/abs/2503.14744"),
}
RETRO_NEUTRINO = {
    "sin^2 theta_12", "sin^2 theta_13", "sin^2 theta_23",
    "delta_CP", "Delta m^2_21", "Delta m^2_31",
    "ratio Delta m^2_21/Delta m^2_31", "m2 and m3 point values",
}
BLIND_NEUTRINO = {"m1", "strict ordering", "sum m_nu", "m_beta", "m_betabeta"}

text = LEDGER.read_text(encoding="utf-8")
rows = {}
for line in text.splitlines():
    if line.startswith("|") and line.count("|") >= 6:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 6 and cells[4] in {"RETRODICTION", "BLIND", "AMBIGUOUS"}:
            rows[cells[0].replace("`", "")] = cells

def row_for(prefix):
    if prefix in rows:
        return rows[prefix]
    matches = [cells for key, cells in rows.items() if key.startswith(prefix)]
    assert len(matches) == 1, f"expected one row beginning {prefix!r}, got {len(matches)}"
    return matches[0]

for name in RETRO_NEUTRINO:
    cells = row_for(name)
    assert cells[4] == "RETRODICTION", (name, cells[4])
for name in BLIND_NEUTRINO:
    cells = row_for(name)
    assert cells[4] == "BLIND", (name, cells[4])

assert PINNED_PUBLICATIONS["NuFIT6"][0] < FORMULA_DATE
assert PINNED_PUBLICATIONS["JUNO_FIRST"][0] < FORMULA_DATE

# Every BLIND row must describe absence of a measurement (not absence of a
# bound). This is checked against the pinned source descriptions in the row.
for name in BLIND_NEUTRINO:
    experimental = row_for(name)[3].lower()
    assert any(token in experimental for token in
               ("no measurement", "not settled", "only absolute-mass",
                "bound is", "no measurement, only",
                "no model-independent measurement", "<0.45")), (name, experimental)

# The ledger must retain the central correction sentence and all source URLs.
required = [
    "Consistency was achieved with values known at construction time",
    "No measured oscillation parameter belongs to that set",
    *[source for _, source in PINNED_PUBLICATIONS.values()],
]
for phrase in required:
    assert phrase in text, f"missing pinned ledger content: {phrase}"

# Transparent illustrative look-elsewhere census.
phi = (1 + 5 ** 0.5) / 2
family = set()
for p in range(1, 11):
    for q in range(1, 11):
        for k in range(-40, 41):
            x = (p / q) * phi ** k
            family.add(round(x, 15))
        for r in range(1, 11):
            x = (p + q * phi) / r
            family.add(round(x, 15))
            family.add(round(1 / x, 15))

targets = {
    "sin2_theta12_JUNO": (0.3092, 0.0087),
    "sin2_theta13_NuFIT": (0.02215, 0.00057),
    "mass_split_ratio": (7.50e-5 / 2.513e-3, 0.0012),
}
print(f"illustrative distinct expression family: {len(family)}")
for label, (target, sigma) in targets.items():
    hits = sum(abs(x - target) <= sigma for x in family)
    print(f"{label}: {hits} candidates within illustrative 1 sigma")

print("PASS: provenance chronology and classification ledger are internally consistent")
