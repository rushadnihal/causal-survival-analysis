"""Save weighted survival / cumulative-incidence curves by arm to data/results/curves.csv (monthly grid, 0-120 months).

Run from the project root:  python3 scripts/export_curves.py     (needs data/processed/cohort_weighted.csv)
"""
import os
import numpy as np
import pandas as pd

from tte.survival import weighted_aj, at

d = pd.read_csv("data/processed/cohort_weighted.csv")
T, EV, A = d.PERMTH_EXM.values, d.event.values, d.acei_arb.values
schemes = [("IPTW", d.sw_trunc.values),
           ("IPTW x survey weight", (d.sw_trunc * d.svy_w).values),
           ("Unweighted", np.ones(len(d)))]

rows = []
for label, w in schemes:
    for a, arm in ((0, "Untreated"), (1, "ACEi/ARB")):
        m = A == a
        u, S, c1, c2 = weighted_aj(T[m], EV[m], w[m])
        for month in range(0, 121):
            rows.append({"weighting": label, "arm": arm, "month": month,
                         "survival": 1.0 if month < u[0] else at(u, S, month),
                         "cif_cvd_death": at(u, c1, month),
                         "cif_other_death": at(u, c2, month)})

os.makedirs("data/results", exist_ok=True)
pd.DataFrame(rows).to_csv("data/results/curves.csv", index=False)
print("Saved data/results/curves.csv with", len(rows), "rows")
