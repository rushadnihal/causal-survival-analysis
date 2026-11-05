"""Weighted cumulative incidence of CVD / other death by arm, with validation checks.

Run from the project root:  python3 scripts/survival_report.py
"""
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter

from tte.survival import weighted_aj, at

d = pd.read_csv("data/processed/cohort_weighted.csv")
t, ev, A = d.PERMTH_EXM.values, d.event.values, d.acei_arb.values

# validation 1: no competing events -> CIF equals 1 - Kaplan-Meier
ev_nc = np.where(ev == 2, 0, ev)
u, S, c1, c2 = weighted_aj(t, ev_nc, np.ones(len(t)))
km = KaplanMeierFitter().fit(t, ev_nc == 1)
gap = max(abs(at(u, c1, x) - (1 - km.predict(x))) for x in (24, 60, 120))
print("Validation 1 (CIF vs 1-KM, max gap):", round(gap, 10))

# validation 2: S + CIF1 + CIF2 == 1
u, S, c1, c2 = weighted_aj(t, ev, d.sw_trunc.values)
print("Validation 2 (max |S+CIF1+CIF2-1|):", round(np.abs(S + c1 + c2 - 1).max(), 10))

for label, w in [("IPTW only", d.sw_trunc.values),
                 ("IPTW x survey weight", (d.sw_trunc * d.svy_w).values)]:
    print(f"\n=== {label} ===")
    res = {}
    for a in (0, 1):
        m = A == a
        u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
        res[a] = {h: (at(u, c1, h), at(u, c2, h)) for h in (60, 120)}
    for h in (60, 120):
        u0, u1 = res[0][h], res[1][h]
        print(f"{h} months | CVD-death risk: untreated {u0[0]*100:.1f}% , ACEi/ARB {u1[0]*100:.1f}% , "
              f"difference {100*(u1[0]-u0[0]):+.1f} pts | other-death risk: {u0[1]*100:.1f}% vs {u1[1]*100:.1f}%")
