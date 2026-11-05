"""Sensitivity of the IPTW cumulative-incidence risk differences to analytic choices.

Run from the project root:  python3 scripts/sensitivity.py [B]     (B = bootstrap replicates, default 200)
"""
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from tte.psmodel import build_X
from tte.survival import weighted_aj, at

B = int(sys.argv[1]) if len(sys.argv) > 1 else 200
d = pd.read_csv("data/processed/cohort_weighted.csv")
X = build_X(d).values.astype(float)
A, T, EV = d.acei_arb.values, d.PERMTH_EXM.values, d.event.values


def make_weights(scheme, ps, a):
    p = a.mean()
    sw = np.where(a == 1, p / ps, (1 - p) / (1 - ps))
    if scheme == "none":
        return sw
    if scheme == "trunc99":
        lo, hi = np.percentile(sw, [1, 99]); return np.clip(sw, lo, hi)
    if scheme == "trunc95":
        lo, hi = np.percentile(sw, [5, 95]); return np.clip(sw, lo, hi)
    if scheme == "overlap":
        return np.where(a == 1, 1 - ps, ps)
    if scheme == "trim":
        return np.where((ps >= 0.05) & (ps <= 0.95), sw, 0.0)
    raise ValueError(scheme)


def estimate(idx, scheme):
    a, t, ev = A[idx], T[idx], EV[idx]
    Xs = StandardScaler().fit_transform(X[idx])
    ps = LogisticRegression(max_iter=2000).fit(Xs, a).predict_proba(Xs)[:, 1]
    w = make_weights(scheme, ps, a)
    out = []
    for h in (60, 120):
        r = {}
        for g in (0, 1):
            m = (a == g) & (w > 0)
            u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
            r[g] = (at(u, c1, h), at(u, c2, h))
        out += [100 * (r[1][0] - r[0][0]), 100 * (r[1][1] - r[0][1])]
    return out          # cvd60, other60, cvd120, other120


all_idx = np.arange(len(d))
not_0506 = np.where(d.cycle.values != "2005_2006")[0]
variants = [("Main (truncate 1/99)", "trunc99", all_idx),
            ("No truncation", "none", all_idx),
            ("Truncate 5/95", "trunc95", all_idx),
            ("Overlap weights*", "overlap", all_idx),
            ("Trim PS outside 0.05-0.95", "trim", all_idx),
            ("Drop 2005-06 cycle", "trunc99", not_0506)]

rng = np.random.default_rng(5)
SENS = []
labels = ["CVD @60", "Other @60", "CVD @120", "Other @120"]
print(f"Risk difference (ACEi/ARB minus untreated), percentage points, IPTW only, bootstrap B={B}\n")
print(f"{'variant':28s}" + "".join(f"{l:>26s}" for l in labels))
for name, scheme, idx in variants:
    pt = estimate(idx, scheme)
    boot = np.array([estimate(rng.choice(idx, len(idx)), scheme) for _ in range(B)])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    SENS.append([name] + [x for i in range(4) for x in (pt[i], lo[i], hi[i])])
    print(f"{name:28s}" + "".join(f"{pt[i]:+7.2f} ({lo[i]:+.1f},{hi[i]:+.1f})".rjust(26) for i in range(4)))
print("\n* overlap weights target a different population (people with similar chances of either treatment)")

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_cols = [f"{l}_{s}" for l in ["cvd60", "other60", "cvd120", "other120"] for s in ("estimate", "ci_lo", "ci_hi")]
_pd.DataFrame(SENS, columns=["variant"] + _cols).to_csv("data/results/sensitivity.csv", index=False)
