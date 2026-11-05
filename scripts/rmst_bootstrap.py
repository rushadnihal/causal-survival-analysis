"""Restricted mean times (months) by arm and their IPTW differences, with bootstrap CIs.

For a horizon tau: months alive = int S, months lost to CVD death = int CIF1, months lost to other death = int CIF2.
They sum to tau for each arm, so the three differences sum to zero.

Run from the project root:  python3 scripts/rmst_bootstrap.py [B]    (default B = 500)
"""
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from tte.psmodel import build_X
from tte.survival import weighted_aj, restricted_mean
from tte.weighting import stabilized_weights

B = int(sys.argv[1]) if len(sys.argv) > 1 else 500
d = pd.read_csv("data/processed/cohort_weighted.csv")
X = build_X(d).values.astype(float)
A, T, EV = d.acei_arb.values, d.PERMTH_EXM.values, d.event.values


def estimate(idx):
    a, t, ev = A[idx], T[idx], EV[idx]
    Xs = StandardScaler().fit_transform(X[idx])
    ps = LogisticRegression(max_iter=2000).fit(Xs, a).predict_proba(Xs)[:, 1]
    w = stabilized_weights(ps, a, trunc=(1, 99))
    out = []
    for tau in (60, 120):
        r = {}
        for g in (0, 1):
            m = a == g
            u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
            r[g] = (restricted_mean(u, S, tau, 1.0),
                    restricted_mean(u, c1, tau, 0.0),
                    restricted_mean(u, c2, tau, 0.0))
        out += [r[1][k] - r[0][k] for k in range(3)]
    return out


n = len(d)
point = estimate(np.arange(n))
rng = np.random.default_rng(3)
boot = np.array([estimate(rng.integers(0, n, n)) for _ in range(B)])
lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
names = [f"{k} @{tau}mo" for tau in (60, 120) for k in ("Months alive", "Months lost to CVD death", "Months lost to other death")]
print(f"Difference in restricted mean time (ACEi/ARB minus untreated), months; IPTW; bootstrap B={B}")
for i, nm in enumerate(names):
    print(f"{nm:34s} {point[i]:+.2f}  (95% CI {lo[i]:+.2f} to {hi[i]:+.2f})")
print("\nCheck: the three differences at each horizon should sum to ~0:",
      round(sum(point[0:3]), 8), round(sum(point[3:6]), 8))

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_pd.DataFrame({"quantity": names, "estimate": list(point), "ci_lo": list(lo), "ci_hi": list(hi)}).to_csv("data/results/rmst.csv", index=False)
