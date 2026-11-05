import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from tte.psmodel import build_X
from tte.survival import weighted_aj, at
import sys
MAXD = float(sys.argv[1]) if len(sys.argv) > 1 else np.inf
d = pd.read_csv("data/processed/cohort_weighted.csv")
d = d[(d.acei_arb == 0) | (d.acei_arb_days <= MAXD)].reset_index(drop=True)
print("N =", len(d), "| treated =", int(d.acei_arb.sum()), "| max treated duration (days):", MAXD)
X = build_X(d).values.astype(float)
A, T = d.acei_arb.values, d.PERMTH_EXM.values
nc = (d.MORTSTAT == 1) & d.UCOD_LEADING.isin([2, 3, 4])
EV = np.where(nc, 1, np.where(d.MORTSTAT == 1, 2, 0))
print('Negative-control deaths (cancer/lung/accident):', int(nc.sum()), '| all other deaths:', int((EV == 2).sum()))
def estimate(idx):
    a, t, ev = A[idx], T[idx], EV[idx]
    Xs = StandardScaler().fit_transform(X[idx])
    ps = LogisticRegression(max_iter=2000).fit(Xs, a).predict_proba(Xs)[:, 1]
    p = a.mean()
    sw = np.where(a == 1, p / ps, (1 - p) / (1 - ps))
    lo, hi = np.percentile(sw, [1, 99]); sw = np.clip(sw, lo, hi)
    out = []
    for h in (60, 120):
        r = {}
        for g in (0, 1):
            m = a == g
            u, S, c1, c2 = weighted_aj(t[m], ev[m], sw[m])
            r[g] = (at(u, c1, h), at(u, c2, h))
        out += [100 * (r[1][0] - r[0][0]), 100 * (r[1][1] - r[0][1])]
    return out   # cvd60, other60, cvd120, other120 (risk differences, percentage points)
rng = np.random.default_rng(1)
n = len(d)
point = estimate(np.arange(n))
boot = np.array([estimate(rng.integers(0, n, n)) for _ in range(500)])
lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
names = ["NegCtrl death @60mo", "All other deaths @60mo", "NegCtrl death @120mo", "All other deaths @120mo"]
print("Risk difference (ACEi/ARB minus untreated), percentage points, IPTW only, B=500:")
for i, nm in enumerate(names):
    print(f"{nm:20s} {point[i]:+.2f}  (95% CI {lo[i]:+.2f} to {hi[i]:+.2f})")

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_pd.DataFrame({"quantity": names, "estimate": list(point), "ci_lo": list(lo), "ci_hi": list(hi)}).to_csv("data/results/negative_control_" + ("all" if MAXD == float("inf") else f"maxdays{int(MAXD)}") + ".csv", index=False)
