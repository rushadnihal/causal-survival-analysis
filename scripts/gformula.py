"""Parametric g-formula for cumulative incidence with a competing event (pooled logistic, 6-month intervals).

Cause-specific discrete-time hazards h1 (CVD death) and h2 (other death) are fit as pooled logistic regressions on
person-intervals with: baseline covariates, treatment, a cubic spline in interval index, and treatment x spline terms.
Each person's risk is then predicted under A = 1 and A = 0 and averaged over the cohort (standardization):
    S_k = S_{k-1} (1 - h1_k - h2_k),   CIF1_K = sum_{k<=K} S_{k-1} h1_k,   CIF2_K likewise.
Assumes no treatment x covariate interactions and correct specification of the hazard models.

Run from the project root:  python3 scripts/gformula.py [B]     (default B = 200; try B = 20 first to check the timing)
"""
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, SplineTransformer

from tte.psmodel import build_X

warnings.filterwarnings("ignore")
B = int(sys.argv[1]) if len(sys.argv) > 1 else 200
W, K = 6, 20                                          # 6-month intervals, 20 intervals = 120 months

d = pd.read_csv("data/processed/cohort_weighted.csv")
X = StandardScaler().fit_transform(build_X(d).values.astype(float))
A, T, EV = d.acei_arb.values, d.PERMTH_EXM.values, d.event.values
spl = SplineTransformer(n_knots=4, degree=3, include_bias=False).fit(np.arange(1, K + 1).reshape(-1, 1))


def design(Xrows, a, k):
    Bk = spl.transform(k.reshape(-1, 1))
    a = a.reshape(-1, 1)
    return np.hstack([Xrows, a, Bk, a * Bk])


def person_period(idx):
    t, ev = T[idx], EV[idx]
    n_orig = np.where(ev > 0, np.maximum(np.ceil(t / W), 1), np.floor(t / W)).astype(int)
    occurs = (ev > 0) & (n_orig <= K)                 # events after 120 months are treated as censored
    n_int = np.minimum(n_orig, K)
    rows = np.repeat(np.arange(len(idx)), n_int)
    k = np.arange(len(rows)) - np.repeat(np.cumsum(n_int) - n_int, n_int) + 1
    last = k == np.repeat(n_int, n_int)
    y1 = (last & np.repeat(occurs & (ev == 1), n_int)).astype(int)
    y2 = (last & np.repeat(occurs & (ev == 2), n_int)).astype(int)
    return rows, k, y1, y2


def risks(idx):
    rows, k, y1, y2 = person_period(idx)
    Xi, Ai = X[idx], A[idx]
    Z = design(Xi[rows], Ai[rows], k)
    m1 = LogisticRegression(max_iter=300).fit(Z, y1)
    m2 = LogisticRegression(max_iter=300).fit(Z, y2)
    n = len(idx)
    rep, ks = np.repeat(np.arange(n), K), np.tile(np.arange(1, K + 1), n)
    out = {}
    for a in (0, 1):
        Zp = design(Xi[rep], np.full(len(rep), a), ks)
        h1 = m1.predict_proba(Zp)[:, 1].reshape(n, K)
        h2 = m2.predict_proba(Zp)[:, 1].reshape(n, K)
        S_prev = np.cumprod(np.hstack([np.ones((n, 1)), np.clip(1 - h1 - h2, 0, 1)]), axis=1)[:, :-1]
        out[a] = (np.cumsum(S_prev * h1, axis=1).mean(0), np.cumsum(S_prev * h2, axis=1).mean(0))
    return out


def summary(res):
    i60, i120 = 60 // W - 1, 120 // W - 1
    return [100 * (res[1][0][i60] - res[0][0][i60]), 100 * (res[1][1][i60] - res[0][1][i60]),
            100 * (res[1][0][i120] - res[0][0][i120]), 100 * (res[1][1][i120] - res[0][1][i120])]


n = len(d)
res = risks(np.arange(n))
i60, i120 = 60 // W - 1, 120 // W - 1
print("g-formula standardized risks (untreated vs ACEi/ARB):")
for nm, j in (("CVD death", 0), ("Other death", 1)):
    print(f"  {nm:12s} 60 mo: {100*res[0][j][i60]:.1f}% vs {100*res[1][j][i60]:.1f}%   |   "
          f"120 mo: {100*res[0][j][i120]:.1f}% vs {100*res[1][j][i120]:.1f}%")
point = summary(res)

rng = np.random.default_rng(21)
boot = np.array([summary(risks(rng.integers(0, n, n))) for _ in range(B)])
lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
names = ["CVD death @60mo", "Other death @60mo", "CVD death @120mo", "Other death @120mo"]
print(f"\nRisk difference (ACEi/ARB minus untreated), percentage points; g-formula; bootstrap B={B}")
for i, nm in enumerate(names):
    print(f"{nm:20s} {point[i]:+.2f}  (95% CI {lo[i]:+.2f} to {hi[i]:+.2f})")

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_pd.DataFrame({"quantity": names, "estimate": list(point), "ci_lo": list(lo), "ci_hi": list(hi)}).to_csv("data/results/gformula.csv", index=False)
