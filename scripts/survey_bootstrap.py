"""Survey-design bootstrap (Rao-Wu rescaled bootstrap) for survey-weighted IPTW cumulative incidence.

Within each stratum (cycle x SDMVSTRA) with n_h >= 2 PSUs, draw n_h - 1 PSUs with replacement and rescale
base survey weights by n_h / (n_h - 1) x (number of times the PSU was drawn). The propensity model is refit
(weighted) on every replicate. Strata with a single observed PSU are left unresampled (they contribute no variance).

Run from the project root:  python3 scripts/survey_bootstrap.py
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from tte.psmodel import build_X
from tte.survival import weighted_aj, at

B = 500
d = pd.read_csv("data/processed/cohort_weighted.csv")
X = build_X(d).values.astype(float)
Xs = StandardScaler().fit_transform(X)
A, T, EV = d.acei_arb.values, d.PERMTH_EXM.values, d.event.values
base_w = (d.WTMEC2YR / 7).values                    # 7 pooled 2-year cycles

stratum = d.cycle.astype(str) + "_" + d.SDMVSTRA.astype(int).astype(str)
psu = stratum + "_" + d.SDMVPSU.astype(int).astype(str)
groups = []                                          # per stratum: list of index arrays, one per PSU
for s in stratum.unique():
    m = (stratum == s).values
    groups.append([np.where(m & (psu == p).values)[0] for p in psu[m].unique()])
single = sum(len(g) == 1 for g in groups)
print(f"{len(groups)} strata, {int(psu.nunique())} PSUs in the cohort, {single} strata with a single observed PSU")


def replicate_weights(rng):
    mult = np.ones(len(d))
    for g in groups:
        nh = len(g)
        if nh < 2:
            continue
        draws = rng.integers(0, nh, nh - 1)
        counts = np.bincount(draws, minlength=nh)
        for idx, c in zip(g, counts):
            mult[idx] = nh / (nh - 1) * c
    return base_w * mult


def estimate(w_svy):
    ps = LogisticRegression(max_iter=2000).fit(Xs, A, sample_weight=w_svy).predict_proba(Xs)[:, 1]
    p = np.average(A, weights=w_svy)
    sw = np.where(A == 1, p / ps, (1 - p) / (1 - ps))
    lo, hi = np.percentile(sw, [1, 99]); sw = np.clip(sw, lo, hi)
    w = sw * w_svy
    out = []
    for h in (60, 120):
        r = {}
        for g in (0, 1):
            m = (A == g) & (w > 0)                   # drop zero-weight rows (unselected PSUs)
            u, S, c1, c2 = weighted_aj(T[m], EV[m], w[m])
            r[g] = (at(u, c1, h), at(u, c2, h))
        out += [100 * (r[1][0] - r[0][0]), 100 * (r[1][1] - r[0][1])]
    return out                                       # cvd60, other60, cvd120, other120


point = estimate(base_w)
rng = np.random.default_rng(11)
boot = np.array([estimate(replicate_weights(rng)) for _ in range(B)])
lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
names = ["CVD death @60mo", "Other death @60mo", "CVD death @120mo", "Other death @120mo"]
print(f"\nSurvey-weighted IPTW risk difference (ACEi/ARB minus untreated), percentage points; Rao-Wu bootstrap B={B}:")
for i, nm in enumerate(names):
    print(f"{nm:20s} {point[i]:+.2f}  (95% CI {lo[i]:+.2f} to {hi[i]:+.2f})")

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_pd.DataFrame({"quantity": names, "estimate": list(point), "ci_lo": list(lo), "ci_hi": list(hi)}).to_csv("data/results/survey_bootstrap.csv", index=False)
