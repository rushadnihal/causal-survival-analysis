"""Multiple imputation sensitivity analysis for the IPTW risk differences.

M completed datasets are drawn with a chained-equations imputer (sklearn IterativeImputer, posterior sampling). The
imputation model uses all covariates, treatment, event indicators and log follow-up time. For each completed dataset the IPTW
estimate is computed with a bootstrap (BW replicates) for its within-imputation variance; results are pooled with Rubin's rules.

Run from the project root:  python3 scripts/multiple_imputation.py [BW]   (default BW = 100 bootstrap reps per imputation)
"""
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge, LogisticRegression
from sklearn.preprocessing import StandardScaler

from tte.survival import weighted_aj, at
from tte.weighting import stabilized_weights

warnings.filterwarnings("ignore")
M = 10
BW = int(sys.argv[1]) if len(sys.argv) > 1 else 100

d = pd.read_csv("data/processed/cohort_full.csv")
d["log_acr"] = np.log1p(d.acr)
MEDS = ["statin", "metformin", "aspirin", "beta_blocker", "diuretic", "sulfonylurea"]
cont = ["RIDAGEYR", "INDFMPIR", "BMXBMI", "sbp", "egfr", "LBXGH", "log_acr", "tchol", "hdl"]
binv = ["female", "insulin", "htn_hx"] + MEDS
A, T, EV = d.acei_arb.values, d.PERMTH_EXM.values, d.event.values

print("Missing share before imputation:")
print(d[["acr", "INDFMPIR", "BMXBMI", "sbp", "htn_hx"]].isna().mean().round(3).to_string())

# imputation matrix: covariates + exposure + outcome information + design dummies
aux = pd.DataFrame({"A": A, "cvd": (EV == 1).astype(int), "oth": (EV == 2).astype(int),
                    "logT": np.log1p(T), "smoke": d.smoke})
dummies = pd.get_dummies(d[["RIDRETH1", "cycle"]].astype(str), drop_first=True).astype(int)
Z = pd.concat([d[cont + binv], aux, dummies], axis=1)
cat = pd.get_dummies(d[["RIDRETH1", "smoke", "cycle"]].astype(str), drop_first=True).astype(int)


def estimate(X, idx):
    a, t, ev = A[idx], T[idx], EV[idx]
    Xs = StandardScaler().fit_transform(X[idx])
    ps = LogisticRegression(max_iter=2000).fit(Xs, a).predict_proba(Xs)[:, 1]
    w = stabilized_weights(ps, a, trunc=(1, 99))
    out = []
    for h in (60, 120):
        r = {}
        for g in (0, 1):
            m = a == g
            u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
            r[g] = (at(u, c1, h), at(u, c2, h))
        out += [100 * (r[1][0] - r[0][0]), 100 * (r[1][1] - r[0][1])]
    return out                                    # cvd60, other60, cvd120, other120


rng = np.random.default_rng(2024)
n = len(d)
points, within = [], []
for m in range(M):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True, max_iter=10,
                           min_value=Z.min().values, max_value=Z.max().values, random_state=m)
    Zi = pd.DataFrame(imp.fit_transform(Z.values), columns=Z.columns)
    Zi["htn_hx"] = (Zi.htn_hx > 0.5).astype(int)
    X = pd.concat([Zi[cont + binv], cat], axis=1).values.astype(float)
    points.append(estimate(X, np.arange(n)))
    boot = np.array([estimate(X, rng.integers(0, n, n)) for _ in range(BW)])
    within.append(boot.var(axis=0, ddof=1))
    print(f"imputation {m + 1}/{M} done")

points, within = np.array(points), np.array(within)
qbar = points.mean(axis=0)                        # pooled estimate
W = within.mean(axis=0)                           # within-imputation variance
Bv = points.var(axis=0, ddof=1)                   # between-imputation variance
Tv = W + (1 + 1 / M) * Bv                         # Rubin's total variance
se = np.sqrt(Tv)
lam = (1 + 1 / M) * Bv / Tv                       # share of variance due to missing data

try:
    ref = pd.read_csv("data/results/iptw_bootstrap_all.csv")
except FileNotFoundError:
    ref = None
names = ["CVD death @60mo", "Other death @60mo", "CVD death @120mo", "Other death @120mo"]
print(f"\nRisk difference (ACEi/ARB minus untreated), percentage points; IPTW; M = {M} imputations (Rubin's rules)")
print(f"{'':20s}{'pooled (95% CI)':>28s}{'median-fill result':>26s}{'missing-data share':>22s}")
for i, nm in enumerate(names):
    ci = f"{qbar[i]:+.2f} ({qbar[i] - 1.96 * se[i]:+.2f}, {qbar[i] + 1.96 * se[i]:+.2f})"
    rf = "" if ref is None else f"{ref.estimate[i]:+.2f} ({ref.ci_lo[i]:+.2f}, {ref.ci_hi[i]:+.2f})"
    print(f"{nm:20s}{ci:>28s}{rf:>26s}{lam[i]:>21.0%}")

import os
os.makedirs("data/results", exist_ok=True)
pd.DataFrame({"quantity": names, "estimate": qbar, "ci_lo": qbar - 1.96 * se, "ci_hi": qbar + 1.96 * se,
              "within_var": W, "between_var": Bv}).to_csv("data/results/multiple_imputation.csv", index=False)
