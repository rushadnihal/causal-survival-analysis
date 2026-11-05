import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from tte.psmodel import build_X
from tte.survival import weighted_aj, at
d = pd.read_csv("data/processed/cohort_weighted.csv")
Xdf = build_X(d); X = Xdf.values.astype(float); n = len(d)
Xs = StandardScaler().fit_transform(X)
C = d.PERMTH_EXM.clip(lower=0.5).values          # real follow-up = administrative censoring
# true treatment model
ps_true = LogisticRegression(max_iter=2000).fit(Xs, d.acei_arb.values).predict_proba(Xs)[:, 1]

# true outcome model: exponential cause-specific hazards; covariate effects from a penalized Cox fit
df = pd.DataFrame(Xs, columns=[f"x{i}" for i in range(Xs.shape[1])]); df["T"] = C
coefs, base = {}, {}
import sys
CONF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
lg = np.log(ps_true / (1 - ps_true)); Z = (lg - lg.mean()) / lg.std()
is_cycle = np.array([c.startswith("cycle_") for c in Xdf.columns])
for k in (1, 2):
    df["E"] = (d.event.values == k).astype(int)
    b = CoxPHFitter(penalizer=0.1).fit(df, "T", "E").params_.values.copy()
    b[is_cycle] = 0.0            # keeps censoring (set by cycle) independent of outcome
    coefs[k] = b
    base[k] = df.E.sum() / (np.exp(Xs @ b + CONF * Z) * C).sum()
GAMMA = {1: np.log(0.8), 2: 0.0}   # TRUE effects: 20% lower CVD hazard, no effect on other deaths
def rates(a):
    return {k: base[k] * np.exp(Xs @ coefs[k] + CONF * Z + GAMMA[k] * a) for k in (1, 2)}
# exact truth for exponential competing risks, averaged over the cohort
def true_cif(a, h):
    r = rates(np.full(n, a)); lam = r[1] + r[2]
    return [np.mean(r[k] / lam * (1 - np.exp(-lam * h))) for k in (1, 2)]
truth = {(h, k): 100 * (true_cif(1, h)[k - 1] - true_cif(0, h)[k - 1]) for h in (60, 120) for k in (1, 2)}
def estimate(A, t, ev, weighted):
    if weighted:
        Xi = StandardScaler().fit_transform(X)
        ps = LogisticRegression(max_iter=2000).fit(Xi, A).predict_proba(Xi)[:, 1]
        p = A.mean(); w = np.where(A == 1, p / ps, (1 - p) / (1 - ps))
        lo, hi = np.percentile(w, [1, 99]); w = np.clip(w, lo, hi)
    else:
        w = np.ones(n)
    out = {}
    for h in (60, 120):
        r = {}
        for g in (0, 1):
            m = A == g
            u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
            r[g] = (at(u, c1, h), at(u, c2, h))
        out[(h, 1)] = 100 * (r[1][0] - r[0][0]); out[(h, 2)] = 100 * (r[1][1] - r[0][1])
    return out
rng = np.random.default_rng(7)
res = {"naive": [], "iptw": []}
for _ in range(200):
    A = rng.binomial(1, ps_true)
    r = rates(A)
    t1, t2 = rng.exponential(1 / r[1]), rng.exponential(1 / r[2])
    t = np.minimum(t1, t2); ev = np.where(t1 < t2, 1, 2)
    ev = np.where(t <= C, ev, 0); t = np.minimum(t, C)
    res["naive"].append(estimate(A, t, ev, False)); res["iptw"].append(estimate(A, t, ev, True))
print("Risk difference (treated - untreated), percentage points; 200 simulated datasets")
print(f"{'outcome @ horizon':22s}{'TRUE':>8s}{'naive mean':>12s}{'naive bias':>12s}{'IPTW mean':>11s}{'IPTW bias':>11s}{'IPTW RMSE':>11s}")
for (h, k), tv in truth.items():
    nv = np.array([x[(h, k)] for x in res["naive"]]); ip = np.array([x[(h, k)] for x in res["iptw"]])
    print(f"{'CVD' if k == 1 else 'Other'} death @{h}mo".ljust(22),
          f"{tv:8.2f}{nv.mean():12.2f}{nv.mean() - tv:12.2f}{ip.mean():11.2f}{ip.mean() - tv:11.2f}{np.sqrt(((ip - tv) ** 2).mean()):11.2f}")

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_rows = []
for (h, k), tv in truth.items():
    nv = np.array([x[(h, k)] for x in res["naive"]]); ip = np.array([x[(h, k)] for x in res["iptw"]])
    _rows.append({"outcome": "CVD death" if k == 1 else "Other death", "horizon_months": h, "true": tv,
                  "naive_mean": nv.mean(), "naive_bias": nv.mean() - tv,
                  "iptw_mean": ip.mean(), "iptw_bias": ip.mean() - tv, "iptw_rmse": float(np.sqrt(((ip - tv) ** 2).mean()))})
_pd.DataFrame(_rows).to_csv(f"data/results/simulation_conf{CONF}.csv", index=False)
