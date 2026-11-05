"""Simulation B1: time-varying treatment-confounder feedback, marginal structural model.

L_k = time-varying 'illness severity' (high = sicker). Sicker people are more likely to be treated
(confounding by indication) AND treatment lowers later severity (treatment -> confounder feedback).
Conventional adjustment for L_k blocks part of the treatment effect and is biased; the IPTW MSM is not.

Truth: risk difference at K intervals under 'always treat' vs 'never treat' from a large
simulation in which treatment is set by intervention.

Run from project root:  python3 scripts/simulate_tv.py [N] [REPS]
"""
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 200
K = 6
BETA = -0.4  # direct log-odds effect of current treatment on death


def expit(x):
    return 1 / (1 + np.exp(-x))


def simulate(n, rng, intervene=None):
    """Return long person-interval frame. intervene = 0/1 sets treatment for everyone."""
    L = rng.normal(size=n)
    Aprev = np.zeros(n)
    alive = np.ones(n, bool)
    rows = []
    dead_ever = np.zeros(n, bool)
    for k in range(K):
        if intervene is None:
            A = rng.binomial(1, expit(-0.3 + 0.8 * L + 1.0 * Aprev)).astype(float)
        else:
            A = np.full(n, float(intervene))
        D = rng.binomial(1, expit(-3.2 + 0.7 * L + BETA * A)).astype(bool)
        idx = np.where(alive)[0]
        rows.append(pd.DataFrame({"id": idx, "k": k, "L": L[idx], "A": A[idx], "Aprev": Aprev[idx],
                                  "D": D[idx].astype(int)}))
        dead_ever |= alive & D
        alive = alive & ~D
        L = 0.7 * L - 0.4 * A + rng.normal(scale=0.7, size=n)
        Aprev = A
    d = pd.concat(rows, ignore_index=True)
    d["cumA"] = d.groupby("id")["A"].cumsum()
    return d, dead_ever.mean()


def logit_fit(X, y, w=None):
    m = LogisticRegression(C=1e8, max_iter=500, fit_intercept=False)
    m.fit(X, y, sample_weight=w)
    return m


def tdum(d):
    return np.eye(K)[d.k.values]


def risk_from(coef, cum_a, extra=0.0):
    """coef = [K time effects, beta_cumA, ...]. cum_a(k) -> cumulative treatment at k."""
    s = 1.0
    for k in range(K):
        h = expit(coef[k] + coef[K] * cum_a(k) + extra)
        s *= 1 - h
    return 1 - s


def estimate(d):
    out = {}
    T = tdum(d)
    always, never = (lambda k: k + 1.0), (lambda k: 0.0)
    # 1 unadjusted
    m = logit_fit(np.c_[T, d.cumA.values], d.D.values)
    c = m.coef_[0]
    out["unadjusted"] = risk_from(c, always) - risk_from(c, never)
    # 2 conventional: adjust for time-varying L_k, evaluated at L=0
    m = logit_fit(np.c_[T, d.cumA.values, d.L.values], d.D.values)
    c = m.coef_[0]
    out["adjusted_for_L"] = risk_from(c, always) - risk_from(c, never)
    # 3 IPTW MSM with stabilized cumulative weights
    den = logit_fit(np.c_[T, d.L.values, d.Aprev.values], d.A.values)
    num = logit_fit(np.c_[T, d.Aprev.values], d.A.values)
    pd_ = den.predict_proba(np.c_[T, d.L.values, d.Aprev.values])[:, 1]
    pn_ = num.predict_proba(np.c_[T, d.Aprev.values])[:, 1]
    a = d.A.values
    f = np.where(a == 1, pn_ / pd_, (1 - pn_) / (1 - pd_))
    d = d.assign(f=f)
    sw = d.groupby("id")["f"].cumprod().values
    lo, hi = np.percentile(sw, [1, 99])
    sw = np.clip(sw, lo, hi)
    m = logit_fit(np.c_[T, d.cumA.values], d.D.values, sw)
    c = m.coef_[0]
    out["IPTW_MSM"] = risk_from(c, always) - risk_from(c, never)
    return out


rng = np.random.default_rng(2026)
big = 400_000
risk1 = simulate(big, rng, intervene=1)[1]
risk0 = simulate(big, rng, intervene=0)[1]
truth = risk1 - risk0
print(f"True risk @K={K}: always-treat {risk1:.4f}, never-treat {risk0:.4f}, RD {100*truth:+.2f} points")

res = []
for r in range(REPS):
    d, _ = simulate(N, np.random.default_rng(1000 + r))
    res.append(estimate(d))
R = pd.DataFrame(res) * 100
tp = 100 * truth
summ = pd.DataFrame({"mean": R.mean(), "bias": R.mean() - tp, "sd": R.std(),
                     "rmse": np.sqrt(((R - tp) ** 2).mean())}).round(2)
summ["true_RD"] = round(tp, 2)
print(f"\nRisk difference at {K} intervals, percentage points; N={N}, {REPS} replications")
print(summ)
import os
os.makedirs("data/results", exist_ok=True)
summ.to_csv("data/results/simulation_tv.csv")
print("saved data/results/simulation_tv.csv")
