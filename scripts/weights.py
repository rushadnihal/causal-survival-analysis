import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
d = pd.read_csv("data/processed/cohort_full.csv")
MEDS = ["statin", "metformin", "aspirin", "beta_blocker", "diuretic", "sulfonylurea"]
d["log_acr"] = np.log1p(d.acr)
num = ["RIDAGEYR", "INDFMPIR", "BMXBMI", "sbp", "egfr", "LBXGH", "log_acr", "tchol", "hdl"]
binv = ["female", "insulin", "htn_hx"] + MEDS
miss = []
for v in ["INDFMPIR", "BMXBMI", "sbp", "log_acr", "htn_hx"]:
    d[v + "_miss"] = d[v].isna().astype(int); miss.append(v + "_miss")
for v in num + ["htn_hx"]:
    d[v] = d[v].fillna(d[v].median())
cat = pd.get_dummies(d[["RIDRETH1", "smoke", "cycle"]].astype(str), drop_first=True).astype(int)
X = pd.concat([d[num + binv + miss], cat], axis=1)
A = d.acei_arb.values
Xs = StandardScaler().fit_transform(X)
ps = LogisticRegression(max_iter=2000).fit(Xs, A).predict_proba(Xs)[:, 1]
p = A.mean()
sw = np.where(A == 1, p / ps, (1 - p) / (1 - ps))
lo, hi = np.percentile(sw, [1, 99])
d["ps"] = ps
d["sw"] = sw
d["sw_trunc"] = np.clip(sw, lo, hi)
d["svy_w"] = d.WTMEC2YR / 7          # 7 pooled cycles
def wmean_var(x, w):
    m = np.average(x, weights=w)
    return m, np.average((x - m) ** 2, weights=w)
def smd(x, w):
    m1, v1 = wmean_var(x[A == 1], w[A == 1])
    m0, v0 = wmean_var(x[A == 0], w[A == 0])
    return (m1 - m0) / np.sqrt((v1 + v0) / 2)
rows = [(c, smd(X[c].values.astype(float), np.ones(len(X))),
         smd(X[c].values.astype(float), d.sw_trunc.values)) for c in X.columns]
bal = pd.DataFrame(rows, columns=["covariate", "smd_raw", "smd_weighted"])
bal["abs_raw"] = bal.smd_raw.abs(); bal["abs_w"] = bal.smd_weighted.abs()
print("PS range  untreated: %.3f-%.3f | treated: %.3f-%.3f" %
      (ps[A == 0].min(), ps[A == 0].max(), ps[A == 1].min(), ps[A == 1].max()))
print("Stabilized weights: mean %.3f | min %.3f | max %.3f | 1st/99th pct %.3f / %.3f" %
      (sw.mean(), sw.min(), sw.max(), lo, hi))
print("\nLargest imbalances (|SMD|):")
print(bal.sort_values("abs_raw", ascending=False).head(8)[["covariate", "smd_raw", "smd_weighted"]].round(3).to_string(index=False))
print("\nMax |SMD| raw: %.3f | weighted: %.3f | covariates still > 0.10 after weighting: %d"
      % (bal.abs_raw.max(), bal.abs_w.max(), (bal.abs_w > 0.10).sum()))
plt.figure(figsize=(6, 4))
plt.hist(ps[A == 0], bins=30, alpha=.6, density=True, label="Untreated")
plt.hist(ps[A == 1], bins=30, alpha=.6, density=True, label="ACEi/ARB")
plt.xlabel("Propensity score"); plt.ylabel("Density"); plt.legend()
plt.title("Propensity score overlap"); plt.tight_layout();
plt.savefig("figures/ps_overlap.png", dpi=150)
b = bal.sort_values("abs_raw")
plt.figure(figsize=(6, 7)); y = np.arange(len(b))
plt.scatter(b.abs_raw, y, label="Unweighted"); plt.scatter(b.abs_w, y, label="Weighted")
plt.yticks(y, b.covariate, fontsize=7); plt.axvline(0.1, ls="--", c="grey")
plt.xlabel("|Standardized mean difference|"); plt.legend(); plt.title("Love plot")
plt.tight_layout(); plt.savefig("figures/love_plot.png", dpi=150)
d.to_csv("data/processed/cohort_weighted.csv", index=False)

# --- export results
import os as _os
_os.makedirs("data/results", exist_ok=True)
bal[["covariate", "smd_raw", "smd_weighted"]].to_csv("data/results/balance.csv", index=False)
d[["SEQN", "acei_arb", "ps", "sw", "sw_trunc"]].to_csv("data/results/propensity_scores.csv", index=False)
