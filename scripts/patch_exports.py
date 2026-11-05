"""One-time patch: make each analysis script also write its estimates to data/results/*.csv.

Run from the project root:  python3 scripts/patch_exports.py
Safe to run twice (scripts that already carry the export block are skipped).
"""
MARK = "# --- export results"

STD = '''

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_pd.DataFrame({"quantity": names, "estimate": list(point), "ci_lo": list(lo), "ci_hi": list(hi)}).to_csv(%s, index=False)
'''
MAXD_TAG = '("all" if MAXD == float("inf") else f"maxdays{int(MAXD)}")'
STD_TARGETS = {
    "scripts/bootstrap.py": '"data/results/iptw_bootstrap_" + ' + MAXD_TAG + ' + ".csv"',
    "scripts/nc_bootstrap.py": '"data/results/negative_control_" + ' + MAXD_TAG + ' + ".csv"',
    "scripts/survey_bootstrap.py": '"data/results/survey_bootstrap.csv"',
    "scripts/gformula.py": '"data/results/gformula.csv"',
    "scripts/rmst_bootstrap.py": '"data/results/rmst.csv"',
}

WEIGHTS = '''

# --- export results
import os as _os
_os.makedirs("data/results", exist_ok=True)
bal[["covariate", "smd_raw", "smd_weighted"]].to_csv("data/results/balance.csv", index=False)
d[["SEQN", "acei_arb", "ps", "sw", "sw_trunc"]].to_csv("data/results/propensity_scores.csv", index=False)
'''

SIMULATE = '''

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
'''

SENS_APPEND = '''

# --- export results
import os as _os, pandas as _pd
_os.makedirs("data/results", exist_ok=True)
_cols = [f"{l}_{s}" for l in ["cvd60", "other60", "cvd120", "other120"] for s in ("estimate", "ci_lo", "ci_hi")]
_pd.DataFrame(SENS, columns=["variant"] + _cols).to_csv("data/results/sensitivity.csv", index=False)
'''


def patch(path, fn):
    s = open(path).read()
    if MARK in s:
        print("skip (already patched):", path); return
    out = fn(s)
    open(path, "w").write(out)
    print("patched:", path)


def append(block):
    return lambda s: s.rstrip("\n") + block


def sens(s):
    a = "rng = np.random.default_rng(5)\n"
    b = '    print(f"{name:28s}" + "".join('
    assert s.count(a) == 1 and s.count(b) == 1, "sensitivity.py is not in the expected form"
    s = s.replace(a, a + "SENS = []\n")
    s = s.replace(b, "    SENS.append([name] + [x for i in range(4) for x in (pt[i], lo[i], hi[i])])\n" + b)
    return s.rstrip("\n") + SENS_APPEND


for path, target in STD_TARGETS.items():
    patch(path, append(STD % target))
patch("scripts/weights.py", append(WEIGHTS))
patch("scripts/simulate.py", append(SIMULATE))
patch("scripts/sensitivity.py", sens)
