"""Draw the causal DAG and the cohort attrition (CONSORT-style) flowchart -> figures/dag.png, figures/consort.png

Run from the project root:  python3 scripts/figures.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tte.mortality import read_mortality

R = "data/raw/"
CYCLES = {"D": "2005_2006", "E": "2007_2008", "F": "2009_2010", "G": "2011_2012",
          "H": "2013_2014", "I": "2015_2016", "J": "2017_2018"}


def load(f):
    return pd.read_sas(R + f + ".xpt", encoding="latin-1")


# ---------------------------------------------------------------- attrition counts
counts = dict(all=0, age=0, dm=0, egfr=0, nohf=0, final=0)
for suf, span in CYCLES.items():
    df = load(f"DEMO_{suf}")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
    for t in (load(f"GHB_{suf}"), load(f"BIOPRO_{suf}")[["SEQN", "LBXSCR"]],
              load(f"DIQ_{suf}")[["SEQN", "DIQ010"]], load(f"MCQ_{suf}")[["SEQN", "MCQ160B"]]):
        df = df.merge(t, on="SEQN", how="left")
    female = df.RIAGENDR == 2
    ratio = df.LBXSCR / np.where(female, 0.7, 0.9)
    df["egfr"] = (142 * np.minimum(ratio, 1) ** np.where(female, -0.241, -0.302)
                  * np.maximum(ratio, 1) ** -1.200 * 0.9938 ** df.RIDAGEYR * np.where(female, 1.012, 1.0))
    mort = read_mortality(f"{R}NHANES_{span}_MORT_2019_PUBLIC.dat")
    df["SEQN"] = df.SEQN.astype(int)
    df = df.merge(mort[["SEQN", "ELIGSTAT", "PERMTH_EXM"]], on="SEQN", how="left")
    age = df.RIDAGEYR >= 45
    dm = age & ((df.LBXGH >= 6.5) | (df.DIQ010 == 1))
    eg = dm & df.egfr.between(30, 89.999)
    nohf = eg & (df.MCQ160B != 1)
    fin = nohf & (df.ELIGSTAT == 1) & df.PERMTH_EXM.notna()
    for k, m in zip(counts, (pd.Series(True, index=df.index), age, dm, eg, nohf, fin)):
        counts[k] += int(m.sum())

coh = pd.read_csv("data/processed/cohort.csv")
assert counts["final"] == len(coh), f"flowchart ({counts['final']}) does not match cohort ({len(coh)})"
arm = coh.acei_arb.value_counts()
print("Attrition counts:", counts, "| untreated", arm.get(0, 0), "| ACEi/ARB", arm.get(1, 0))

steps = [(f"NHANES 2005-2018 examined participants\nn = {counts['all']:,}", None),
         (f"Aged 45 or older\nn = {counts['age']:,}", f"Excluded: under 45 (n = {counts['all'] - counts['age']:,})"),
         (f"Diabetes (HbA1c >= 6.5% or diagnosed)\nn = {counts['dm']:,}", f"Excluded: no diabetes / missing (n = {counts['age'] - counts['dm']:,})"),
         (f"eGFR 30-89 mL/min/1.73m2\nn = {counts['egfr']:,}", f"Excluded: eGFR outside 30-89 or no creatinine (n = {counts['dm'] - counts['egfr']:,})"),
         (f"No history of heart failure\nn = {counts['nohf']:,}", f"Excluded: prior heart failure (n = {counts['egfr'] - counts['nohf']:,})"),
         (f"Linked to mortality follow-up (analytic cohort)\nn = {counts['final']:,}", f"Excluded: not eligible for linkage / no follow-up time (n = {counts['nohf'] - counts['final']:,})"),
         (f"ACEi/ARB users n = {arm.get(1, 0):,}      Non-users n = {arm.get(0, 0):,}", None)]
fig, ax = plt.subplots(figsize=(9, 11)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ys = np.linspace(0.93, 0.07, len(steps))
for i, (txt, excl) in enumerate(steps):
    ax.text(0.33, ys[i], txt, ha="center", va="center", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.5", fc="#eef3fa", ec="#2a4d7a"))
    if i > 0:
        ax.annotate("", xy=(0.33, ys[i] + 0.04), xytext=(0.33, ys[i - 1] - 0.04),
                    arrowprops=dict(arrowstyle="->", color="#2a4d7a"))
    if excl:
        ym = (ys[i] + ys[i - 1]) / 2
        ax.annotate("", xy=(0.58, ym), xytext=(0.33, ym), arrowprops=dict(arrowstyle="->", color="grey"))
        ax.text(0.59, ym, excl, ha="left", va="center", fontsize=8.5, color="#444",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey"), wrap=True)
plt.tight_layout(); plt.savefig("figures/consort.png", dpi=150); plt.close()

print("Saved figures/consort.png")
