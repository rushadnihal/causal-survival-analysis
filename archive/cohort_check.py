import pandas as pd, numpy as np
R = "data/raw/"
def load(f):
    return pd.read_sas(R + f + ".xpt", encoding="latin-1")
demo = load("DEMO_H")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
df = demo
for t in (load("GHB_H"),
          load("BIOPRO_H")[["SEQN", "LBXSCR"]],
          load("DIQ_H")[["SEQN", "DIQ010"]],
          load("MCQ_H")[["SEQN", "MCQ160B"]]):
    df = df.merge(t, on="SEQN", how="left")
# CKD-EPI 2021 (race-free)
female = df.RIAGENDR == 2
k = np.where(female, 0.7, 0.9)
a = np.where(female, -0.241, -0.302)
ratio = df.LBXSCR / k
df["egfr"] = (142 * np.minimum(ratio, 1) ** a * np.maximum(ratio, 1) ** -1.200
              * 0.9938 ** df.RIDAGEYR * np.where(female, 1.012, 1.0))

print("All participants:", len(df))
df = df[df.RIDAGEYR >= 45];                         print("Age >= 45:", len(df))
df = df[(df.LBXGH >= 6.5) | (df.DIQ010 == 1)];      print("Diabetes (HbA1c>=6.5 or diagnosed):", len(df))
df = df[df.egfr.between(30, 89.999)];               print("eGFR 30-89:", len(df))
df = df[df.MCQ160B != 1];                           print("No heart failure:",
len(df))
# ---- treatment groups ----
rx = load("RXQ_RX_H")[["SEQN", "RXDDRUG", "RXDDAYS"]]
rx["RXDDRUG"] = rx["RXDDRUG"].astype(str)
rx = rx[rx.SEQN.isin(df.SEQN)]
def group(pattern):
    m = rx[rx.RXDDRUG.str.contains(pattern, regex=True)]
    return m.groupby("SEQN").RXDDAYS.min()   # shortest duration among matching drugs
acei_arb = group(r"PRIL\b|SARTAN\b")
ccb = group(r"DIPINE\b|DILTIAZEM|VERAPAMIL")
for name, g in [("ACEi/ARB", acei_arb), ("CCB", ccb)]:
    print(name, "users:", len(g),
          "| duration missing:", g.isna().sum(),
          "| <=180 days:", (g <= 180).sum(),
          "| <=365 days:", (g <= 365).sum())
print("Both ACEi/ARB and CCB:", len(set(acei_arb.index) & set(ccb.index)))
