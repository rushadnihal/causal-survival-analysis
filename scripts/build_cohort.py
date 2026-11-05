import pandas as pd, numpy as np

from tte.mortality import read_mortality
R = "data/raw/"
CYCLES = {"D": "2005_2006", "E": "2007_2008", "F": "2009_2010", "G": "2011_2012",
          "H": "2013_2014", "I": "2015_2016", "J": "2017_2018"}
ACEI_ARB = r"PRIL\b|SARTAN\b"
def load(f):
    return pd.read_sas(R + f + ".xpt", encoding="latin-1")
parts = []
for suf, span in CYCLES.items():
    df = load(f"DEMO_{suf}")[["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "INDFMPIR",
                              "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]]
    for t in (load(f"GHB_{suf}"),
              load(f"BIOPRO_{suf}")[["SEQN", "LBXSCR"]],
              load(f"DIQ_{suf}")[["SEQN", "DIQ010"]],
              load(f"MCQ_{suf}")[["SEQN", "MCQ160B"]]):
        df = df.merge(t, on="SEQN", how="left")

    female = df.RIAGENDR == 2
    ratio = df.LBXSCR / np.where(female, 0.7, 0.9)
    df["egfr"] = (142 * np.minimum(ratio, 1) ** np.where(female, -0.241, -0.302)
                  * np.maximum(ratio, 1) ** -1.200 * 0.9938 ** df.RIDAGEYR
                  * np.where(female, 1.012, 1.0))

    df = df[(df.RIDAGEYR >= 45)
            & ((df.LBXGH >= 6.5) | (df.DIQ010 == 1))
            & df.egfr.between(30, 89.999)
            & (df.MCQ160B != 1)].copy()

    rx = load(f"RXQ_RX_{suf}")[["SEQN", "RXDDRUG", "RXDDAYS"]]
    rx["RXDDRUG"] = rx["RXDDRUG"].astype(str).str.upper()
    m = rx[rx.RXDDRUG.str.contains(ACEI_ARB, regex=True)]
    df["acei_arb"] = df.SEQN.isin(m.SEQN).astype(int)
    df["acei_arb_days"] = df.SEQN.map(m.groupby("SEQN").RXDDAYS.min())

    mort = read_mortality(f"{R}NHANES_{span}_MORT_2019_PUBLIC.dat")
    df["SEQN"] = df.SEQN.astype(int)
    df = df.merge(mort, on="SEQN", how="left")
    df["cycle"] = span
    parts.append(df)
c = pd.concat(parts, ignore_index=True)
n0 = len(c)
c = c[(c.ELIGSTAT == 1) & c.PERMTH_EXM.notna()].copy()
cvd = (c.MORTSTAT == 1) & c.UCOD_LEADING.isin([1, 5])
c["event"] = np.where(cvd, 1, np.where(c.MORTSTAT == 1, 2, 0))  # 0 alive, 1 CVD, 2 other
c.to_csv("data/processed/cohort.csv", index=False)
print("Eligible before mortality linkage:", n0, "| with usable follow-up:", len(c))
print(c.groupby("cycle").agg(n=("SEQN", "size"), acei_arb=("acei_arb", "sum"),
                             cvd_deaths=("event", lambda e: (e == 1).sum()),
                             other_deaths=("event", lambda e: (e == 2).sum())))
print("\nBy treatment group:")
print(c.groupby("acei_arb").agg(n=("SEQN", "size"),
                                cvd_deaths=("event", lambda e: (e == 1).sum()),
                                other_deaths=("event", lambda e: (e == 2).sum()),
                                median_followup_mo=("PERMTH_EXM", "median")))
