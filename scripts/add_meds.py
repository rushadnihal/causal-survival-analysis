import pandas as pd
R = "data/raw/"
SUF = {"2005_2006": "D", "2007_2008": "E", "2009_2010": "F", "2011_2012": "G",
       "2013_2014": "H", "2015_2016": "I", "2017_2018": "J"}
PAT = {"statin": r"STATIN\b",
       "metformin": r"METFORMIN",
       "aspirin": r"ASPIRIN",
       "beta_blocker": r"OLOL\b",
       "diuretic": r"THIAZIDE|CHLORTHALIDONE|FUROSEMIDE|SPIRONOLACTONE|TORSEMIDE|BUMETANIDE|INDAPAMIDE",
       "sulfonylurea": r"GLIPIZIDE|GLYBURIDE|GLIMEPIRIDE"}
d = pd.read_csv("data/processed/cohort_full.csv")
parts = []
for span, s in SUF.items():
    sub = d[d.cycle == span].copy()
    rx = pd.read_sas(f"{R}RXQ_RX_{s}.xpt", encoding="latin-1")
    rx["SEQN"] = rx.SEQN.astype(int)
    rx["RXDDRUG"] = rx.RXDDRUG.astype(str).str.upper()
    for name, pat in PAT.items():
        ids = set(rx[rx.RXDDRUG.str.contains(pat, regex=True)].SEQN)
        sub[name] = sub.SEQN.isin(ids).astype(int)
    parts.append(sub)
out = pd.concat(parts, ignore_index=True)
out.to_csv("data/processed/cohort_full.csv", index=False)
print(out.groupby("acei_arb")[list(PAT)].mean().round(3).T)
