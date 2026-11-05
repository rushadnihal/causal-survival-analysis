import pandas as pd, numpy as np
R = "data/raw/"
SUF = {"2005_2006": "D", "2007_2008": "E", "2009_2010": "F", "2011_2012": "G",
       "2013_2014": "H", "2015_2016": "I", "2017_2018": "J"}
def get(f):
    d = pd.read_sas(R + f + ".xpt", encoding="latin-1")
    d["SEQN"] = d.SEQN.astype(int)
    return d
def pick(d, names):
    for n in names:
        if n in d.columns:
            return d[n]
    return pd.Series(np.nan, index=d.index)
c = pd.read_csv("data/processed/cohort.csv")
parts = []
for span, s in SUF.items():
    sub = c[c.cycle == span].copy()
    bmx = get(f"BMX_{s}")[["SEQN", "BMXBMI"]]
    smq = get(f"SMQ_{s}")
    smq = smq[["SEQN"]].assign(SMQ020=pick(smq, ["SMQ020"]), SMQ040=pick(smq, ["SMQ040"]))
    alb = get(f"ALB_CR_{s}")
    alb = alb[["SEQN"]].assign(acr=pick(alb, ["URDACT"]))
    hdl = get(f"HDL_{s}")
    hdl = hdl[["SEQN"]].assign(hdl=pick(hdl, ["LBDHDD", "LBXHDD"]))
    tc = get(f"TCHOL_{s}")
    tc = tc[["SEQN"]].assign(tchol=pick(tc, ["LBXTC"]))
    bpq = get(f"BPQ_{s}")
    bpq = bpq[["SEQN"]].assign(htn_hx=pick(bpq, ["BPQ020"]))
    bp = get(f"BPXO_{s}" if s == "J" else f"BPX_{s}")
    cols = [f"BPXOSY{i}" for i in (1, 2, 3)] if s == "J" else [f"BPXSY{i}" for i in (1, 2, 3, 4)]
    bp["sbp"] = bp[[x for x in cols if x in bp.columns]].replace(0, np.nan).mean(axis=1)
    bp = bp[["SEQN", "sbp"]]
    rx = get(f"RXQ_RX_{s}")
    rx["ins"] = rx.RXDDRUG.astype(str).str.upper().str.contains("INSULIN")
    ins = rx.groupby("SEQN").ins.max().astype(int).rename("insulin").reset_index()
    for t in (bmx, smq, alb, hdl, tc, bpq, bp, ins):
        sub = sub.merge(t, on="SEQN", how="left")
    parts.append(sub)
d = pd.concat(parts, ignore_index=True)
d["insulin"] = d.insulin.fillna(0)
d["female"] = (d.RIAGENDR == 2).astype(int)
d["smoke"] = np.select([d.SMQ020 == 2, (d.SMQ020 == 1) & d.SMQ040.isin([1, 2]),
                        (d.SMQ020 == 1) & (d.SMQ040 == 3)], [0, 2, 1], default=np.nan)
d["htn_hx"] = d.htn_hx.map({1: 1, 2: 0})
d.to_csv("data/processed/cohort_full.csv", index=False)
vars_ = ["BMXBMI", "sbp", "smoke", "acr", "hdl", "tchol", "htn_hx", "INDFMPIR", "insulin"]
print("Missing share:\n", d[vars_].isna().mean().round(3).to_string())
print("\nBaseline by treatment:")
print(d.groupby("acei_arb").agg(n=("SEQN", "size"), age=("RIDAGEYR", "mean"),
      female=("female", "mean"), bmi=("BMXBMI", "mean"), sbp=("sbp", "mean"),
      egfr=("egfr", "mean"), hba1c=("LBXGH", "mean"), acr_median=("acr", "median"),
      insulin=("insulin", "mean"), htn_hx=("htn_hx", "mean"),
      current_smoker=("smoke", lambda x: (x == 2).mean())).round(2).T)
