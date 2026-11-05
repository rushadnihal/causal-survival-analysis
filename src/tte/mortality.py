import pandas as pd, glob
COLS = [(0, 14), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
NAMES = ["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES", "HYPERTEN",
         "PERMTH_INT", "PERMTH_EXM"]
def read_mortality(path):
    return pd.read_fwf(path, colspecs=COLS, names=NAMES, na_values=["."])
if __name__ == "__main__":
    for p in sorted(glob.glob("data/raw/NHANES_*_MORT_2019_PUBLIC.dat")):
        m = read_mortality(p)
        e = m[m.ELIGSTAT == 1]
        dead = e[e.MORTSTAT == 1]
        print(p.split("NHANES_")[1][:9],
              "| rows", len(m), "| eligible", len(e), "| deaths", len(dead),
              "| deaths w/o cause", dead.UCOD_LEADING.isna().sum(),
              "| alive w/ cause", e[(e.MORTSTAT == 0) & e.UCOD_LEADING.notna()].shape[0],
              "| max follow-up (mo)", e.PERMTH_EXM.max())
