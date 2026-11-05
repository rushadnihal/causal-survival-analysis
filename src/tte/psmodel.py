import numpy as np, pandas as pd
def build_X(d):
    num = ["RIDAGEYR", "INDFMPIR", "BMXBMI", "sbp", "egfr", "LBXGH", "log_acr", "tchol", "hdl"]
    binv = ["female", "insulin", "htn_hx", "statin", "metformin", "aspirin", "beta_blocker", "diuretic", "sulfonylurea"]
    miss = ["INDFMPIR_miss", "BMXBMI_miss", "sbp_miss", "log_acr_miss", "htn_hx_miss"]
    cat = pd.get_dummies(d[["RIDRETH1", "smoke", "cycle"]].astype(str), drop_first=True).astype(int)
    return pd.concat([d[num + binv + miss], cat], axis=1)
