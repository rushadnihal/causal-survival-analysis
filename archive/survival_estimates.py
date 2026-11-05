import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter
def weighted_aj(t, ev, w):
    """Weighted Aalen-Johansen. ev: 0 censored, 1 CVD death, 2 other death.
    Returns times, S(t), CIF1(t), CIF2(t) as step functions at event times."""
    o = np.argsort(t, kind="stable"); t, ev, w = t[o], ev[o], w[o]
    u, idx = np.unique(t, return_index=True)
    at_risk = np.cumsum(w[::-1])[::-1][idx]
    d1 = np.add.reduceat(w * (ev == 1), idx)
    d2 = np.add.reduceat(w * (ev == 2), idx)
    S, c1, c2 = 1.0, 0.0, 0.0
    S_out, c1_out, c2_out = [], [], []
    for n, a, b in zip(at_risk, d1, d2):
        c1 += S * a / n; c2 += S * b / n; S *= 1 - (a + b) / n
        S_out.append(S); c1_out.append(c1); c2_out.append(c2)
    return u, np.array(S_out), np.array(c1_out), np.array(c2_out)
def at(times, vals, x):
    i = np.searchsorted(times, x, side="right") - 1
    return 0.0 if i < 0 else vals[i]
if __name__ == "__main__":
    d = pd.read_csv("data/processed/cohort_weighted.csv")
    t, ev, A = d.PERMTH_EXM.values, d.event.values, d.acei_arb.values

    # ---- validation 1: no competing events -> CIF == 1 - KM (lifelines)
    ev_nc = np.where(ev == 2, 0, ev)
    u, S, c1, c2 = weighted_aj(t, ev_nc, np.ones(len(t)))
    km = KaplanMeierFitter().fit(t, ev_nc == 1)
    gap = max(abs(at(u, c1, x) - (1 - km.predict(x))) for x in (24, 60, 120))
    print("Validation 1 (CIF vs 1-KM, max gap):", round(gap, 10))
    # ---- validation 2: probabilities sum to 1
    u, S, c1, c2 = weighted_aj(t, ev, d.sw_trunc.values)
    print("Validation 2 (max |S+CIF1+CIF2-1|):", round(np.abs(S + c1 + c2 - 1).max(), 10))

    for label, w in [("IPTW only", d.sw_trunc.values),
                     ("IPTW x survey weight", (d.sw_trunc * d.svy_w).values)]:
        print(f"\n=== {label} ===")
        res = {}
        for a in (0, 1):
            m = A == a
            u, S, c1, c2 = weighted_aj(t[m], ev[m], w[m])
            res[a] = {h: (at(u, c1, h), at(u, c2, h)) for h in (60, 120)}
        for h in (60, 120):
            u0, u1 = res[0][h], res[1][h]
            print(f"{h} months | CVD-death risk: untreated {u0[0]*100:.1f}% , ACEi/ARB {u1[0]*100:.1f}% , "
                  f"difference {100*(u1[0]-u0[0]):+.1f} pts | other-death risk: {u0[1]*100:.1f}% vs {u1[1]*100:.1f}%")
