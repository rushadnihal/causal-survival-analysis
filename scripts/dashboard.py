"""Streamlit dashboard for the target-trial-emulation project.

Run from the project root:
    python3 -m pip install streamlit
    python3 -m streamlit run scripts/dashboard.py
"""
from pathlib import Path
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "data" / "results"
FIG = ROOT / "figures"

st.set_page_config(page_title="ACEi/ARB target trial emulation", layout="wide")
st.title("ACEi/ARB and mortality in type 2 diabetes with CKD")
st.caption("NHANES 2005-2018 linked to NCHS mortality. Prevalent-user contrast; "
           "results are exploratory, not causal conclusions.")


def load(name):
    p = RES / f"{name}.csv"
    return pd.read_csv(p) if p.exists() else None


def show_csv(name, note=None):
    df = load(name)
    st.subheader(name.replace("_", " "))
    if note:
        st.caption(note)
    if df is None:
        st.warning(f"{name}.csv not found. Run `bash scripts/run_all.sh` first.")
        return
    st.dataframe(df, use_container_width=True)


def show_fig(name, caption):
    p = FIG / name
    if p.exists():
        st.image(str(p), caption=caption, use_container_width=True)
    else:
        st.warning(f"{name} not found in figures/")


tab1, tab2, tab3, tab4 = st.tabs(["Summary", "Effect estimates", "Bias diagnostics", "Simulations"])

with tab1:
    st.markdown("""
**Question.** Does current ACEi/ARB use lower cardiovascular death among adults with type 2 diabetes and CKD?

**Take-away.** The CVD effect is null and inconclusive: intervals span zero at 60 and 120 months.
Other-cause mortality looks lower in users too, which should not be caused by the drug, so residual
healthy-user bias (probably unmeasured confounding) remains. Treat all estimates as associations.

**Robustness.** Median-fill and multiple imputation agree within 0.6 points; survey, g-formula and IPTW
estimates tell the same story.
""")
    c1, c2 = st.columns(2)
    with c1:
        show_fig("consort.png", "Cohort flow")
    with c2:
        show_fig("dag.png", "Causal diagram")

with tab2:
    st.markdown("Risk difference, ACEi/ARB minus untreated (percentage points unless stated).")
    show_csv("iptw_bootstrap_all", "IPTW with bootstrap CIs (PS refit each replicate)")
    show_csv("survey_bootstrap", "Rao-Wu survey bootstrap")
    show_csv("gformula", "Parametric g-formula")
    show_csv("rmst", "Restricted mean survival time / months lost")
    show_csv("multiple_imputation", "Multiple imputation pooled with Rubin's rules")
    cur = load("curves")
    if cur is not None:
        st.subheader("Cumulative incidence curves")
        num = cur.select_dtypes("number").columns.tolist()
        if len(num) >= 2:
            x = st.selectbox("x axis", num, index=0)
            ys = st.multiselect("series", [c for c in num if c != x],
                                default=[c for c in num if c != x][:4])
            if ys:
                st.line_chart(cur.set_index(x)[ys])
        with st.expander("Curve data"):
            st.dataframe(cur, use_container_width=True)

with tab3:
    c1, c2 = st.columns(2)
    with c1:
        show_fig("love_plot.png", "Covariate balance (SMD) before and after weighting")
    with c2:
        show_fig("ps_overlap.png", "Propensity score overlap")
    show_csv("balance")
    show_csv("negative_control_all", "Negative control outcome (cancer, lung, accident deaths): "
             "a non-zero effect signals residual bias")
    show_csv("sensitivity", "Sensitivity analyses (trimming, overlap weights, follow-up restriction)")
    show_csv("iptw_bootstrap_maxdays730", "Restricted to short treatment duration (<=730 days)")

with tab4:
    st.markdown("Estimator checks against a known truth (true CVD hazard ratio 0.8).")
    sims = sorted(p.stem for p in RES.glob("simulation*.csv"))
    if sims:
        for s in sims:
            show_csv(s)
    else:
        st.warning("No simulation results found.")
    st.markdown("""
**Time-varying feedback (`simulation_tv`).** Unadjusted and severity-adjusted estimates are heavily biased
when treatment lowers later severity; the IPTW marginal structural model recovers most of the true effect.
""")
