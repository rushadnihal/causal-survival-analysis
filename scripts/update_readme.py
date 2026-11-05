"""Append the multiple-imputation, time-varying simulation and dashboard sections to README.md.
Run once from the project root: python3 scripts/update_readme.py  (safe to re-run)."""
from pathlib import Path

MARK = "<!-- extra-sections -->"
TEXT = f"""
{MARK}
## Additional analyses

### Multiple imputation (`scripts/multiple_imputation.py`)
Urine ACR is 27% missing (income 10%, SBP 5%, BMI 2.5%). The main analysis fills missing values with the
median plus missing flags. As a check, values were imputed 10 times (chained regression using treatment,
outcomes and follow-up time as auxiliary variables), IPTW was rerun in each, and results were pooled with
Rubin's rules (100 bootstrap replicates per imputation).

| Risk difference, ACEi/ARB minus untreated (points) | Multiple imputation | Median-fill |
|---|---|---|
| CVD death @60 mo | +0.69 (-1.64, +3.03) | +0.74 (-1.42, +2.74) |
| Other death @60 mo | -3.89 (-7.56, -0.23) | -3.32 (-7.12, +0.01) |
| CVD death @120 mo | -2.85 (-7.24, +1.54) | -2.31 (-6.66, +1.94) |
| Other death @120 mo | -2.27 (-7.77, +3.22) | -1.77 (-7.98, +3.79) |

Missing-data handling does not change the conclusions. The apparent protection against other-cause death is
the same residual healthy-user bias seen elsewhere, not a drug effect.

### Time-varying treatment-confounder feedback (`scripts/simulate_tv.py`)
Simulation with a known truth: sicker people are more likely to be treated, and treatment lowers later
severity. True risk difference at 6 intervals: -12.76 points (5,000 people, 200 replications).

| Method | Mean | Bias | RMSE |
|---|---|---|---|
| Unadjusted | -0.63 | +12.13 | 12.26 |
| Adjusted for current severity | -5.66 | +7.10 | 7.38 |
| IPTW marginal structural model | -10.45 | +2.31 | 3.29 |

Conventional adjustment fails under feedback; weighting recovers most of the effect. This supports the
weighting machinery but does not remove the single-baseline limitation of NHANES.

### Dashboard
```
python3 -m pip install streamlit
python3 -m streamlit run scripts/dashboard.py
```
Four tabs (Summary, Effect estimates, Bias diagnostics, Simulations) read `data/results/*.csv` and `figures/`.
"""

p = Path("README.md")
s = p.read_text()
if MARK in s:
    print("README already updated")
else:
    p.write_text(s.rstrip() + "\n" + TEXT)
    print("README updated")
