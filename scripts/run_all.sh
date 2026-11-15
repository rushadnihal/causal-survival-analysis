#!/usr/bin/env bash
# Re-run the whole analysis from the weighting step onward and write every result to data/results/.
# Prerequisites: raw data downloaded and cohort built (build_cohort / build_covariates / add_meds).
# Takes roughly 20-40 minutes. Run from anywhere:  bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

run() { echo; echo "=== [$1] $2"; shift 2; "$@"; }

run "1/11" "propensity model and balance"        python3 scripts/weights.py
run "2/11" "cumulative incidence curves"         python3 scripts/export_curves.py
run "3/11" "IPTW bootstrap (all treated)"        python3 scripts/bootstrap.py
run "4/11" "IPTW bootstrap (<= 730 days)"        python3 scripts/bootstrap.py 730
run "5/11" "negative control outcome"            python3 scripts/nc_bootstrap.py
run "6/11" "survey-design bootstrap"             python3 scripts/survey_bootstrap.py
run "7/11" "sensitivity analyses"                python3 scripts/sensitivity.py
run "8/11" "restricted mean times"               python3 scripts/rmst_bootstrap.py
run "9/11" "parametric g-formula"                python3 scripts/gformula.py 200
run "10/11" "simulation, no built-in confounding" python3 scripts/simulate.py 0
run "11/11" "simulation, healthy-user confounding" python3 scripts/simulate.py -0.5

echo; echo "Done. Results are in data/results/:"; ls -1 data/results
