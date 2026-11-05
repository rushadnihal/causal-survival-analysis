#!/usr/bin/env bash
# One-time tidy-up: moves the loose scripts into src/ and scripts/, fixes their imports,
# installs the package in editable mode, runs the tests, and re-runs the weighting step.
# Run from the project root:   bash reorganize.sh
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate

mkdir -p src/tte scripts tests archive

# library code
[ -f load_mortality.py ] && mv load_mortality.py src/tte/mortality.py
[ -f psmodel.py ]        && mv psmodel.py src/tte/psmodel.py

# runnable scripts
for f in add_meds bootstrap build_cohort build_covariates download_covariates \
         download_mortality download_nhanes nc_bootstrap simulate weights; do
  [ -f "$f.py" ] && mv "$f.py" scripts/
done

# superseded scratch files (kept, not deleted)
[ -f cohort_check.py ]        && mv cohort_check.py archive/
[ -f survival_estimates.py ]  && mv survival_estimates.py archive/

# imports now come from the installed package
sed -i '' 's/from load_mortality import read_mortality/from tte.mortality import read_mortality/' scripts/build_cohort.py
sed -i '' 's/from psmodel import build_X/from tte.psmodel import build_X/' scripts/bootstrap.py scripts/nc_bootstrap.py scripts/simulate.py
sed -i '' 's/from survival_estimates import weighted_aj, at/from tte.survival import weighted_aj, at/' scripts/bootstrap.py scripts/nc_bootstrap.py scripts/simulate.py

# make the pipeline reproducible from a fresh checkout: medication flags live in cohort_full.csv
sed -i '' 's#cohort_weighted.csv#cohort_full.csv#g' scripts/add_meds.py
sed -i '' '/\[\["SEQN"\] + MEDS\]/d' scripts/weights.py

python3 -m pip install -q -e ".[dev]"
echo "--- tests ---"
python3 -m pytest -q
echo "--- re-running medication flags + weights (should match earlier: max |SMD| weighted 0.051) ---"
python3 scripts/add_meds.py
python3 scripts/weights.py
echo "Done."
