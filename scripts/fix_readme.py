"""Replace the outdated 'treated <= 2 years' paragraph (computed before the medication flags were added) in the README.

Run from the project root:  python3 scripts/fix_readme.py
"""
import os

path = "README.md" if os.path.exists("README.md") else "README_new.md"
s = open(path).read()

old1 = ("**Restricting treated people to <= 2 years on the drug** (477 treated): other-cause death at 60 months moved to "
        "-3.01 (-6.60, +1.15), CVD death +0.13 (-2.89, +3.13): consistent with prevalent-user bias, not proof of it.")
new1 = ("**Restricting treated people to <= 2 years on the drug** (477 treated, medication-adjusted model): CVD death "
        "-0.15 (-3.18, +2.58) at 60 months; other-cause death -3.28 (-7.36, +0.90) at 60 months and +2.20 (-4.79, +8.63) at 120 months. "
        "The other-cause estimate did not move toward zero (it is -3.32 in the full cohort), so duration of use does not explain the "
        "bias. Unmeasured confounding such as frailty or quality of care is the more likely source, but these data cannot confirm that.")

old2 = "This is the signature of healthy-user / prevalent-user bias;"
new2 = "This is the signature of healthy-user bias (prevalent-user bias is also possible, but restricting to recent users did not reduce it);"

for old, new in ((old1, new1), (old2, new2)):
    if old in s:
        s = s.replace(old, new)
        print("updated one passage")
    else:
        print("WARNING: passage not found (already fixed or README edited):", old[:60])
open(path, "w").write(s)
print("wrote", path)
