"""Draw the causal DAG -> figures/dag.png

Run from the project root:  python3 scripts/draw_dag.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))
ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

nodes = {
    "L": (0.50, 0.88, "Measured baseline covariates (L)\nage, sex, race, income, BMI, SBP, eGFR, HbA1c,\n"
                      "urine ACR, lipids, smoking, insulin, co-medications", "-"),
    "U": (0.50, 0.10, "Unmeasured (U)\nfrailty, health-seeking behaviour,\nquality of care", "--"),
    "A": (0.13, 0.50, "ACEi / ARB use\n(A)", "-"),
    "Y": (0.87, 0.66, "Cardiovascular death\n(Y)", "-"),
    "C": (0.87, 0.34, "Other-cause death\n(competing event)", "-"),
}
boxes = {}
for k, (x, y, txt, ls) in nodes.items():
    boxes[k] = ax.text(x, y, txt, ha="center", va="center", fontsize=9.5, zorder=3,
                       bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", ls=ls))


def arrow(a, b, bold=False, dashed=False, rad=0.0):
    (x1, y1, _, _), (x2, y2, _, _) = nodes[a], nodes[b]
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=1, arrowprops=dict(
        arrowstyle="-|>", lw=2.8 if bold else 1.2, color="#b22222" if bold else "#555",
        ls="--" if dashed else "-", shrinkA=0, shrinkB=0,
        connectionstyle=f"arc3,rad={rad}",
        patchA=boxes[a].get_bbox_patch(), patchB=boxes[b].get_bbox_patch()))   # clip at the box outlines


for src, dst in (("L", "A"), ("L", "Y")):
    arrow(src, dst)
arrow("L", "C", rad=-0.25)                      # curve around the right-hand side to limit crossings
for src, dst in (("U", "A"), ("U", "C")):
    arrow(src, dst, dashed=True)
arrow("U", "Y", dashed=True, rad=0.25)
arrow("A", "Y", bold=True)
arrow("A", "C")

ax.text(0.32, 0.585, "effect of interest", color="#b22222", fontsize=9, ha="center", rotation=7)
ax.text(0.5, 0.01, "Weighting adjusts for L. The dashed paths through U cannot be blocked with these data.",
        ha="center", fontsize=9, style="italic")
plt.tight_layout(); plt.savefig("figures/dag.png", dpi=150); plt.close()
print("Saved figures/dag.png")
