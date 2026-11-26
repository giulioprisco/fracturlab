import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
S = json.load(open("pathwise/summary.json"))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for a, r in S.items():
    d = r["one"]["diffs"]; e = np.array([x["eps"] for x in d])
    m = np.array([x["mc"] for x in d]); se = np.array([x["se"] for x in d])
    l = ax[0].errorbar(e, m / r["one"]["EY2"], yerr=1.96 * se / r["one"]["EY2"], marker="o", ms=3, capsize=2, label=f"$\\alpha={float(a):.1f}$")
    if "exact" in d[0]:
        ax[0].plot(e, np.array([x["exact"] for x in d]) / r["one"]["EY2"], ls=":", color=l[0].get_color())
    import os
    gf = f"pathwise/gridfree_{float(a):.2f}.json"
    if os.path.exists(gf) and float(a) <= 0.6:
        g = json.load(open(gf)); ax[0].plot(e, np.array([g[f"{x}"] for x in e]) / r["one"]["EY2"], ls="--", lw=0.9, color=l[0].get_color())
    c = np.array([x["ctrl"] for x in d]); ax[1].plot(e, c, marker="s", ms=3, label=f"$\\alpha={float(a):.1f}$")
for k in (0, 1):
    ax[k].set_xscale("log"); ax[k].set_yscale("log"); ax[k].set_xticks([0.002, 0.004, 0.008, 0.016, 0.032]); ax[k].set_xticklabels(["0.002", "0.004", "0.008", "0.016", "0.032"]); ax[k].minorticks_off(); ax[k].set_xlabel("$\\varepsilon$"); ax[k].grid(alpha=0.2)
ax[0].set_ylabel("$E(Y_\\varepsilon-Y_{\\varepsilon/2})^2\\,/\\,E Y_{0.001}^2$")
ax[0].set_title("fBm: MC (points), exact on grid (dotted), grid-free (dashed)", fontsize=10)
ax[1].set_ylabel("$E(Y_\\varepsilon-Y_{\\varepsilon/2})^2$"); ax[1].set_title("control $u^H\\xi$ (same marginals)", fontsize=10)
ax[0].legend(fontsize=7, ncol=2); plt.tight_layout(); plt.savefig("plots/pathwise.png", dpi=130)
