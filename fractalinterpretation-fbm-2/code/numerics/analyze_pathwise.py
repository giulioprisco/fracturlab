"""Summarise the pathwise test: pathwise/summary.json and pathwise/table.md."""
import numpy as np, json, glob, os
from fbmlt import discrete_target, SQ2PI
out = {}
def mse(x):  # mean of squares with SE
    m = np.mean(x); return m, np.std(x, ddof=1) / np.sqrt(len(x))
for fnm in sorted(glob.glob("pathwise/pw_*.npz")):
    d = np.load(fnm); a = float(d['a']); H = 1 - a; EPS = d['EPS']; N = int(d['N'])
    ex = json.load(open(f"pathwise/exact_{a:.2f}.json")) if os.path.exists(f"pathwise/exact_{a:.2f}.json") else None
    r = {"N": N, "EPS": EPS.tolist()}
    for fi, fname in enumerate(d['names']):
        Y = d['Y'][:, :, fi]; C = d['C'][:, :, fi]; Yv = d['Yv'][:, :, fi]
        diffs = [];  # E (Y_e - Y_{e/2})^2
        for i in range(len(EPS) - 1):
            m, se = mse((Y[:, i] - Y[:, i + 1]) ** 2)
            mv, sev = mse((Yv[:, i] - Yv[:, i + 1]) ** 2)
            row = {"eps": float(EPS[i]), "mc": m, "se": se, "ctrl": mv, "ctrl_se": sev,
                   "var": float(np.var(Y[:, i] - Y[:, i + 1], ddof=1))}
            if ex is not None and fname == "one":
                e1, e2 = EPS[i], EPS[i + 1]
                row["exact"] = ex[f"{e1}_{e1}"] + ex[f"{e2}_{e2}"] - 2 * ex[f"{e1}_{e2}"]
            diffs.append(row)
        m2 = float(np.mean(Y[:, -1] ** 2))
        cr = []
        for j, k in enumerate(d['KS']):
            m, se = mse((C[:, j] - Y[:, -1]) ** 2)
            cr.append({"k": int(k), "mse": m, "se": se, "rel": float(np.sqrt(m / m2)),
                       "corr": float(np.corrcoef(C[:, j], Y[:, -1])[0, 1]),
                       "meanC": float(C[:, j].mean()), "meanY": float(Y[:, -1].mean())})
        cy = [[float(np.mean((C[:, j] - Y[:, i]) ** 2)) for i in range(len(EPS))] for j in range(len(d['KS']))]
        r[str(fname)] = {"cy": cy, "diffs": diffs, "EY2": m2, "cross": cr,
                         "ctrl_mean": [float(x) for x in Yv.mean(0)], "Y_mean": [float(x) for x in Y.mean(0)],
                         "ctrl_EY2": [float(x) for x in (Yv ** 2).mean(0)], "Y_EY2": [float(x) for x in (Y ** 2).mean(0)],
                         "frac_ctrl_zero": float(np.mean(Yv[:, -1] < 1e-12))}
    # replicate of Table 1, f = 1, full [0,1], eps = 0.001
    tgt = discrete_target("one", 1.0, H, 0.001, int(d['n']))
    yf = d['Yfull']; m = yf.mean(); se = yf.std(ddof=1) / np.sqrt(N)
    yv = d['Yvfull']
    r["replicate"] = {"mean": float(m), "se": float(se), "target_discrete": float(tgt), "z": float((m - tgt) / se),
                      "ctrl_mean": float(yv.mean()), "ctrl_se": float(yv.std(ddof=1) / np.sqrt(N))}
    out[f"{a:.2f}"] = r
json.dump(out, open("pathwise/summary.json", "w"), indent=1)
for a, r in out.items():
    o = r["one"]
    print(f"alpha={a} N={r['N']}  E Y^2(eps=.001)={o['EY2']:.4f}")
    for row in o["diffs"]:
        print(f"  eps={row['eps']:.3f}  MC E(Ye-Ye/2)^2={row['mc']:.5f}+-{row['se']:.5f}  exact={row.get('exact',float('nan')):.5f}  z={(row['mc']-row.get('exact',np.nan))/row['se']:+.2f}  ctrl={row['ctrl']:.3g}+-{row['ctrl_se']:.2g}")
    for c in o["cross"]:
        print(f"  cross k={c['k']} relL2={c['rel']:.3f} corr={c['corr']:.4f} meanC={c['meanC']:.4f} meanY={c['meanY']:.4f}")
    print("  E(C_k - Y_eps)^2 rows k, cols eps:", [[round(x,4) for x in row] for row in o["cy"]])
    print("  ctrl means", np.round(o['ctrl_mean'],4), "Y means", np.round(o['Y_mean'],4), "ctrl zero frac", o['frac_ctrl_zero'])
    print("  ctrl EY2", np.round(o['ctrl_EY2'],3), "Y EY2", np.round(o['Y_EY2'],3))
    s = r["sin"]; print("  sin diffs", [f"{x['mc']:.2e}" for x in s["diffs"]], "cross rel", [round(c['rel'],3) for c in s['cross']])
    rp = r["replicate"]; print(f"  replicate f=1 eps=.001: {rp['mean']:.5f}+-{rp['se']:.5f} target {rp['target_discrete']:.5f} z={rp['z']:+.2f}; ctrl mean {rp['ctrl_mean']:.4f}+-{rp['ctrl_se']:.4f}")
