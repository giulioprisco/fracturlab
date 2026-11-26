"""Analysis of the alpha = 0.3 refinement run (results/refine_p*.npz) vs the main run."""
import glob, json
import numpy as np
from fbmlt import FUNCS, LABELS, exact_target, smoothed_target, discrete_target
a = 0.3; H = 0.7; z = 1.96
Ys = [np.load(f, allow_pickle=True) for f in sorted(glob.glob("results/refine_p*.npz"))]
names = list(Ys[0]["names"]); EPS = list(Ys[0]["EPS"]); n = int(Ys[0]["n"])
Y = np.concatenate([d["Y"] for d in Ys]); N = len(Y)
runtime = [float(d["runtime"]) for d in Ys]
main = np.load("results/alpha_0.30.npz", allow_pickle=True); Ym = main["Y"]; EPSm = list(main["EPS"]); nm = int(main["n"])
rows = []
for f in ["one", "lin", "sin", "exp", "sq"]:
    ex = exact_target(f, 1.0, a); k = names.index(f); r = dict(f=LABELS[f], exact=ex)
    for lab, yy, nn, e in [("main", Ym[:, EPSm.index(0.001), names.index(f)], nm, 0.001),
                           ("fine_e1e-3", Y[:, EPS.index(0.001), k], n, 0.001),
                           ("fine", Y[:, EPS.index(0.00025), k], n, 0.00025)]:
        se = yy.std(ddof=1) / np.sqrt(len(yy)); St = smoothed_target(f, 1.0, H, e); Dn = discrete_target(f, 1.0, H, e, nn)
        r[lab] = dict(N=len(yy), n=nn, eps=e, raw=(yy.mean() - ex) / ex * 100, ci=z * se / ex * 100,
                      eps_bias=(St - ex) / ex * 100, grid_bias=(Dn - St) / ex * 100, pred=(Dn - ex) / ex * 100,
                      resid=(yy.mean() - Dn) / ex * 100, zs=(yy.mean() - Dn) / se)
    r["drop_obs"] = r["fine"]["raw"] - r["main"]["raw"]; r["drop_pred"] = r["fine"]["pred"] - r["main"]["pred"]
    r["drop_ci"] = np.hypot(r["fine"]["ci"], r["main"]["ci"])
    rows.append(r)
json.dump(dict(N=N, n=n, runtime=runtime, rows=rows), open("results/refine_summary.json", "w"), indent=1, default=float)
with open("results/refine_table.md", "w") as fh:
    fh.write(f"# alpha = 0.3 refinement: n = {n}, eps = 0.00025, N = {N} (vs main run n = {nm}, eps = 0.001, N = {len(Ym)})\n\n")
    fh.write("| f | main raw % | main predicted bias % | main z | fine raw % | fine eps-bias % | fine grid % | fine predicted % | fine residual % | fine z | observed drop % | predicted drop % | fine grid, eps=0.001: raw % / pred % / z |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        m, fi, g = r["main"], r["fine"], r["fine_e1e-3"]
        fh.write(f"| {r['f']} | {m['raw']:+.2f} ± {m['ci']:.2f} | {m['pred']:+.2f} | {m['zs']:+.2f} | {fi['raw']:+.2f} ± {fi['ci']:.2f} | {fi['eps_bias']:+.3f} | {fi['grid_bias']:+.3f} | {fi['pred']:+.2f} | {fi['resid']:+.2f} ± {fi['ci']:.2f} | {fi['zs']:+.2f} | {r['drop_obs']:+.2f} ± {r['drop_ci']:.2f} | {r['drop_pred']:+.2f} | {g['raw']:+.2f} / {g['pred']:+.2f} / {g['zs']:+.2f} |\n")
    fh.write(f"\nruntimes (s, 4 parallel processes): {runtime}\n")
print(open("results/refine_table.md").read())
