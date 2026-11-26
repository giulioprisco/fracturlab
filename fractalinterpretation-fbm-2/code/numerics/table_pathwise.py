import json, numpy as np
S = json.load(open("pathwise/summary.json"))
def f3(x):  # in units of 1e-3
    return f"{x*1e3:.3g}" if x*1e3 >= 0.01 else f"{x*1e3:.1e}".replace("e-0", "\\text{e-}")
rows = []
for a, r in S.items():
    o = r["one"]; d = {round(x["eps"], 3): x for x in o["diffs"]}
    g = json.load(open(f"pathwise/gridfree_{float(a):.2f}.json"))
    cells = [f"{float(a):.1f}"]
    for e in (0.032, 0.008, 0.002):
        x = d[e]
        cells += [f"${x['mc']*1e3:.3g}\\pm{1.96*x['se']*1e3:.2g}$", f"{x['exact']*1e3:.3g}"]
        gv = g[str(e)]
        cells += [f"{gv*1e3:.3g}" if float(a) <= 0.6 else "$<10^{-2}$" if e == 0.032 else "$\\ll10^{-2}$"]
    cells += [f"{d[0.032]['ctrl']:.2f}", f"{d[0.002]['ctrl']:.1f}"]
    rows.append(" & ".join(cells) + "\\\\")
print("\n".join(rows))
print()
rows = []
for a, r in S.items():
    o = r["one"]; cy = np.array(o["cy"]); EY2 = o["EY2"]; EPS = r["EPS"]
    best = []
    for j, k in enumerate((0, 2, 4)):
        i = int(np.argmin(cy[j])); best.append((np.sqrt(cy[j, i] / EY2), EPS[i]))
    rp = r["replicate"]
    rows.append(f"{float(a):.1f} & " + " & ".join(f"{b*100:.1f} ({e:g})" for b, e in best) +
                f" & {o['cross'][0]['meanC']:.4f} & {o['cross'][0]['meanY']:.4f} & ${rp['mean']:.5f}\\pm{1.96*rp['se']:.5f}$ & {rp['target_discrete']:.5f} & ${rp['z']:+.2f}$\\\\")
print("\n".join(rows))
zs = [r["replicate"]["z"] for r in S.values()]; print("sum z", sum(zs), "comb", sum(zs)/np.sqrt(7), "ssq", sum(z*z for z in zs))
