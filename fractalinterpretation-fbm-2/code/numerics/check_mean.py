"""Independent replicate of the mean test with extra low-variance marginal checks.
Usage: python check_mean.py ALPHA N STREAM. Functionals per path (u_j = j/n, j=1..n, Delta = 1/n):
F_th = Delta sum cos(th B_u) (exact Delta sum exp(-th^2 u^{2H}/2)), th = 1, 4, 16; Y_eps for eps = 0.01, 0.001."""
import numpy as np, sys, json
from fbmlt import fbm_paths, discrete_target, SQ2PI
a = float(sys.argv[1]); N = int(sys.argv[2]); si = int(sys.argv[3]); H = 1 - a; n = 2 ** 15; D = 1 / n
rng = np.random.default_rng(np.random.SeedSequence(99991).spawn(16)[si])
u = np.arange(1, n + 1) / n; TH = [1.0, 4.0, 16.0]; EP = [0.01, 0.001]
vals = []
for b in range(N // 125):
    B = fbm_paths(n, H, 125, rng)[:, 1:]
    cols = [np.cos(th * B).sum(1) * D for th in TH] + [(np.exp(-0.5 * (B / e) ** 2) / (e * SQ2PI)).sum(1) * D for e in EP]
    vals.append(np.stack(cols, 1))
V = np.concatenate(vals); m = V.mean(0); se = V.std(0, ddof=1) / np.sqrt(len(V))
ex = [np.sum(np.exp(-0.5 * th ** 2 * u ** (2 * H))) * D for th in TH] + [discrete_target("one", 1.0, H, e, n) for e in EP]
names = ["cos1", "cos4", "cos16", "Y0.01", "Y0.001"]
out = {k: {"mc": float(m[i]), "se": float(se[i]), "exact": float(ex[i]), "z": float((m[i] - ex[i]) / se[i])} for i, k in enumerate(names)}
json.dump(out, open(f"pathwise/checkmean_{a:.2f}.json", "w"), indent=1)
print(a, len(V), {k: (round(v["mc"], 5), round(v["exact"], 5), round(v["z"], 2)) for k, v in out.items()}, flush=True)
