"""Refinement run: n = 2^LOG2N grid, eps in EPS. Usage: python run_refine.py ALPHA N LOG2N SEED TAG
Saves results/refine_<tag>.npz with per-path Y (N, len(EPS), F)."""
import sys, time, os
import numpy as np
from fbmlt import fbm_paths, FUNCS, SQ2PI
a = float(sys.argv[1]); N = int(sys.argv[2]); L2N = int(sys.argv[3]); seed = int(sys.argv[4]); tag = sys.argv[5]
H = 1 - a; n = 2 ** L2N; batch = 100
EPS = [0.00025, 0.001]
names = list(FUNCS)
rng = np.random.default_rng(seed)
u = np.arange(1, n + 1) / n
G = np.stack([FUNCS[k](1 - u) for k in names], axis=1) / n
Y = np.zeros((N, len(EPS), len(names)))
t0 = time.time(); tgen = 0.0
for b in range(N // batch):
    t1 = time.time(); B = fbm_paths(n, H, batch, rng); tgen += time.time() - t1
    Bi = B[:, 1:].astype(np.float32); del B
    sl = slice(b * batch, (b + 1) * batch)
    for i, e in enumerate(EPS):
        W = np.exp(-0.5 * (Bi / np.float32(e)) ** 2) * np.float32(1 / (e * SQ2PI))
        Y[sl, i, :] = W.astype(np.float64) @ G
el = time.time() - t0
os.makedirs("results", exist_ok=True)
np.savez(f"results/refine_{tag}.npz", a=a, H=H, n=n, N=N, EPS=EPS, names=names, Y=Y, runtime=el, tgen=tgen, seed=seed)
print(f"tag={tag} alpha={a} N={N} n={n} runtime {el:.1f}s (fBm {tgen:.1f}s)", flush=True)
