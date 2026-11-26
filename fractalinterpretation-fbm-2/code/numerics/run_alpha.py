"""Monte Carlo run for one alpha. Usage: python run_alpha.py ALPHA [N] [LOG2N] [SEED]
Saves results/alpha_<a>.npz with per-path functionals."""
import sys, time, os
import numpy as np
from fbmlt import fbm_paths, FUNCS, SQ2PI

a = float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 16000
L2N = int(sys.argv[3]) if len(sys.argv) > 3 else 15
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 12345
H = 1 - a; n = 2 ** L2N; batch = 250
EPS = np.array([0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001])
EPS_MAIN = 0.002                        # used for the t-curves
LEVELS = [0, 1, 2, 3, 4, 5, 6]          # grid coarsening 2^k (grid scan at eps = 0.01)
EPS_GRID = 0.01
names = list(FUNCS)
rng = np.random.default_rng(seed + int(round(a * 1000)))
u = np.arange(1, n + 1) / n             # right-endpoint rule, t = 1
G = np.stack([FUNCS[k](1 - u) for k in names], axis=1) / n          # (n, F)
Gl = {k: np.stack([FUNCS[f](1 - u[2**k - 1::2**k]) for f in names], axis=1) * (2**k / n) for k in LEVELS}
Y = np.zeros((N, len(EPS), len(names)))
Yg = np.zeros((N, len(LEVELS), len(names)))
nb = N // batch
wbar_b = np.zeros((nb, n + 1))          # batch means of phi_eps_main(B_u) for curves
t0 = time.time(); tgen = 0.0
for b in range(nb):
    t1 = time.time(); B = fbm_paths(n, H, batch, rng); tgen += time.time() - t1
    Bi = B[:, 1:].astype(np.float32)
    sl = slice(b * batch, (b + 1) * batch)
    for i, e in enumerate(EPS):
        W = np.exp(-0.5 * (Bi / np.float32(e)) ** 2) * np.float32(1 / (e * SQ2PI))
        Y[sl, i, :] = W.astype(np.float64) @ G
        if e == EPS_MAIN:
            wbar_b[b, 1:] = W.mean(axis=0)
            wbar_b[b, 0] = 1 / (e * SQ2PI)
        if e == EPS_GRID:
            for j, k in enumerate(LEVELS):
                Yg[sl, j, :] = W[:, 2**k - 1::2**k].astype(np.float64) @ Gl[k]
el = time.time() - t0
os.makedirs("results", exist_ok=True)
np.savez(f"results/alpha_{a:.2f}.npz", a=a, H=H, n=n, N=nb * batch, EPS=EPS, EPS_MAIN=EPS_MAIN,
         LEVELS=LEVELS, EPS_GRID=EPS_GRID, names=names, Y=Y[:nb*batch], Yg=Yg[:nb*batch],
         wbar_b=wbar_b, runtime=el, tgen=tgen)
print(f"alpha={a} N={nb*batch} n={n} runtime {el:.1f}s (fBm generation {tgen:.1f}s)")
