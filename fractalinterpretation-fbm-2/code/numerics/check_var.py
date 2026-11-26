"""Sanity check of the fBm sampler: mean of B_u^2 / u^{2H} and of phi_eps(B_u) vs exact, new stream."""
import numpy as np, sys
from fbmlt import fbm_paths, discrete_target, SQ2PI
a = float(sys.argv[1]); N = int(sys.argv[2]); H = 1 - a; n = 2 ** 15
rng = np.random.default_rng(np.random.SeedSequence(777).spawn(16)[int(sys.argv[3])])
u = np.arange(1, n + 1) / n; w = u ** (2 * H)
R = []; Yf = []
idx = [n // 1024, n // 32, n // 4, n]   # u = 2^-10, 2^-5, 1/4, 1
for b in range(N // 125):
    B = fbm_paths(n, H, 125, rng)[:, 1:]
    R.append((B[:, np.array(idx) - 1] ** 2) / w[np.array(idx) - 1])
    Yf.append((B ** 2 / w).mean(1))
    e = 0.001
R = np.concatenate(R); Yf = np.concatenate(Yf)
print("alpha", a, "N", len(R))
print(" E B_u^2/u^2H at u=2^-10,2^-5,1/4,1:", np.round(R.mean(0), 4), "+-", np.round(R.std(0) / np.sqrt(len(R)), 4))
print(" path-averaged B_u^2/u^2H:", round(Yf.mean(), 5), "+-", round(Yf.std() / np.sqrt(len(Yf)), 5))
