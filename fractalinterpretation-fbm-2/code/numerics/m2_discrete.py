"""Exact E[Y_eps^2] of the discrete right-endpoint estimator (f = 1, t = 1, grid n):
Delta^2 sum_{i,j} (2 pi)^-1 det(Sigma(u_i,u_j) + eps^2 I)^-1/2, no Monte Carlo."""
import numpy as np, json, sys, os, time
n = 2 ** 15; D = 1.0 / n
u = np.arange(1, n + 1) / n
def m2d(H, eps, rows=256):
    p = u ** (2 * H); a = p + eps ** 2
    lagp = (np.arange(n) / n) ** (2 * H)            # |u_i - u_j|^{2H} by lag
    idx = np.arange(n)
    tot = 0.0
    for s in range(0, n, rows):
        i = idx[s:s + rows, None]
        c = 0.5 * (p[s:s + rows, None] + p[None, :] - lagp[np.abs(i - idx[None, :])])
        det = a[s:s + rows, None] * a[None, :] - c * c
        tot += np.sum(det ** -0.5)
    return tot * D * D / (2 * np.pi)
alphas = [float(x) for x in sys.argv[1:]]
fn = "results/m2_discrete_exact.json"
for al in alphas:
    t0 = time.time()
    r = {e: m2d(1 - al, e) for e in [0.01, 0.001]}
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    out[str(al)] = r; json.dump(out, open(fn, "w"))
    print(al, r, f"{time.time()-t0:.0f}s", flush=True)
