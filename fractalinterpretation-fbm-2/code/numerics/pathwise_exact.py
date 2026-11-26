"""Exact on-grid second moments for the pathwise test (f = 1, window (2^-10, 1], n = 2^15):
E[Y_a Y_b] = Delta^2/(2 pi) sum_{i,j in window} det(Sigma(u_i,u_j) + diag(a^2, b^2))^{-1/2}
for a = b in EPS and (a, b) = (eps, eps/2). Hence the exact E(Y_eps - Y_eps/2)^2. Usage: ALPHA..."""
import numpy as np, json, sys, os, time
n = 2 ** 15; D = 1.0 / n; J0 = 2 ** 5
EPS = [0.032, 0.016, 0.008, 0.004, 0.002, 0.001]
idx = np.arange(J0 + 1, n + 1); u = idx / n
def run(H, rows=128):
    p = u ** (2 * H); lagp = (np.arange(n + 1) / n) ** (2 * H)
    pairs = [(e, e) for e in EPS] + [(EPS[i], EPS[i + 1]) for i in range(len(EPS) - 1)]
    tot = np.zeros(len(pairs))
    for s in range(0, len(u), rows):
        ii = idx[s:s + rows, None]
        c2 = (0.5 * (p[s:s + rows, None] + p[None, :] - lagp[np.abs(ii - idx[None, :])])) ** 2
        for q, (a, b) in enumerate(pairs):
            # symmetrise the cross term: average of diag(a^2,b^2) and diag(b^2,a^2)
            d1 = (p[s:s + rows, None] + a * a) * (p[None, :] + b * b) - c2
            if a == b:
                tot[q] += np.sum(d1 ** -0.5)
            else:
                tot[q] += np.sum(d1 ** -0.5)
    tot *= D * D / (2 * np.pi)
    return {f"{a}_{b}": float(v) for (a, b), v in zip(pairs, tot)}
pass
for al in [float(x) for x in sys.argv[1:]]:
    t0 = time.time(); r = run(1 - al)
    json.dump(r, open(f"pathwise/exact_{al:.2f}.json", "w"), indent=1)
    print(al, f"{time.time()-t0:.0f}s", flush=True)
