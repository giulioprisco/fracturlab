"""v0.5 independent-sampler cross-check of the one-sided z-scores (referee report on v0.4, point 4).
Exact fBm on the uniform grid u_j=j/n, j=1..n, n=2^12, by dense Cholesky factorisation of the fBm covariance
(no FFT / circulant embedding, no shared code with fbmlt.fbm_paths). Same test statistic as Table 1 and
check_mean.py: Y_eps = Delta sum_j phi_eps(B_{u_j}) (f=1, t=1, right-end-point rule), compared with its exact
discrete mean Delta sum_j (2pi(u_j^{2H}+eps^2))^{-1/2}; plus low-variance checks Delta sum cos(th B_u).
Usage: python cholesky_check.py ALPHA N_TOTAL BATCH [STREAM] [LOG2N]  (STREAM default 0 = the v0.5 alpha=0.8 run; LOG2N default 12)"""
import numpy as np, sys, json, time
a = float(sys.argv[1]); NT = int(sys.argv[2]); BS = int(sys.argv[3])
si = int(sys.argv[4]) if len(sys.argv) > 4 else 0; lg = int(sys.argv[5]) if len(sys.argv) > 5 else 12
H = 1 - a; n = 2 ** lg; D = 1.0 / n
u = np.arange(1, n + 1) / n
C = 0.5 * (u[:, None] ** (2 * H) + u[None, :] ** (2 * H) - np.abs(u[:, None] - u[None, :]) ** (2 * H))
L = np.linalg.cholesky(C)  # exact, no jitter needed (u_j > 0)
del C
rng = np.random.default_rng(np.random.SeedSequence(20261125).spawn(8)[si])  # new entropy, unrelated to earlier runs
EP = [0.01, 0.001]; TH = [1.0, 4.0]
cols = []; t0 = time.time()
for b in range(NT // BS):
    B = L @ rng.standard_normal((n, BS))          # n x BS, column = path
    row = [(np.exp(-0.5 * (B / e) ** 2) / (e * np.sqrt(2 * np.pi))).sum(0) * D for e in EP]
    row += [np.cos(th * B).sum(0) * D for th in TH]
    cols.append(np.stack(row, 1))
    if b % 20 == 0: print(b * BS, round(time.time() - t0, 1), flush=True)
V = np.concatenate(cols); N = len(V); m = V.mean(0); se = V.std(0, ddof=1) / np.sqrt(N)
ex = [np.sum(1 / np.sqrt(2 * np.pi * (u ** (2 * H) + e ** 2))) * D for e in EP] + [np.sum(np.exp(-0.5 * th ** 2 * u ** (2 * H))) * D for th in TH]
names = ["Y0.01", "Y0.001", "cos1", "cos4"]
out = {"alpha": a, "H": H, "n": n, "N": N, "stream": si, "seconds": time.time() - t0,
       **{k: {"mc": float(m[i]), "se": float(se[i]), "exact": float(ex[i]), "z": float((m[i] - ex[i]) / se[i]),
              "rel_resid_pct": float((m[i] - ex[i]) / ex[i] * 100), "rel_se_pct": float(se[i] / ex[i] * 100)} for i, k in enumerate(names)}}
json.dump(out, open(f"pathwise/cholesky_{a:.2f}.json", "w"), indent=1)
print(json.dumps(out, indent=1))
