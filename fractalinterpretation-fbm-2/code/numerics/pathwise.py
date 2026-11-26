"""Pathwise test of Theorem item 3 (referee point 2), plus a fresh independent replicate of the
f = 1 mean test (referee point 5). Usage: python pathwise.py ALPHA N STREAM_INDEX
For each fBm path (exact Davies-Harte, n = 2^15, t = 1):
  * kernel estimators Y_eps on the window A = (2^-10, 1] for eps in EPS (f = 1 and f = sin t),
  * an independent crossing estimator C_k on the same window, at grid steps 2^k Delta (k = 0, 2, 4):
    each sign change of the sampled path in (u_a, u_b] gets weight int_{u_a}^{u_b} p_u(0) du / P(sign change),
    P(sign change) = arccos(rho)/pi exactly, so E C_k = E int_A g dL (up to g varying in a step),
  * Y_full at eps = 0.001 on all of [0,1] for f = 1 (replicate of Table 1, new independent stream).
Negative control: the same kernel quantities for V_u = u^H xi (one Gaussian xi per path), which has the same
one-point marginals as fBm (so the same E Y_eps) but no local time limit.
Streams: numpy SeedSequence(20261123).spawn(16)[STREAM_INDEX], independent by construction."""
import sys, time, os
import numpy as np
from fbmlt import fbm_paths, SQ2PI

a = float(sys.argv[1]); N = int(sys.argv[2]); si = int(sys.argv[3])
H = 1 - a; n = 2 ** 15; D = 1.0 / n; batch = 125
EPS = np.array([0.032, 0.016, 0.008, 0.004, 0.002, 0.001])
J0 = 2 ** 5                                # window starts at u = J0/n = 2^-10
KS = [0, 2, 4]
ss = np.random.SeedSequence(20261123).spawn(16)[si]
rng = np.random.default_rng(ss)
u = np.arange(0, n + 1) / n                # u_0 = 0 ... u_n = 1
fn = {"one": lambda s: np.ones_like(s), "sin": np.sin}
names = list(fn)
win = np.arange(J0 + 1, n + 1)             # right-endpoint indices with u_j in (2^-10, 1]
Gw = np.stack([fn[k](1 - u[win]) for k in names], 1) * D          # (nw, F)

def corr(x, y):
    cxy = 0.5 * (x ** (2 * H) + y ** (2 * H) - np.abs(x - y) ** (2 * H))
    return cxy / (x ** H * y ** H)
cross = {}
for k in KS:
    st = 2 ** k
    ia = np.arange(J0, n, st); ib = ia + st                   # intervals (u_a, u_b]
    P = np.arccos(np.clip(corr(u[ia], u[ib]), -1, 1)) / np.pi
    mass = (u[ib] ** (1 - H) - u[ia] ** (1 - H)) / ((1 - H) * SQ2PI)  # int p_u(0) du
    mid = 0.5 * (u[ia] + u[ib])
    W = np.stack([fn[f](1 - mid) for f in names], 1) * (mass / P)[:, None]
    cross[k] = (ia, ib, W)

Y = np.zeros((N, len(EPS), len(names))); C = np.zeros((N, len(KS), len(names)))
Yfull = np.zeros(N); Yv = np.zeros((N, len(EPS), len(names))); Yvfull = np.zeros(N)
nb = N // batch; t0 = time.time()
for b in range(nb):
    sl = slice(b * batch, (b + 1) * batch)
    B = fbm_paths(n, H, batch, rng)
    Bw = B[:, win].astype(np.float32)
    for i, e in enumerate(EPS):
        Wk = np.exp(-0.5 * (Bw / np.float32(e)) ** 2) * np.float32(1 / (e * SQ2PI))
        Y[sl, i, :] = Wk.astype(np.float64) @ Gw
    e = 0.001
    Yfull[sl] = (np.exp(-0.5 * (B[:, 1:] / e) ** 2) / (e * SQ2PI)).sum(1) * D
    for j, k in enumerate(KS):
        ia, ib, W = cross[k]
        sc = (B[:, ia] * B[:, ib] < 0).astype(np.float64)
        C[sl, j, :] = sc @ W
    # negative control: V_u = u^H xi
    xi = rng.standard_normal(batch)
    V = (u[None, :] ** H) * xi[:, None]
    Vw = V[:, win]
    for i, e in enumerate(EPS):
        Wk = np.exp(-0.5 * (Vw / e) ** 2) / (e * SQ2PI)
        Yv[sl, i, :] = Wk @ Gw
    Yvfull[sl] = (np.exp(-0.5 * (V[:, 1:] / 0.001) ** 2) / (0.001 * SQ2PI)).sum(1) * D
el = time.time() - t0
os.makedirs("pathwise", exist_ok=True)
np.savez(f"pathwise/pw_{a:.2f}.npz", a=a, H=H, n=n, N=nb * batch, EPS=EPS, KS=KS, J0=J0, names=names,
         Y=Y, C=C, Yfull=Yfull, Yv=Yv, Yvfull=Yvfull, runtime=el, stream=si)
print(f"alpha={a} N={nb*batch} runtime {el:.1f}s", flush=True)
