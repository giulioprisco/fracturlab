"""Grid-free exact E(Y_a - Y_b)^2, f = 1, window (2^-10, 1], b = a/2, by nested quadrature of the
difference kernel (2 pi)^-1 [d(a,a)^-1/2 + d(b,b)^-1/2 - d(a,b)^-1/2 - d(b,a)^-1/2], d(x,y) = det(Sigma + diag(x^2, y^2))."""
import numpy as np, json, sys
from scipy.integrate import quad
w0 = 2.0 ** -10
EPS = [0.032, 0.016, 0.008, 0.004, 0.002, 0.001]
def diff2(H, a, b):
    def k(u, v):
        p, q = u ** (2 * H), v ** (2 * H); c = 0.5 * (p + q - (u - v) ** (2 * H)); c2 = c * c
        d = lambda x, y: (p + x * x) * (q + y * y) - c2
        return d(a, a) ** -0.5 + d(b, b) ** -0.5 - d(a, b) ** -0.5 - d(b, a) ** -0.5
    sc = a ** (1 / H)
    def inner(u):
        pts = sorted(set(p for p in [u - sc * m for m in (1e-3, 1e-2, 0.1, 1, 10, 100, 1e3, 1e4)] if w0 < p < u))
        return quad(lambda v: k(u, v), w0, u, points=pts or None, limit=400, epsabs=1e-12, epsrel=1e-9)[0]
    val = quad(inner, w0, 1, limit=400, epsabs=1e-11, epsrel=1e-7)[0]
    return 2 * val / (2 * np.pi)
for al in [float(x) for x in sys.argv[1:]]:
    H = 1 - al
    r = {f"{EPS[i]}": diff2(H, EPS[i], EPS[i + 1]) for i in range(len(EPS) - 1)}
    json.dump(r, open(f"pathwise/gridfree_{al:.2f}.json", "w"), indent=1); print(al, r, flush=True)
