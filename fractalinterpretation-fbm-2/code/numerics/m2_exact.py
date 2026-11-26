"""Exact E[Y_eps^2] for f = 1, t = 1 (no grid), by 2-D quadrature of (2 pi)^-1 det(Sigma + eps^2 I)^-1/2."""
import numpy as np, sys, json
from scipy.integrate import quad
def m2(H, eps):
    def inner(u):
        def g(v):
            a = u**(2*H) + eps**2; b = v**(2*H) + eps**2
            c = 0.5*(u**(2*H) + v**(2*H) - (u - v)**(2*H))
            return (a*b - c*c) ** -0.5
        pts = sorted(set([p for p in [eps**(1/H), u - eps**(1/H), u/2] if 0 < p < u]))
        return quad(g, 0, u, points=pts, limit=200, epsrel=1e-8)[0]
    s = eps**(1/H)
    pts = [p for p in s*np.logspace(0, 6, 13) if p < 1]
    return 2/(2*np.pi) * quad(inner, 0, 1, points=pts, limit=200, epsrel=1e-7)[0]
out = {}
import os
prev = json.load(open("results/m2_smoothed_exact.json")) if os.path.exists("results/m2_smoothed_exact.json") else {}
for a in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9]:
    if str(a) in prev:
        out[a] = {float(k): v for k, v in prev[str(a)].items()}; continue
    out[a] = {e: m2(1-a, e) for e in [0.01, 0.001]}
    print(a, out[a], flush=True)
json.dump({str(k): v for k, v in out.items()}, open("results/m2_smoothed_exact.json", "w"))
