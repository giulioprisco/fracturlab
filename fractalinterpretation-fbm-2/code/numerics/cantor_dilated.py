"""Randomly dilated middle-thirds Cantor measure: M(A) = 2^{-U} mu_inf(3^U A), U ~ Unif[0,1).
Check that G(u) = E M([0,u]) = G(1) u^beta, beta = log2/log3, and compute G(1).
Also the undilated Cantor function F(u)/u^beta (log-periodic)."""
import numpy as np
from scipy.integrate import quad
beta = np.log(2) / np.log(3)
def F(x, depth=40):  # Cantor function on [0,1]
    y = 0.0; s = 0.5
    for _ in range(depth):
        x *= 3
        if x < 1: pass
        elif x < 2: return y + s
        else: y += s; x -= 2
        s /= 2
    return y
def Finf(x):
    k = 0
    while x > 1: x /= 3; k += 1
    return 2 ** k * F(x)
def G(u):
    # integrand has jumps in slope only; split finely
    pts = np.linspace(0, 1, 401)
    return sum(quad(lambda v: 2 ** (-v) * Finf(3 ** v * u), pts[i], pts[i + 1], limit=50)[0] for i in range(400))
for u in [0.01, 0.1, 0.37, 1.0, 2.5]:
    print(f"u={u}: G(u)/u^beta = {G(u)/u**beta:.6f}   undilated F_inf(u)/u^beta = {Finf(u)/u**beta:.4f}")
# Monte Carlo of E[(2+X)^{-beta}]/log(2), X ~ Cantor distribution
# (standard ternary expansion with digits 0 and 2). Matches Example 4.
rng = np.random.default_rng(0)
Nmc = 400_000
bits = rng.integers(0, 2, size=(Nmc, 40))
powers = 3.0 ** -np.arange(1, 41)
X = (2 * bits * powers).sum(axis=1)
vals = (2 + X) ** (-beta) / np.log(2)
est = vals.mean()
se = vals.std(ddof=1) / np.sqrt(Nmc)
print(f"MC E[(2+X)^{{-beta}}]/log(2) = {est:.5f} +/- {se:.5f}")
