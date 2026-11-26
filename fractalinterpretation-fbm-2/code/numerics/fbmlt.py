"""Core routines: exact fBm simulation (Davies-Harte / circulant embedding),
Gaussian-kernel local-time functionals, and exact Riemann-Liouville integrals.

Conventions: B = fBm with Hurst H on [0,1], B_0 = 0, Var B_u = u^{2H}.
alpha = 1 - H. Claim tested:  E int_0^t f(t-u) L(0,du) = Gamma(alpha)/sqrt(2 pi) * I^alpha f(t).
"""
import numpy as np
from scipy.special import gamma, gammainc
from scipy.integrate import quad

SQ2PI = np.sqrt(2 * np.pi)

# ---------------------------------------------------------------- fBm
_eig_cache = {}

def fgn_eigs(n, H):
    """Eigenvalues of the 2n circulant embedding of the fGn autocovariance."""
    key = (n, H)
    if key not in _eig_cache:
        k = np.arange(n + 1, dtype=float)
        g = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
        c = np.concatenate([g, g[-2:0:-1]])           # length 2n
        lam = np.fft.fft(c).real
        if lam.min() < -1e-10 * lam.max():
            raise ValueError("circulant embedding not nonnegative definite")
        _eig_cache[key] = np.sqrt(np.clip(lam, 0, None) / (2 * n))
    return _eig_cache[key]

def fbm_paths(n, H, m, rng):
    """m exact fBm paths on the grid u_j = j/n, j = 0..n; returns array (m, n+1).
    Each complex FFT yields two independent paths (real and imaginary parts)."""
    s = fgn_eigs(n, H)
    mc = (m + 1) // 2
    Z = rng.standard_normal((mc, 2 * n)) + 1j * rng.standard_normal((mc, 2 * n))
    Y = np.fft.fft(Z * s, axis=1)[:, :n]
    G = np.concatenate([Y.real, Y.imag], axis=0)[:m]  # fGn, unit step
    B = np.empty((m, n + 1))
    B[:, 0] = 0.0
    np.cumsum(G, axis=1, out=B[:, 1:])
    B *= n ** (-H)
    return B

def phi_eps(x, eps):
    return np.exp(-0.5 * (x / eps) ** 2) / (eps * SQ2PI)

# ---------------------------------------------------------------- test functions and exact I^alpha
def I_alpha(name, t, a):
    """Exact left RL integral I^a f(t) (lower limit 0)."""
    t = np.asarray(t, float)
    if name == "one":
        return t ** a / gamma(a + 1)
    if name == "lin":        # f = 1 + 2 s
        return t ** a / gamma(a + 1) + 2 * t ** (a + 1) / gamma(a + 2)
    if name == "sq":         # f = s^2
        return 2 * t ** (a + 2) / gamma(a + 3)
    if name == "exp":        # f = e^s  ->  e^t * P(a, t)
        return np.exp(t) * gammainc(a, t)
    if name == "sin":        # f = sin s -> sum_k (-1)^k t^{2k+1+a} / Gamma(2k+2+a)
        out = np.zeros_like(t)
        for k in range(40):
            out += (-1) ** k * t ** (2 * k + 1 + a) / gamma(2 * k + 2 + a)
        return out
    raise KeyError(name)

FUNCS = {"one": lambda s: np.ones_like(np.asarray(s, float)), "lin": lambda s: 1 + 2 * s, "sq": lambda s: s ** 2, "exp": np.exp, "sin": np.sin}
LABELS = {"one": "1", "lin": "1+2t", "sin": "sin t", "exp": "exp t", "sq": "t^2"}

def I_alpha_quad(name, t, a):
    """Independent check of I_alpha by adaptive quadrature with algebraic weight."""
    f = FUNCS[name]
    v, _ = quad(lambda s: f(s), 0, t, weight="alg", wvar=(0, a - 1), limit=200)
    return v / gamma(a)

def exact_target(name, t, a):
    """K * I^a f(t), K = Gamma(a)/sqrt(2 pi): the claimed E int f dL."""
    return gamma(a) / SQ2PI * I_alpha(name, t, a)

def smoothed_target(name, t, H, eps):
    """E int_0^t f(t-u) phi_eps(B_u) du = int_0^t f(t-u) / sqrt(2 pi (u^{2H}+eps^2)) du (exact)."""
    f = FUNCS[name]
    g = lambda u: f(t - u) / np.sqrt(2 * np.pi * (u ** (2 * H) + eps ** 2))
    # split on a geometric grid around the kernel scale eps^{1/H}
    sc = eps ** (1 / H)
    pts = [0.0] + [p for p in sc * np.logspace(-3, 6, 28) if p < t] + [t]
    return sum(quad(g, pts[i], pts[i + 1], limit=200, epsabs=0, epsrel=1e-12)[0]
               for i in range(len(pts) - 1))

def discrete_target(name, t, H, eps, n):
    """Exact expectation of the right-endpoint Riemann-sum estimator on grid step 1/n."""
    u = np.arange(1, int(round(t * n)) + 1) / n
    f = FUNCS[name]
    return np.sum(f(t - u) / np.sqrt(2 * np.pi * (u ** (2 * H) + eps ** 2))) / n

def second_moment_L(H):
    """Exact E[L(0,[0,1])^2] = J / (pi (2-2H)),  J = int_0^1 D(x)^{-1/2} dx,
    D(x) = x^{2H} - ((x^{2H} + 1 - (1-x)^{2H})/2)^2  (2x2 Gaussian determinant, scaled)."""
    def D(x):  # cancellation-free form of x^{2H} - c^2, c = (x^{2H}+1-(1-x)^{2H})/2
        y = 1.0 - x
        c = 0.5 * (x ** (2 * H) + 1 - y ** (2 * H))
        one_minus_xH = -np.expm1(H * np.log1p(-y))
        return (x ** H + c) * 0.5 * (y ** (2 * H) - one_minus_xH ** 2)
    J = sum(quad(lambda x: D(x) ** -0.5, lo, hi, limit=400, epsrel=1e-11)[0]
            for lo, hi in [(0, 1e-6), (1e-6, 0.5), (0.5, 1 - 1e-6), (1 - 1e-6, 1)])
    return J / (np.pi * (2 - 2 * H))
