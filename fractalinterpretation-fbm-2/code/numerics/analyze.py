"""Analysis of results/alpha_*.npz: table with 95% CIs and plots."""
import glob, time, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from fbmlt import (FUNCS, LABELS, exact_target, smoothed_target, discrete_target,
                   second_moment_L, SQ2PI)
t0 = time.time()
R = {}
for fn in sorted(glob.glob("results/alpha_*.npz")):
    d = np.load(fn, allow_pickle=True); R[float(d["a"])] = d
alphas = sorted(R)
names = list(R[alphas[0]]["names"])
FN = ["lin", "sin", "exp", "sq"]
EPS = R[alphas[0]]["EPS"]; n = int(R[alphas[0]]["n"])
z = 1.96
ci = lambda y: z * y.std(ddof=1) / np.sqrt(len(y))
E1, E2 = 0.01, 0.001          # Richardson pair
i1, i2 = list(EPS).index(E1), list(EPS).index(E2)

rows = []
for a in alphas:
    d = R[a]; H = 1 - a; Y = d["Y"]; kap = (1 - H) / H
    for f in ["one"] + FN:
        k = names.index(f); ex = exact_target(f, 1.0, a)
        y = Y[:, i2, k]
        Dn = discrete_target(f, 1.0, H, E2, n); St = smoothed_target(f, 1.0, H, E2)
        row = dict(alpha=a, f=LABELS[f], exact=ex, est=y.mean(), ci=ci(y),
                   relerr=(y.mean() - ex) / ex * 100, relci=ci(y) / ex * 100,
                   eps_bias=(St - ex) / ex * 100, grid_bias=(Dn - St) / ex * 100,
                   z_vs_expected=(y.mean() - Dn) / (y.std(ddof=1) / np.sqrt(len(y))),
                   resid=(y.mean() - Dn) / ex * 100)
        if a <= 0.5:   # leading bias ~ eps^kappa (proof.md sec. 4); extrapolate per path
            w1, w2 = E1 ** kap, E2 ** kap
            yr = (Y[:, i2, k] * w1 - Y[:, i1, k] * w2) / (w1 - w2)
            row.update(rich=yr.mean(), rich_ci=ci(yr), rich_relerr=(yr.mean() - ex) / ex * 100,
                       rich_relci=ci(yr) / ex * 100)
        rows.append(row)

# second moment of L(0,[0,1]) (f = 1)
sm = {}
for a in alphas:
    H = 1 - a; k = names.index("one"); Y = R[a]["Y"][:, :, k]
    m2 = (Y ** 2).mean(axis=0); m2ci = z * (Y ** 2).std(axis=0, ddof=1) / np.sqrt(len(Y))
    cauchy = [((Y[:, i] - Y[:, i + 1]) ** 2).mean() for i in range(len(EPS) - 1)]
    sm[a] = dict(exact=second_moment_L(H), mc=m2.tolist(), mcci=m2ci.tolist(), cauchy=cauchy)

with open("results/table.md", "w") as fh:
    fh.write(f"# Results at t = 1, grid n = {n} (step {1/n:.2e}), eps = {E2}, N = {len(R[alphas[0]]['Y'])} paths per alpha\n\n")
    fh.write("Target: E int_0^1 f(1-u) L(0,du) = Gamma(alpha)/sqrt(2pi) * I^alpha f(1). 95% CIs.\n")
    fh.write("eps-bias, grid-bias: exact deterministic biases of the estimator (proof.md sec. 4), in % of target.\n")
    fh.write("z: (MC mean - exact mean of the discretised estimator)/SE; |z|<2 expected. Richardson: eps pair (0.01, 0.001), exponent alpha/(1-alpha), alpha<=0.5 only.\n\n")
    fh.write("| alpha | f | exact | estimate ± 95% CI | rel. err % (± CI) | eps-bias % | grid-bias % | MC − E[estimator] % (± CI) | z | Richardson rel. err % (± CI) |\n|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        rich = f"{r['rich_relerr']:+.2f} ± {r['rich_relci']:.2f}" if "rich" in r else "n/a"
        fh.write(f"| {r['alpha']} | {r['f']} | {r['exact']:.5f} | {r['est']:.5f} ± {r['ci']:.5f} | {r['relerr']:+.2f} ± {r['relci']:.2f} | "
                 f"{r['eps_bias']:+.3f} | {r['grid_bias']:+.3f} | {r['resid']:+.2f} ± {r['relci']:.2f} | {r['z_vs_expected']:+.2f} | {rich} |\n")
    m2s = json.load(open("results/m2_smoothed_exact.json"))
    m2d = json.load(open("results/m2_discrete_exact.json"))
    fh.write("\n## Second moment of local time, f = 1, t = 1\n\nexact E L^2: proof.md sec. 8. grid-free E Y_eps^2: 2-D quadrature (m2_exact.py). on-grid E Y_eps^2: exact double sum over the n = 2^15 grid (m2_discrete.py), the right reference for the MC. z = (MC - on-grid exact)/SE.\n\n| alpha | exact E[L(0,[0,1])^2] | grid-free E[Y_0.01^2] | on-grid E[Y_0.01^2] | MC E[Y_0.01^2] | z | grid-free E[Y_0.001^2] | on-grid E[Y_0.001^2] | MC E[Y_0.001^2] | z | MC E(Y_0.002 - Y_0.001)^2 |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
    for a in alphas:
        s = sm[a]; q = m2s[str(a)]; qd = m2d[str(a)]
        z1 = (s['mc'][i1] - qd['0.01']) / (s['mcci'][i1] / z); z2 = (s['mc'][i2] - qd['0.001']) / (s['mcci'][i2] / z)
        s.update(grid_exact={'0.01': qd['0.01'], '0.001': qd['0.001']}, z01=z1, z001=z2)
        fh.write(f"| {a} | {s['exact']:.4f} | {q['0.01']:.4f} | {qd['0.01']:.4f} | {s['mc'][i1]:.4f} ± {s['mcci'][i1]:.4f} | {z1:+.2f} | {q['0.001']:.4f} | {qd['0.001']:.4f} | {s['mc'][i2]:.4f} ± {s['mcci'][i2]:.4f} | {z2:+.2f} | {s['cauchy'][-1]:.2e} |\n")
json.dump(dict(rows=rows, second_moment=sm), open("results/summary.json", "w"), indent=1, default=float)

# ---------------- plots
cols = plt.cm.viridis(np.linspace(0, 0.9, len(alphas)))
# 1. curves over t
fig, ax = plt.subplots(len(alphas), 4, figsize=(15, 3 * len(alphas)), sharex=True)
fig2, ax2 = plt.subplots(2, 4, figsize=(16, 7.5))
u = np.arange(n + 1) / n; tt = u[1:]; sel = np.unique(np.geomspace(1, n, 300).astype(int)) - 1
for r_, a in enumerate(alphas):
    d = R[a]; wb = d["wbar_b"]; eps = float(d["EPS_MAIN"])
    for c_, f in enumerate(FN):
        F = FUNCS[f](u[:-1])                       # F[m] = f(m/n), m = 0..n-1
        curves = np.array([fftconvolve(F, w[1:])[:n] / n for w in wb])   # est(t_k), k = 1..n
        mu = curves.mean(0); band = z * curves.std(0, ddof=1) / np.sqrt(len(curves))
        ex = exact_target(f, tt, a)
        A = ax[r_, c_]
        A.plot(tt, ex, "b-", lw=2.5, label="exact Γ(α)/√(2π)·I^α f(t)")
        A.fill_between(tt[sel], (mu - band)[sel], (mu + band)[sel], color="r", alpha=0.3)
        A.plot(tt[sel], mu[sel], "r--", lw=1.2, label=f"MC, ε={eps}")
        A.set_title(f"α={a}, f={LABELS[f]}", fontsize=9)
        if r_ == 0 and c_ == 0: A.legend(fontsize=7)
        pe = np.concatenate([[0.0], 1 / np.sqrt(2 * np.pi * (u[1:] ** (2 * (1 - a)) + eps ** 2))])
        Ecur = fftconvolve(F, pe[1:])[:n] / n          # exact mean of the discretised estimator
        m = tt >= 0.01; sm_ = sel[tt[sel] >= 0.01]
        ax2[0, c_].plot(tt[sm_], ((mu - ex) / ex * 100)[sm_], color=cols[r_], label=f"α={a}")
        ax2[0, c_].plot(tt[sm_], ((Ecur - ex) / ex * 100)[sm_], color=cols[r_], ls=":", lw=1)
        ax2[0, c_].set_title(f"f = {LABELS[f]}: rel. error % vs t (ε={eps})\nsolid MC, dotted exact ε+grid bias", fontsize=9)
        ax2[1, c_].plot(tt[sm_], ((mu - Ecur) / ex * 100)[sm_], color=cols[r_])
        ax2[1, c_].fill_between(tt[sm_], ((-band) / ex * 100)[sm_], (band / ex * 100)[sm_], color=cols[r_], alpha=0.15)
        ax2[1, c_].set_title(f"f = {LABELS[f]}: (MC − E[estimator]) / exact, %, ±95% band", fontsize=9)
for A in ax[-1]: A.set_xlabel("t")
for A in ax2.ravel(): A.set_xscale("log"); A.axhline(0, color="k", lw=0.5); A.set_xlabel("t")
ax2[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig("plots/curves_vs_t.png", dpi=110)
fig2.tight_layout(); fig2.savefig("plots/relerr_vs_t.png", dpi=110)

# 2. error vs N
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for r_, a in enumerate(alphas):
    H = 1 - a; k = names.index("lin"); y = R[a]["Y"][:, i2, k]
    Ns = np.unique(np.geomspace(10, len(y), 60).astype(int))
    cm = np.cumsum(y)[Ns - 1] / Ns
    Dn = discrete_target("lin", 1.0, H, E2, n); ex = exact_target("lin", 1.0, a)
    ax[0].loglog(Ns, np.abs(cm - Dn) / ex, color=cols[r_], label=f"α={a}")
    ax[0].loglog(Ns, y.std() / np.sqrt(Ns) / ex, color=cols[r_], ls=":", lw=1)
    ax[1].loglog(Ns, np.abs(cm - ex) / ex, color=cols[r_], label=f"α={a}")
ax[0].set_title("f=1+2t, t=1: |MC mean − E[estimator]| / target (solid)\nσ/√N (dotted) — pure Monte Carlo error ∝ N^-1/2", fontsize=9)
ax[1].set_title(f"|MC mean − exact| / exact (ε={E2}): MC error until the ε-bias floor", fontsize=9)
for A in ax: A.set_xlabel("N (realizations)"); A.legend(fontsize=8); A.grid(alpha=0.3, which="both")
NN = np.array([10, 40000.]); ax[0].loglog(NN, 0.3 * NN ** -0.5, "k--", lw=1); ax[0].text(100, 0.3 * 100 ** -0.5 * 1.3, "N^-1/2")
fig.tight_layout(); fig.savefig("plots/error_vs_N.png", dpi=110)

# 3. error vs eps
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
eg = np.geomspace(1e-4, 0.1, 40)
for r_, a in enumerate(alphas):
    H = 1 - a; k = names.index("lin"); ex = exact_target("lin", 1.0, a); Y = R[a]["Y"][:, :, k]
    m = Y.mean(0); c = z * Y.std(0, ddof=1) / np.sqrt(len(Y))
    ax[0].errorbar(EPS, (m - ex) / ex * 100, yerr=c / ex * 100, fmt="o", color=cols[r_], ms=4, capsize=2, label=f"α={a} MC")
    th = [(smoothed_target("lin", 1.0, H, e) - ex) / ex * 100 for e in eg]
    ax[0].plot(eg, th, color=cols[r_], lw=1)
    ax[1].loglog(eg, np.abs(th), color=cols[r_], label=f"α={a}: |bias| (exact), slope ε^min(α/(1-α),2)")
    ax[1].loglog(EPS, np.abs(m - ex) / ex * 100, "o", color=cols[r_], ms=4)
ax[0].set_xscale("log"); ax[0].set_title("f=1+2t, t=1: relative error % vs kernel width ε\n(points: MC ± 95% CI; lines: exact ε-bias, proof.md §4)", fontsize=9)
ax[1].set_title("|relative error| %, log-log", fontsize=9)
for A in ax: A.set_xlabel("ε"); A.grid(alpha=0.3, which="both"); A.legend(fontsize=7)
fig.tight_layout(); fig.savefig("plots/error_vs_eps.png", dpi=110)

# 4. error vs grid step (eps = 0.01)
fig, ax = plt.subplots(figsize=(7, 5))
LEV = R[alphas[0]]["LEVELS"]; eg_ = float(R[alphas[0]]["EPS_GRID"])
for r_, a in enumerate(alphas):
    H = 1 - a; k = names.index("lin"); Yg = R[a]["Yg"][:, :, k]; St = smoothed_target("lin", 1.0, H, eg_)
    dts = np.array([2 ** l / n for l in LEV]); m = Yg.mean(0); c = z * Yg.std(0, ddof=1) / np.sqrt(len(Yg))
    ax.errorbar(dts, (m - St) / St * 100, yerr=c / St * 100, fmt="o", color=cols[r_], ms=4, capsize=2, label=f"α={a} MC")
    det = [(discrete_target("lin", 1.0, H, eg_, n // 2 ** l) - St) / St * 100 for l in LEV]
    ax.plot(dts, det, color=cols[r_], lw=1)
ax.set_xscale("log"); ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("grid step Δ"); ax.set_ylabel("% deviation from ε-smoothed target")
ax.set_title(f"Grid-step error at ε={eg_} (right-endpoint rule)\npoints: MC ± 95% CI, lines: exact discrete expectation", fontsize=9)
ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both"); fig.tight_layout(); fig.savefig("plots/error_vs_grid.png", dpi=110)

# 5. second moment
fig, ax = plt.subplots(figsize=(7, 5))
for r_, a in enumerate(alphas):
    s = sm[a]; ax.errorbar(EPS, np.array(s["mc"]) / s["exact"], yerr=np.array(s["mcci"]) / s["exact"], fmt="o-", color=cols[r_], ms=4, capsize=2, label=f"α={a} (exact E L²={s['exact']:.3f})")
ax.axhline(1, color="k", lw=0.8); ax.set_xscale("log"); ax.set_xlabel("ε"); ax.set_ylabel("MC E[Y_ε²] / exact E[L(0,[0,1])²]")
ax.set_title("Second moment of local time (tests the 2-point fBm structure)", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig("plots/second_moment.png", dpi=110)
print(f"analysis done in {time.time()-t0:.1f}s")
