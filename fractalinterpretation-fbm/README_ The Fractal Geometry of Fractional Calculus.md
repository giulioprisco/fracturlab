---
title: 'README: The Fractal Geometry of Fractional Calculus'

---

# README: The Fractal Geometry of Fractional Calculus (Gemini)

This document describes a bridge between **Fractal Geometry** and **Fractional Calculus**. By sampling a linear function only at the points belonging to the **fractal zero-set** of a stochastic process, the deterministic Riemann-Liouville integral is recovered. The Python script shows this numerically for the linear case. This document and the script have been developed with AI assistance from Gemini.

---

## 1. Mathematical Foundation

The Riemann-Liouville fractional integral of a function $f(t)$ can be represented as the **Expectation Value ($\mathbb{E}$)** of an integral weighted by the "Local Time" ($\ell_t^0$)—the density of zeroes—of a Fractional Brownian Motion (fBm). 

The use of the expectation value is crucial; while any single fractal path is a noisy, unpredictable realization, the "long-run average" across a high number of realizations converges to a perfectly smooth, deterministic calculus result. In our simulation, the average of 1,000 realizations serves as the numerical approximation of this expected value.



For an fBm process $B_t^H$ with a Hurst exponent $H = 1 - \alpha$, the Riemann-Liouville fractional integral $J^\alpha f(x)$ is:
$$J^\alpha f(x) = K \cdot \mathbb{E} \left[ \int_0^x f(t) \delta(B_t^H) dt \right]$$

### Scaling ($K$ and $C$)
The proportionality factor **$K$** acts as the bridge between raw stochastic sampling and analytical calculus:
* **Calculus Normalization:** $1/\Gamma(\alpha)$, ensuring the operator adheres to the standard definition.
* **Stochastic Density Constant ($C$):** For an fBm process, the expected density of zeroes scales as $t^{H-1}$. By choosing the mapping $H = 1 - \alpha$, we ensure the density scales as $t^{-\alpha}$, which matches the precise power-law kernel $(x-t)^{\alpha-1}$ required by the integral. 

The constant $C$ depends on the standard deviation $\sigma$ of the process ($C = 1/(\sqrt{2\pi}\sigma)$). In our "unit case" ($\sigma=1.0$), $C$ is a universal statistical property of the fractal. If $\sigma$ were changed, $K$ would adjust proportionally to compensate, but the **shape** of the result (the calculus) would remain invariant.

---

## 2. Visual Interpretation: The Geometry of the Mist

The visualization collapses thousands of random universes into a single deterministic result through three distinct layers:

* **The Mist (Fractal Realizations):** Each faint path (rendered at **2%** opacity) represents $f(t)$ "seen" only when the path $B_t^H$ is near zero. For persistent paths ($H > 0.5$), the mist appears "clumpy" and concentrated near the origin, reflecting the long-range memory of the process.
* **The Blue Line (Analytical Ground Truth):** This represents the "deterministic target" calculated using the Gamma function. Its vertical "kick" at the start represents the infinite memory of the origin inherent in fractional kernels.
* **The Red Line (Stochastic Average):** The empirical result of our simulation. According to the **Law of Large Numbers**, as the number of realizations increases, the red line converges to the blue line.

---

## 3. Generalization: The "Fractal Comb"

While the current simulation uses a linear function, this method generalizes to **nonlinear functions** (e.g., $f(t) = \sin(t)$). Intuitively, this works because any smooth function is locally linear; formally, it works because **Strong Local Nondeterminism** ensures the increments of the process are asymptotically independent. This allows the fractal zero-set to act as a **Universal Stencil** or **Fractal Comb**.



Imagine a comb where the "teeth" (points where the process hits zero) are distributed according to the fractal geometry of the path. Because this comb is constructed independently of the function $f(t)$, it will "harvest" values from any curve you slide under it at the precise rates required by the Riemann-Liouville identity. The memory is a property of the **Comb's Geometry**, not the function being combed.

---

## 4. Numeric Challenges

Translating infinite fractal properties into a **Macbook Pro with M1 processor and 32GB memory** revealed the physical boundaries of silicon:

* **The $O(n^2)$ Computational Wall:** While $\alpha = 0.5$ (Standard Brownian Motion) runs in linear time ($O(n)$), any deviation forces the switch to full covariance matrix calculations ($O(n^2)$), pushing runtimes from seconds to ~15 minutes.
* **The "S-Curve" Lag:** For small $\alpha$, the theoretical density at $t=0$ is infinite. On a discrete grid, the simulation "starves" for data points at the very start because persistent paths wander away from zero. This creates a characteristic S-shaped lag in the red line.
* **The "Jitter Tax":** A narrower sensing window ($\epsilon$) improves mathematical accuracy but introduces high statistical noise, requiring more memory and CPU cycles to smooth the result.

---

## 5. Implementation Details

### Approximating the Fractal Set
In a discrete simulation, a path almost never hits exactly $0.0000$. We approximate the **Dirac delta ($\delta$)** using a **Gaussian Kernel Density Estimate (KDE)**. By giving the delta a small width ($\epsilon$), we create a "glow" around zero that allows us to capture the **fractal set of zeroes** numerically. Each realization is weighted by $w = \exp(-B_t^2 / 2\epsilon^2)$.

### Required Libraries
The Python script requires:
* `numpy`
* `matplotlib`
* `scipy`
* `fbm`

---

## 6. References

1. **Chen, W., Sun, H., & Li, X. (2022).** *Fractional Derivative Modeling in Mechanics and Engineering.*
2. **Nigmatullin, R. R., Zhang, W., & Gubaidullin, I. (2017).** *Accurate Relationship between Fractals and Fractional Integrals: New Approaches and Evaluations.*
3. **Mandelbrot, B. B. (1982).** *The Fractal Geometry of Nature.*
4. **Mandelbrot, B. B., & Van Ness, J. W. (1968).** *Fractional Brownian Motions, Fractional Noises and Applications.*
5. **Xiao, Y. (2007).** *Strong local nondeterminism and sample path properties of Gaussian random fields.* In: ***Asymptotic Theory in Probability and Statistics with Applications*** (T. L. Lai, Q. Shao, L. Wang, Eds.), pp. 136–176. World Scientific, Singapore.
