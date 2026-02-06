---
title: readme

---

# readme

A MATLAB implementation for generating high-resolution 3D fractal surfaces using **Fractional Brownian Motion (fBm)**. This project utilizes the Spectral Synthesis (Fourier-based) method to create realistic, self-similar terrains.

## 🚀 Features
- **High-Resolution Support:** Optimized for large grids (tested up to $4096 \times 4096$).
- **Spectral Synthesis:** Efficient terrain generation via the Fast Fourier Transform (FFT).
- **Adjustable Roughness:** Fine-tune surface texture using the **Hurst exponent ($H$)**.
- **Procedural Visualization:** Includes a custom topographic colormap (Sea-to-Snow) and realistic lighting setup.

## 🛠 How It Works
The algorithm constructs a surface by summing thousands of sine waves in the frequency domain. The power spectrum follows a specific power law defined by the Hurst exponent:

$$P(f) \propto \frac{1}{f^{2H+2}}$$

- **$H \approx 0.2$:** Produces "anticorrelated" noise, resulting in jagged, prickly, and highly fragmented surfaces.
- **$H \approx 0.8$:** Produces "persistent" noise, resulting in smooth, rolling hills and broad geological features.

## 📋 Prerequisites
- MATLAB (Compatible with R2021a and later).
- No additional toolboxes are strictly required (uses standard built-in functions).

## 📖 Usage
1. Open the file in MATLAB.
2. Adjust the following parameters as needed:
   - `N`: Resolution (e.g., `1024`, `2048`, `4096`).
   - `H`: Hurst Exponent (valid between `0` and `1`).
3. Run the script.

## 🖼 Examples
| Hurst Exponent ($H$) | Fractal Dimension ($D$) | Visual Style |
| :--- | :--- | :--- |
| **0.2** (picture h024096) | $D \approx 2.8$ | Volcanic, rocky, highly textured |
| **0.5** | $D \approx 2.5$ | Standard Brownian motion roughness |
| **0.8** (picture h084096) | $D \approx 2.2$ | Smooth mountain ridges and valleys |

## ⚠️ Performance Notes
Generating surfaces at $N = 4096$ creates ~16.7 million data points. For $N = 8192$ (the "White Whale" of resolution), ensure your system has at least 16GB of RAM. If you experience crashes, try casting variables to `single` precision.

---
*Generated with help from Gemini.*