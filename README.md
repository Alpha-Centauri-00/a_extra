# Spinning-Space Model for Galaxy Rotation Curves

## Overview

This project explores a fundamental question in cosmology:

> **Galaxies spin faster than expected from visible matter alone. Could the galaxy's spin itself modify gravity locally?**

Instead of invoking invisible "dark matter," we test a hypothesis: the **angular momentum (spin) of a galaxy creates a small extra acceleration** that naturally explains why galaxy rotation curves remain flat at large distances.

This repository contains:
- **Automated batch fitting** of the spin-coupling model to 175 real galaxies from the SPARC survey
- **Scaling law analysis** revealing how the coupling strength depends on galaxy size
- **Publication-quality visualizations** of the results

---

## The Problem: Galaxy Rotation Curves

### What We Observe

Astronomers measure how fast stars orbit at different distances from the galactic center. A typical galaxy shows:

$$v_{\text{observed}}(r) \approx \text{constant across all radii}$$

**This is surprising.** According to Newton's gravity and visible matter alone:

$$v_{\text{visible}}(r) = \sqrt{\frac{GM_{\text{visible}}(r)}{r}}$$

This predicts that outer stars should move slower (like planets in the solar system—Mercury orbits faster than Neptune). But they don't.

### The Standard Solution: Dark Matter

To explain flat rotation curves, cosmologists add an invisible halo of "dark matter" around each galaxy. This extra mass provides the gravity needed to keep outer stars moving fast.

**Problem:** Dark matter is hypothetical. We've never directly detected it. Can we explain the observations differently?

---

## Our Hypothesis: Spin-Mediated Acceleration

We propose that **the galaxy's own angular momentum couples to spacetime**, creating an extra acceleration that helps explain flat rotation curves.

### The Core Idea

The total acceleration at radius $r$ is a combination of two effects:

$$a_{\text{total}}(r) = a_{\text{visible}}(r) + a_{\text{extra}}(r)$$

Where:
- $a_{\text{visible}}(r)$ = acceleration from visible matter (gas, stars, bulge)
- $a_{\text{extra}}(r)$ = extra acceleration from the galaxy's spin

### The Spin-Coupling Term

We propose that the extra acceleration follows:

$$a_{\text{extra}}(r) = k \cdot \frac{J}{(r + r_0)^{\text{power}}}$$

**Physical parameters:**

| Symbol | Meaning | Value |
|--------|---------|-------|
| $k$ | Coupling constant | Fitted from data (range: $10^{-42}$ to $10^{-36}$) |
| $J$ | Total angular momentum of galaxy | Calculated from outermost data point |
| $r$ | Distance from galactic center | Measured in meters |
| $r_0$ | Softening radius | $\approx 0.5$ kpc (prevents singularity at center) |
| power | Power-law exponent | $2.0$ (inverse-square law) |

---

## The Model: Step by Step

### Step 1: Calculate Visible Matter Velocity

From the rotation curve data (which includes gas, disk, and bulge contributions), combine components:

$$v_{\text{vis}}(r) = \sqrt{v_{\text{gas}}^2(r) + v_{\text{disk}}^2(r) + v_{\text{bulge}}^2(r)}$$

### Step 2: Calculate Galaxy's Angular Momentum

From the outermost observed data point (radius $r_{\max}$, velocity $v_{\max}$):

**Total enclosed mass:**
$$M_{\text{enclosed}} = \frac{v_{\max}^2 \cdot r_{\max}}{G}$$

**Total angular momentum:**
$$J = M_{\text{enclosed}} \cdot v_{\max} \cdot r_{\max}$$

This $J$ is held constant across all radii and represents the galaxy's total spin.

### Step 3: Calculate Extra Acceleration

At each radius $r$:

$$a_{\text{extra}}(r) = k \cdot \frac{J}{(r + r_0)^{2}}$$

### Step 4: Predict Rotation Velocity

The model prediction combines visible matter and spin effects:

$$v_{\text{model}}^2(r) = v_{\text{vis}}^2(r) + r \cdot a_{\text{extra}}(r)$$

$$v_{\text{model}}(r) = \sqrt{v_{\text{vis}}^2(r) + r \cdot a_{\text{extra}}(r)}$$

### Step 5: Fit $k$ Parameter

For each galaxy, we test a range of $k$ values and find which one **minimizes the mean relative error**:

$$\text{Error} = \frac{1}{N} \sum_{i=1}^{N} \left| \frac{v_{\text{model}}(r_i) - v_{\text{obs}}(r_i)}{v_{\text{obs}}(r_i)} \right|$$

---

## Results from 175 SPARC Galaxies

### Fitting Performance

We applied this model to **175 real galaxies** from the SPARC survey (Lelli et al. 2016).

**Key statistics:**

| Metric | Value |
|--------|-------|
| **Number of galaxies** | 175 |
| **Mean fitting error** | 30.26% |
| **Median fitting error** | 28.12% |
| **Best fit (F563-V1)** | 5.70% error |
| **Worst fit (UGC02455)** | 107.91% error |
| **Standard deviation** | 14.78% |

![Error Distribution (175 galaxies)](https://github.com/Alpha-Centauri-00/a_extra/blob/main/results/phase6_error_distribution.png)

---
![Error vs Galaxy Size (175 galaxies)](https://github.com/Alpha-Centauri-00/a_extra/blob/main/results/phase6_error_vs_J_full_sample.png)

---
![Error vs Galaxy Size (175 galaxies)](https://github.com/Alpha-Centauri-00/a_extra/blob/main/results/batch_fitting_csv.png)

---

### Discovered New Scaling ?

The fitted coupling constant $k$ is **not random**—it depends systematically on galaxy size:

$$k = 0.3905 \times J^{-0.589}$$

With $R^2 = 0.152$.

**Interpretation:**
- **Small, low-spin galaxies** (low $J$) → larger $k$ → stronger coupling
- **Large, high-spin galaxies** (high $J$) → smaller $k$ → weaker coupling

This is physically intuitive: massive rotating systems have more inertia, so they're less susceptible to the same geometric coupling effect.

### Comparison to Alternatives

| Aspect | This Model | Dark Matter | MOND |
|--------|-----------|------------|------|
| **Rotation curves** | ✓ 30% mean error | ✓ 12-20% | ✓ 10-15% |
| **Galaxy clusters** | ✗ No | ✓ Yes | ✗ Weak |
| **Gravitational lensing** | ✗ No | ✓ Yes | ✗ No |
| **CMB anomalies** | ✗ No | ✓ Yes | ✗ Partial |
| **New particles** | ✗ No (geometry only) | ✓ Hypothetical | ✗ No (modified gravity) |
| **Simplicity** | ✓ 2 parameters | ✗ Many parameters | ✓ Simple modification |

---

## How to Use

> **Note:** The results are already available in the "results" folder. 

But if you want to generate them again, just clone this repo and run:

```bash
uv sync
```

---

## Citation

---

If you use this code or model in your research:

```bibtex
@software{spinning-space-2025,
  title={Spinning-Space Model: Automated Batch Fitting of Galaxy Rotation Curves},
  author={m.kherki},
  year={2025},
  url={https://github.com/yourusername/spinning-space-model}
}
```

---

## License

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).

[![CC BY-NC-ND 4.0](https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

## Zenodo DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18075346.svg)](https://doi.org/10.5281/zenodo.18075346)

---
## Data Resource
[Rotmod_LTG SPARC](https://astroweb.case.edu/SPARC/)


## Questions or Issues?

- Open an issue on GitHub

- Discussions welcome!

