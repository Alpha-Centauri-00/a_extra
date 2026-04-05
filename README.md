# Spin-Mediated Acceleration as a Phenomenological Fitting Function for Galaxy Rotation

## Overview

This project applies a simple phenomenological fitting function to galaxy rotation curves. The model adds an angular-momentum-dependent acceleration term to the baryonic contribution and tests whether this improves fits to 175 SPARC disk galaxies.

---

## Model

The total circular velocity is modelled as:

$$
v_\text{model}^2(r) = v_\text{vis}^2(r) + r \cdot a_\text{extra}(r)
$$

where the extra acceleration term is:

$$
a_\text{extra}(r) = k \cdot \frac{J_\text{vis}}{(r + r_0)^2}
$$

- $k$ — coupling constant [kg⁻¹ m⁻¹ s²], fitted per galaxy or globally
- $J_\text{vis}$ — visible-matter angular momentum computed from baryonic components only ($v_\text{gas}$, $v_\text{disk}$, $v_\text{bulge}$), evaluated at the half-mass radius — **not** from the observed curve
- $r_0 = 2 \times r_\text{half}$ per galaxy (robustness-scan optimum; previously fixed at 0.5 kpc)

---

## Updated Results (Version 2)

**Per-galaxy fit — primary result:**

| Metric | Value |
|---|---|
| Galaxies | 175 SPARC |
| Mean error | **18.52% ± 14.96%** |
| Median error | **15.76%** |
| Best fit | UGC 00731 (1.86%) |
| Worst fit | UGC 02455 (107.91%) |

**Global model comparison:**

| | Model A (universal $k$) | Model B (power-law) |
|---|---|---|
| Parameters | 1 | 2 |
| $k_0$ | $1.49 \times 10^{-36}$ | $9.42 \times 10^{7}$ |
| $\beta$ | 0 | **−0.652** |
| Mean error | 32.60% | 34.52% |
| Median error | 29.62% | 23.34% |
| ΔAIC | — | **−141,371** (B preferred) |

**Two-population structure:** ~11% of galaxies (20/175) converge to $k \approx 0$ — the baryonic mass model already reproduces those curves without coupling. The global power-law cannot fit both populations simultaneously.

**Morphological correlations (weak):**
- $\log_{10} k$ vs $\log_{10} J_\text{vis}$: $r = -0.405$, $R^2 \approx 0.16$
- $\log_{10} k$ vs disk scale length $h$: $r = -0.249$
- $\log_{10} k$ vs B/T: not significant

**Cluster test:** Mean error 72.5% on 6 galaxy clusters — model correctly fails for pressure-supported systems.

**Lense-Thirring comparison:** Fitted $k$ values ($10^{-36}$ to $10^{-32}$) exceed bare LT predictions ($10^{-46}$ to $10^{-48}$) by **10–14 orders of magnitude**. No rigorous derivation of this gap exists. The LT functional form motivates the $J/r^2$ scaling but cannot be used to claim a physical mechanism.

---

## How to Use

Clone and install:

```bash
uv sync
```

---

## Citation

If you use this code or model in your research:

```bibtex
@software{spinning-space-2026,
  title={Spin-Mediated Acceleration as a Phenomenological Fitting Function for Galaxy Rotation},
  author={Kherki, M.},
  year={2026},
  version={2.0},
  url={https://github.com/Alpha-Centauri-00/a_extra}
}
```

---

## License

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)


This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).

[![CC BY-NC-ND 4.0](https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

## Zenodo DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19428420.svg)](https://doi.org/10.5281/zenodo.19428420)

---
## Data Resource
[Rotmod_LTG SPARC](https://astroweb.case.edu/SPARC/)


## Questions or Issues?

- Open an issue on GitHub

- Discussions welcome!


