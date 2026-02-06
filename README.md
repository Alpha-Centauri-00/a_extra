# Spin-Mediated Acceleration in Galaxy Rotation Curves: A Phenomenological Analysis

## Overview

This project explores a fundamental question in cosmology:

> **Galaxies spin faster than expected from visible matter alone. Could the galaxy's spin itself modify gravity locally?**

Instead of invoking invisible "dark matter," we test a hypothesis: the **angular momentum (spin) of a galaxy creates a small extra acceleration** that naturally explains why galaxy rotation curves remain flat at large distances.


---

The flatness of galaxy rotation curves at large radii constitutes a fundamental challenge to Newtonian dynamics and luminous matter distributions. We propose a **phenomenological spin-coupling model** wherein orbital acceleration receives a geometric contribution from total galactic angular momentum:

$$
a_\text{extra}(r) = k \cdot \frac{J_\text{vis}}{(r + r_0)^2}
$$

with coupling constant $k$ [kg$^{-1}$m$^{-1}$s$^2$], softening radius $r_0 = 0.5$ kpc, and $J_\text{vis}$ constructed exclusively from SPARC baryonic velocity decomposition ($v_\text{gas}$, $v_\text{disk}$, $v_\text{bulge}$) to eliminate circularity.

**Fitting results for 175 SPARC galaxies** (Lelli et al. 2016):
- **Mean error**: **28.39%** (±14.26%), **median**: **26.03%**
- **Systematic trend**: $\log_{10}k$ vs $\log_{10}J_\text{vis}$ correlation **r = -0.447** (p$<$0.001)
- **Power-law fit**: $\log_{10}k = -1.306 \times \log_{10}J_\text{vis} + C$ (**R$^2$ = 0.200**)
- **Best**: UGC09992 (**4.60%**), **worst**: UGC02455 (**107.91%**)

The model exhibits **morphological selectivity**, performing best for pure-disk systems ($\sim$5-10% error) and degrading substantially in bulge-dominated galaxies ($>$80%). Fitted $k$ values ($10^{-45}$ to $10^{-34}$) exceed naive Lense-Thirring predictions by 4-6 dex but follow the expected $J/r^2$ functional form, suggesting disk-integration and secular coherence as plausible amplification mechanisms.

While less precise than multi-parameter MOND/CDM fits (10-20% errors), this **one-parameter-per-galaxy model** demonstrates coherent systematics across seven orders of magnitude in $J_\text{vis}$, motivating further investigation of baryonic spin-geometry couplings as contributors to galactic dynamics.

**Keywords:** galaxy rotation curves, frame-dragging, SPARC survey, angular momentum, geometric gravity




---

## How to Use

Clone repo and run:

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
  url={https://github.com/Alpha-Centauri-00/a_extra}
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


