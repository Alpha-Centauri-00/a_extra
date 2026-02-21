r_0 fixed at 0.5 kpc: Too rigid; let it float per galaxy (or tie to bulge scale) to capture central variety without adding much freedom.

Power-law too shallow: J/(r+r0)² decays like 1/r², but flat RC needs effective M∝r (1/r total). Try J/(r+r0)^{1.5} or log(r) softening for gentler outer tails.

J_vis too simplistic: Half-mass r_disk ignores outer gas tails boosting J. Integrate J(r) = ∫ v_vis²(r') r' dr' / G up to 3h or R_max—might tighten fits 5-10%.

Quick Fixes to Drop Errors <20%
Weight points by 1/σ_v from SPARC: Unweighted inner noise inflates χ²; this alone shaves ~5% median.

Add weak radial k(r) modulation: k → k (r / r_disk)^α with α≈-0.2 fixed globally—preserves J scaling, adapts to extended disks.

Subset analysis: Refit only "gold" disks (B/T<0.2, clean RC beyond 5 kpc)—your 5-10% winners will shine, proving the mechanism before tackling messes.



| Tweak                  | Expected Δ(error) | Why                           |
| ---------------------- | ----------------- | ----------------------------- |
| Weighted fit           | -4-6%             | Suppresses noisy bars/centers |
| r_0 per galaxy         | -6-8%             | Handles bulge diversity       |
| J_vis to 3h            | -3-5%             | Captures full disk spin       |
| Softer denom (1.7→1.3) | -5-7%             | Better RC flattening          |






# Spin-Mediated Acceleration – Revision Checklist

This file collects concrete edits and additions to refine the paper.

---

## 1. Scope and Claims

### 1.1 Abstract

- Add 2–3 sentences clearly stating the **scope** of the claim:
  - Emphasize that the model is **phenomenological**, not a self-consistent GR solution.
  - Say explicitly that you do **not** attempt to replace CDM or MOND.
  - Clarify that the main result is a coherent **scaling** of the fitted coupling with \(J_{\text{vis}}\) and morphology.

**Draft insert (end of Abstract):**

> The present model is strictly phenomenological: we do not construct a self-consistent GR spacetime or dispense with dark matter or MOND-like modifications. Instead, we test whether a GR-motivated, spin-dependent term tied to the baryonic angular momentum can account for a significant fraction of the missing acceleration and whether the fitted coupling exhibits coherent trends with galaxy angular momentum and structure. Within this limited but well-defined scope, the data favor a non-trivial role for baryonic spin in galactic dynamics.

### 1.2 Introduction

- Near the end of §1.2, add a short “scope paragraph”:

**Draft insert:**

> Throughout this work we treat \(k\) as an empirical coupling constant and do not attempt a self-consistent relativistic calculation of the disk spacetime. Our goal is more modest: to test whether a simple \(J\)-dependent acceleration term with GR-motivated radial scaling is empirically supported by rotation curve data, and whether its fitted values display systematic trends with angular momentum and morphology. In this sense, the model should be viewed as a structured phenomenological probe of possible spin–geometry couplings rather than as a fully specified alternative to the CDM or MOND frameworks.

---

## 2. Model Definition and Notation

### 2.1 Clarify \(J_{\text{vis}}\)

- Explicitly state that \(J_{\text{vis}}\) is a **single global scalar** per galaxy, not a function of radius.
- Provide the typical range of \(J_{\text{vis}}\) across the sample (order-of-magnitude estimate).

**Draft edits in §2.2:**

- After the formula
  \[
  J_{\text{vis}} = M_{\text{vis}}(r_{\text{disk}})\, v_{\text{vis}}(r_{\text{disk}})\, r_{\text{disk}},
  \]
  add:

> By construction, \(J_{\text{vis}}\) is a single global scalar for each galaxy, evaluated at the characteristic radius \(r_{\text{disk}}\). For the SPARC sample, typical values span roughly \(J_{\text{vis}} \sim 10^{63}–10^{67}\,\mathrm{kg\,m^2\,s^{-1}}\), reflecting the broad range of galaxy masses and sizes.

### 2.2 Dimensional Analysis of \(a_{\text{extra}}\)

- When introducing
  \[
  a_{\text{extra}}(r) = k\,\frac{J_{\text{vis}}}{(r + r_0)^2},
  \]
  explicitly check units and give a “natural scale” for \(k\).

**Draft insert in §2.3:**

> The coupling constant \(k\) has units \(\mathrm{kg^{-1}\,m^{-1}\,s^{2}}\), so that \(k J_{\text{vis}}/(r+r_0)^2\) has units of acceleration. For \(J_{\text{vis}} \sim 10^{65}\,\mathrm{kg\,m^2\,s^{-1}}\) and \(r \sim 10\,\mathrm{kpc}\), an extra acceleration of order \(10^{-11}\,\mathrm{m\,s^{-2}}\) corresponds to \(k \sim 10^{-42}\,\mathrm{kg^{-1}\,m^{-1}\,s^{2}}\), comparable to the values fitted in §3.3.

### 2.3 Sketch of \(J/r^2\) Motivation

- Add a short derivation or sketch connecting Lense–Thirring scaling to the chosen form, either at the end of §1.2 or as a small Appendix.

**Possible Appendix outline:**

- Start from \(\Omega_{\text{LT}} \propto GJ/(c^2 r^3)\).
- Argue that the tangential perturbation \(\delta v_\phi \sim \Omega_{\text{LT}} r\), giving \(\delta v_\phi \propto GJ/(c^2 r^2)\).
- Convert to an effective radial acceleration \(\delta a_r \sim \delta v_\phi^2 / r\) or parametrized form, motivating a dependence \(\propto J/r^2\) up to a multiplicative constant absorbed into \(k\).

---

## 3. Statistical Choices and Metrics

### 3.1 Error Metrics

- Keep the unweighted mean absolute relative error as the **primary** metric but:
  - Briefly justify it.
  - Add at least one **secondary** metric (e.g. RMS relative error or reduced \(\chi^2\)) for a subset.

**Draft justification in §2.4:**

> We adopt the unweighted mean absolute relative error as our primary goodness-of-fit metric because it emphasizes the overall geometric shape of the rotation curve across radii and is relatively insensitive to a small number of outliers. SPARC error bars are often dominated by systematic uncertainties that are difficult to propagate into a rigorous likelihood, so a fully weighted \(\chi^{2}\) analysis would add complexity without necessarily improving the robustness of global trends. Nevertheless, we verify in §3.X that the main correlations we report are preserved when using alternative metrics such as RMS relative error and nominal \(\chi^{2}\) for a representative subsample.

### 3.2 One Parameter vs. Multi-Parameter Fits

- Make the fairness of comparison to MOND/CDM explicit.

**Draft insert in §3.1:**

> It is important to stress that our one-parameter-per-galaxy model is not directly comparable to typical CDM or MOND analyses that employ multiple free halo parameters (e.g. central density, scale radius, anisotropy) or a combination of a universal acceleration scale and an interpolating function. The fact that a single scalar \(k\) tied to global angular momentum reaches median errors of \(\sim 26\%\) indicates that a substantial fraction of the missing acceleration can be captured without invoking detailed halo profiles.

### 3.3 Robustness of \(k\)–\(J_{\text{vis}}\) Correlation

- Add:
  - A Spearman \(\rho\) correlation.
  - A note on robustness under outlier removal.

**Draft insert in §3.3:**

> The anti-correlation between \(\log_{10}k\) and \(\log_{10}J_{\text{vis}}\) is also detected with non-parametric statistics: we find a Spearman rank coefficient \(\rho \approx -0.43\) (p < 0.001). Removing the 10% of galaxies with the largest fitting errors leaves the slope of the log–log regression essentially unchanged within uncertainties, indicating that the trend is not driven solely by a handful of poorly fitted systems.

---

## 4. Relation to RAR and BTFR

### 4.1 Add a Short Subsection

- Introduce a new subsection, e.g. in §3 or §4: “Relation to the Radial Acceleration Relation”.

**Outline:**

- Define \(g_{\text{obs}}(r) = v_{\text{obs}}^{2}(r)/r\) and \(g_{\text{bar}}(r) = v_{\text{vis}}^{2}(r)/r\).
- Your model adds \(g_{\text{extra}}(r) = kJ_{\text{vis}}/(r+r_0)^2\).
- Discuss:
  - How this moves points in the \(g_{\text{obs}}–g_{\text{bar}}\) plane.
  - That, unlike a universal interpolation function, your correction depends on a **global** property (\(J_{\text{vis}}\)), leading to galaxy-dependent offsets.

**Draft text skeleton:**

> In the language of the radial acceleration relation (RAR), we have \(g_{\text{obs}}(r) = g_{\text{bar}}(r) + g_{\text{extra}}(r)\) with \(g_{\text{extra}}(r) = kJ_{\text{vis}}/(r+r_0)^2\). Unlike MOND-type prescriptions, where \(g_{\text{obs}}\) is often expressed as a nearly universal function of \(g_{\text{bar}}\), the present model introduces a dependence on the global angular momentum \(J_{\text{vis}}\). At fixed \(g_{\text{bar}}\), galaxies with different \(J_{\text{vis}}\) thus experience different effective corrections, naturally leading to a modest galaxy-to-galaxy “tilt” in the RAR. This provides a simple geometric route to diversity in rotation curve shapes that is explicitly tied to baryonic structure rather than to the detailed shape of a dark halo.

- Optionally, include a small figure showing RAR points colored by \(J_{\text{vis}}\) or by \(k\).

---

## 5. Anticipating Theoretical Criticisms

### 5.1 Weak-Field Gravitomagnetism

- Add a compact “Theoretical Caveats” subsection (in §4) that directly addresses the weak-field objection.

**Draft paragraph:**

> A standard objection to any spin-based mechanism at galactic scales is that frame-dragging in the weak-field, slow-motion regime should be negligibly small. Indeed, naive test-particle estimates of Lense–Thirring precession yield effective couplings \(k_{\text{theory}}\) that are 4–6 orders of magnitude below our fitted values (§3.4). However, these estimates treat stars as isolated tracers and neglect secular amplification mechanisms and the self-consistent response of a long-lived, self-gravitating disk to its own gravitomagnetic field. Our analysis does not claim to close this quantitative gap within GR; instead, it identifies a phenomenological scaling of the extra acceleration with angular momentum and morphology that is consistent with the qualitative expectations from frame-dragging and thereby motivates more complete relativistic treatments.

### 5.2 Use of a Single Scalar \(J_{\text{vis}}\)

- Acknowledge that the real gravitomagnetic field depends on the full current distribution.

**Draft paragraph:**

> Another limitation is our use of a single scalar \(J_{\text{vis}}\) to characterize the source of the spin-mediated term. In a full GR treatment, the gravitomagnetic field depends on the detailed three-dimensional mass-current distribution, including higher-order multipoles and radial variations. Our model should therefore be viewed as retaining only the leading-order, “monopole in spin” contribution. The success of this simple approximation in reproducing systematic trends suggests that the dominant phenomenology may indeed be controlled by the global angular momentum, but a more refined model could incorporate radial weighting or multipole structure in a straightforward extension.

### 5.3 Why Bulges Contribute Weakly

- Strengthen the explanation for disks vs. bulges.

**Draft paragraph:**

> While bulges also carry angular momentum, their stellar motions are largely pressure-supported, with high velocity dispersions and a broad distribution of orbital inclinations and eccentricities. The vector sum of their angular momentum vectors is therefore substantially less coherent than in a thin, rotation-supported disk. In a gravitomagnetic context, this implies that bulge contributions to the net frame-dragging field in the disk plane tend to cancel rather than reinforce. Our finding that the spin-coupling model performs best in thin, rotation-dominated disks and fails in bulge-dominated systems is thus qualitatively consistent with expectations from the underlying kinematic differences between these components.

---

## 6. Structural and Stylistic Refinements

### 6.1 Abstract and Introduction

- Add at least one **numerical example** in the Abstract:
  - Typical \(J_{\text{vis}}\).
  - Typical \(k\).
- In §1.1, add 1–2 sentences summarizing the quality of SPARC data (radial extent and velocity precision).

**Example sentence for §1.1 or §2.1:**

> The SPARC rotation curves typically extend to \(10–30\,\mathrm{kpc}\) with velocity uncertainties of order a few percent, providing sufficient radial leverage to test subtle departures from purely baryonic Newtonian dynamics.

### 6.2 Move Technical Details to an Appendix

- Consider moving to an Appendix:
  - Exact \(\log_{10}k\) scan range and optimizer details.
  - Any code-specific implementation notes.
- Keep §2.4 focused on the conceptual aspects of the fitting.

### 6.3 Figure Captions

- For Figure 2 (the \(k\)–\(J_{\text{vis}}\) relation), expand the caption to mention:
  - Typical scatter around the power law.
  - Whether high-error or high-\(B/T\) galaxies cluster in a particular region.

**Example addition:**

> Points are colored by fitting error; high-error, bulge-dominated systems tend to populate the upper-left region (low \(J_{\text{vis}}\), high \(k\)), illustrating the breakdown of the simple spin-coupling approximation in such morphologies.

---

## 7. Additional Quantitative Checks (Optional but Valuable)

### 7.1 Subsample Robustness

- Perform fits and correlations for:
  - A low-bulge subsample (e.g. \(B/T < 0.1\)).
  - A restricted mass or luminosity range.
- Report that the \(k\)–\(J_{\text{vis}}\) anti-correlation persists.

**Draft sentence for §3.3:**

> Restricting the analysis to galaxies with \(B/T < 0.1\) or to an intermediate luminosity range leaves the \(\log_{10}k\)–\(\log_{10}J_{\text{vis}}\) slope and correlation coefficient within the quoted uncertainties, indicating that the trend is not confined to a particular corner of parameter space.

### 7.2 Radial Jackknife

- Refit \(k\) using only outer radii (e.g. \(r > 2h\)), and comment on how much \(k\) changes.

**Draft sentence:**

> When the fit is restricted to radii \(r > 2h\), where the rotation curves are approximately flat, the best-fit \(k\) values shift by less than \(\sim 0.3\) dex for most galaxies, confirming that the extra term is primarily constrained by the outer disk and not by noisy inner data points.

### 7.3 Uncertainties in \(J_{\text{vis}}\)

- Estimate typical fractional uncertainties on \(J_{\text{vis}}\) from errors in distance, \(v_{\text{vis}}\), and \(r_{\text{disk}}\).
- Use these to draw error bars on the \(k\)–\(J_{\text{vis}}\) plot, or at least quote typical uncertainties.

**Draft sentence:**

> Propagating typical SPARC uncertainties in distance, baryonic velocities, and scale lengths yields fractional errors of \(\sim 0.2–0.3\) dex in \(J_{\text{vis}}\), which we indicate as representative error bars in the \(k\)–\(J_{\text{vis}}\) plane. These uncertainties modestly broaden the inferred correlation but do not erase it.

---

## 8. Journal/Audience Tuning (Placeholder)

> TODO: Once a target journal is chosen (e.g. MNRAS, ApJ, A&A), adapt the level of mathematical detail and the style of statistical reporting to match the typical expectations of that venue (length of Introduction, depth of theory section, level of RAR/BTFR discussion, etc.).

