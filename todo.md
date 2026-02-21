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