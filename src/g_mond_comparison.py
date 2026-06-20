"""
MOND Comparison: RAR interpolating function vs spin-mediated model
==================================================================

Uses McGaugh et al. (2016) / Li et al. (2018) simple interpolating function:
    g_obs = g_bar / (1 - exp(-sqrt(g_bar / g_dagger)))
    g_dagger = 1.20e-10 m/s^2  (fixed)

Upsilon_disk and Upsilon_bul are optimised per galaxy via L-BFGS-B with
Gaussian log-space priors centred at 0.5 and 0.7 M_sun/L_sun (sigma = 0.1 dex),
following Li et al. (2018, A&A, 615, A3).

Also performs two-population Model B analysis: reports Model B stats separately
for k>0 vs k≈0 galaxies (threshold: log10(k) <= K_ZERO_THRESHOLD_LOG10).

Outputs
-------
results/mond_comparison_results.csv  -- per-galaxy MOND fit results (+ log10_k, k_zero, spin_mean_error_pct)
results/mond_error_histogram.png     -- overlaid error distributions
results/mond_vs_spin_scatter.png     -- MOND error vs spin error per galaxy, coloured by k≈0 status
results/mond_rotation_curves.png     -- rotation curve overlays (fill PLOT_GALAXIES first)
results/two_population_model_b.txt   -- Model B stats split by population

Run AFTER c_fit_k_continuous.py.
"""

import os
import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
from scipy.optimize import minimize
from tools.config import constants


# =============================================================================
# CONFIGURATION
# =============================================================================

G_DAGGER               = 1.20e-10        # RAR acceleration scale [m/s^2]
LOG_UD_PRIOR           = math.log10(0.5) # -0.30103  Li et al. (2018)
LOG_UB_PRIOR           = math.log10(0.7) # -0.15490  Li et al. (2018)
PRIOR_SIGMA_DEX        = 0.1             # dex
UPSILON_LOG_BOUNDS     = (-1.0, 1.0)     # Υ in [0.1, 10] in log10 space
K_ZERO_THRESHOLD_LOG10 = -45.0           # log10(k) below this → k≈0 population

# Fill in after running the script and reviewing mond_comparison_results.csv.
# Example: PLOT_GALAXIES = ["NGC 0024", "UGC 00731", "F563-V1", "NGC 2976"]
PLOT_GALAXIES: List[str] = []


# =============================================================================
# DATA LOADING
# =============================================================================

class RotationCurveLoader:
    """Load and parse rotation curve .dat files (SPARC format)."""

    @staticmethod
    def load(filepath: str) -> Optional[Dict]:
        result = {
            "distance_mpc": None,
            "radii": [], "v_obs": [], "v_err": [],
            "v_gas": [], "v_disk": [], "v_bul": [],
        }
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if line.startswith("# Distance") and "=" in line:
                    try:
                        result["distance_mpc"] = float(
                            line.split("=")[1].strip().split()[0])
                    except Exception:
                        pass
                    break

            for line in lines:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                try:
                    p = s.split()
                    if len(p) < 3:
                        continue
                    result["radii"].append(float(p[0]))
                    result["v_obs"].append(float(p[1]))
                    result["v_err"].append(max(float(p[2]) if len(p) > 2 else 1.0, 1.0))
                    result["v_gas"].append(float(p[3])  if len(p) > 3 else 0.0)
                    result["v_disk"].append(float(p[4]) if len(p) > 4 else 0.0)
                    result["v_bul"].append(float(p[5])  if len(p) > 5 else 0.0)
                except (ValueError, IndexError):
                    continue

        except FileNotFoundError:
            return None

        return result if result["radii"] else None


# =============================================================================
# MOND MODEL
# =============================================================================

def mond_g_obs(g_bar: float) -> float:
    """
    RAR simple interpolating function (McGaugh et al. 2016):
        g_obs = g_bar / (1 - exp(-sqrt(g_bar / g_dagger)))

    Deep-MOND limit (g_bar → 0): g_obs → sqrt(g_bar * g_dagger)
    Newtonian limit (g_bar >> g_dagger): g_obs → g_bar
    """
    if g_bar <= 0.0:
        return 0.0
    x = math.sqrt(g_bar / G_DAGGER)
    if x < 1e-7:
        return math.sqrt(g_bar * G_DAGGER)
    denom = 1.0 - math.exp(-x)
    if denom < 1e-15:
        return math.sqrt(g_bar * G_DAGGER)
    return g_bar / denom


def mond_velocity(r_kpc: float,
                  v_gas_kms: float, v_disk_kms: float, v_bul_kms: float,
                  upsilon_disk: float, upsilon_bul: float) -> float:
    """
    MOND circular velocity at radius r [km/s].

    g_bar = (Υ_d * v_disk^2 + Υ_b * v_bul^2 + v_gas^2) / r
    V_MOND = sqrt(g_obs(g_bar) * r)

    SPARC .dat files store v_disk and v_bul at Υ = 1, so the mass-to-light
    ratios enter as direct multipliers on v^2.  v_gas^2 already includes
    the 1.33 helium correction baked into the SPARC photometry pipeline.
    """
    r_m = r_kpc * constants.KPC_TO_M
    if r_m <= 0.0:
        return 0.0
    vd = v_disk_kms * 1e3
    vb = v_bul_kms  * 1e3
    vg = v_gas_kms  * 1e3
    g_bar = (upsilon_disk * vd**2 + upsilon_bul * vb**2 + vg**2) / r_m
    return math.sqrt(max(mond_g_obs(g_bar) * r_m, 0.0)) / 1e3


# =============================================================================
# SPIN MODEL HELPER  (mirrors c_fit_k_continuous.py, used for overlay plot)
# =============================================================================

def spin_velocity(r_kpc: float,
                  v_gas_kms: float, v_disk_kms: float, v_bul_kms: float,
                  k: float, J: float, r0_kpc: float) -> float:
    """Spin-mediated model circular velocity at radius r [km/s]."""
    r_m      = r_kpc  * constants.KPC_TO_M
    r0_m     = r0_kpc * constants.KPC_TO_M
    v_vis_ms = math.sqrt(v_gas_kms**2 + v_disk_kms**2 + v_bul_kms**2) * 1e3
    a_extra  = k * J / (r_m + r0_m) ** 2
    v_sq     = v_vis_ms**2 + r_m * a_extra
    if v_sq <= 0:
        return v_vis_ms / 1e3
    return math.sqrt(v_sq) / 1e3


# =============================================================================
# FIT METRICS  (no prior penalties — for reporting)
# =============================================================================

def mond_chi2_red(galaxy_data: Dict, ud: float, ub: float) -> float:
    """Reduced chi2 without prior penalties."""
    total, n = 0.0, 0
    for i, r in enumerate(galaxy_data["radii"]):
        vo = galaxy_data["v_obs"][i]
        ve = galaxy_data["v_err"][i]
        if vo <= 0 or ve <= 0:
            continue
        vm = mond_velocity(r, galaxy_data["v_gas"][i], galaxy_data["v_disk"][i],
                           galaxy_data["v_bul"][i], ud, ub)
        total += ((vm - vo) / ve) ** 2
        n += 1
    return total / n if n > 0 else 1e10


def mond_mare(galaxy_data: Dict, ud: float, ub: float) -> float:
    """Unweighted mean absolute relative error (for reporting)."""
    errors = []
    for i, r in enumerate(galaxy_data["radii"]):
        vo = galaxy_data["v_obs"][i]
        if vo <= 0:
            continue
        vm = mond_velocity(r, galaxy_data["v_gas"][i], galaxy_data["v_disk"][i],
                           galaxy_data["v_bul"][i], ud, ub)
        errors.append(abs(vm - vo) / vo)
    return float(np.mean(errors)) if errors else 1e10


# =============================================================================
# PER-GALAXY MOND FITTER
# =============================================================================

class MondFitter:
    """
    Optimise Υ_disk (and Υ_bul where a bulge exists) per galaxy using L-BFGS-B
    with Gaussian log-space priors, following Li et al. (2018).
    """

    def fit(self, galaxy_data: Dict) -> Dict:
        has_bulge = any(v > 0 for v in galaxy_data["v_bul"])

        if has_bulge:
            x0     = np.array([LOG_UD_PRIOR, LOG_UB_PRIOR])
            bounds = [UPSILON_LOG_BOUNDS, UPSILON_LOG_BOUNDS]
        else:
            x0     = np.array([LOG_UD_PRIOR])
            bounds = [UPSILON_LOG_BOUNDS]

        def objective(x: np.ndarray) -> float:
            log_ud = x[0]
            log_ub = x[1] if has_bulge else LOG_UB_PRIOR
            ud_    = 10.0 ** log_ud
            ub_    = 10.0 ** log_ub
            chi2   = mond_chi2_red(galaxy_data, ud_, ub_)
            pen    = ((log_ud - LOG_UD_PRIOR) / PRIOR_SIGMA_DEX) ** 2
            if has_bulge:
                pen += ((log_ub - LOG_UB_PRIOR) / PRIOR_SIGMA_DEX) ** 2
            return chi2 + pen

        res    = minimize(objective, x0=x0, method="L-BFGS-B", bounds=bounds,
                          options={"ftol": 1e-12, "maxiter": 5000})
        log_ud = float(res.x[0])
        log_ub = float(res.x[1]) if has_bulge else LOG_UB_PRIOR
        ud     = 10.0 ** log_ud
        ub     = 10.0 ** log_ub if has_bulge else 0.0

        chi2 = mond_chi2_red(galaxy_data, ud, ub if has_bulge else 0.0)
        mare = mond_mare(galaxy_data, ud, ub if has_bulge else 0.0)

        return {
            "upsilon_disk":       ud,
            "log10_upsilon_disk": log_ud,
            "upsilon_bul":        ub,
            "log10_upsilon_bul":  log_ub if has_bulge else None,
            "has_bulge":          has_bulge,
            "chi2_reduced":       chi2,
            "mean_error":         mare,
            "mean_error_pct":     mare * 100.0,
        }


# =============================================================================
# BATCH PROCESSOR
# =============================================================================

class BatchProcessor:

    def __init__(self, data_dir: str):
        self.data_dir = data_dir
        self.results: List[Dict] = []

    def _find_files(self) -> List[Tuple[str, str]]:
        files = []
        for fname in sorted(os.listdir(self.data_dir)):
            if fname.endswith("_rotmod.dat"):
                gname = fname.replace("_rotmod.dat", "")
                files.append((os.path.join(self.data_dir, fname), gname))
        return files

    def run(self) -> List[Dict]:
        files  = self._find_files()
        fitter = MondFitter()
        print(f"\n  {len(files)} galaxies found\n")
        for idx, (fp, gname) in enumerate(files, 1):
            print(f"  ({idx:3d}/{len(files)}) {gname} ...", end=" ")
            data = RotationCurveLoader.load(fp)
            if data is None:
                print("LOAD FAILED")
                continue
            fit = fitter.fit(data)
            row = {"galaxy_name": gname, **fit}
            self.results.append(row)
            bulge_str = f"  Υ_b={fit['upsilon_bul']:.3f}" if fit["has_bulge"] else ""
            print(f"Υ_d={fit['upsilon_disk']:.3f}{bulge_str}"
                  f"  chi2={fit['chi2_reduced']:.3f}"
                  f"  err={fit['mean_error_pct']:.2f}%")
        return self.results

    def save_csv(self, path: str):
        if not self.results:
            return
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fieldnames = [
            "galaxy_name",
            "upsilon_disk", "log10_upsilon_disk",
            "upsilon_bul",  "log10_upsilon_bul",
            "has_bulge",
            "chi2_reduced", "mean_error", "mean_error_pct",
            "log10_k", "k_zero", "spin_mean_error_pct",
        ]
        with open(path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.results)
        print(f"\n  ✓ MOND results → {path}")

    def print_summary(self):
        if not self.results:
            return
        errs  = [r["mean_error_pct"] for r in self.results]
        chi2s = [r["chi2_reduced"]   for r in self.results]
        print("\n" + "="*60)
        print("MOND FIT SUMMARY")
        print("="*60)
        print(f"  Galaxies     : {len(self.results)}")
        print(f"  Mean  error  : {np.mean(errs):.2f}% ± {np.std(errs):.2f}%")
        print(f"  Median error : {np.median(errs):.2f}%")
        print(f"  Best / worst : {min(errs):.2f}% / {max(errs):.2f}%")
        print(f"  Mean  chi2   : {np.mean(chi2s):.3f}")
        print(f"  Median chi2  : {np.median(chi2s):.3f}")
        print("="*60)


# =============================================================================
# LOAD SPIN MODEL RESULTS
# =============================================================================

def load_spin_results(per_galaxy_csv: str) -> Dict[str, Dict]:
    """Load per-galaxy k-fit results produced by c_fit_k_continuous.py."""
    data = {}
    try:
        with open(per_galaxy_csv, 'r') as f:
            for row in csv.DictReader(f):
                gname = row["galaxy_name"]
                data[gname] = {
                    "k":              float(row["k_optimal"]),
                    "log10_k":        float(row["log10_k_optimal"]),
                    "J":              float(row["J_new"]),
                    "r0_kpc":         float(row["r0_kpc"]),
                    "mean_error_pct": float(row["mean_error_pct"]),
                    "chi2_reduced":   float(row["chi2_reduced"]),
                }
    except FileNotFoundError:
        print(f"  Warning: {per_galaxy_csv} not found — run c_fit_k_continuous.py first.")
    return data


def load_model_b_errors(global_csv: str) -> Dict[str, float]:
    """Load per-galaxy Model B errors from global_fit_results.csv."""
    data = {}
    try:
        with open(global_csv, 'r') as f:
            for row in csv.DictReader(f):
                try:
                    data[row["galaxy_name"]] = float(row["error_pct_modelB"])
                except (KeyError, ValueError):
                    pass
    except FileNotFoundError:
        print(f"  Warning: {global_csv} not found — run c_fit_k_continuous.py first.")
    return data


# =============================================================================
# COMPARISON SUMMARY PRINTER
# =============================================================================

def print_comparison_summary(mond_results: List[Dict],
                              spin_results: Dict[str, Dict],
                              model_b_errors: Optional[Dict[str, float]] = None):
    common      = [r["galaxy_name"] for r in mond_results
                   if r["galaxy_name"] in spin_results]
    mond_errs   = [r["mean_error_pct"] for r in mond_results
                   if r["galaxy_name"] in common]
    spin_errs   = [spin_results[n]["mean_error_pct"] for n in common]
    mond_chi2s  = [r["chi2_reduced"]   for r in mond_results
                   if r["galaxy_name"] in common]
    spin_chi2s  = [spin_results[n]["chi2_reduced"]   for n in common]

    has_b = bool(model_b_errors)
    b_errs = [model_b_errors[n] for n in common if n in (model_b_errors or {})] if has_b else []
    # model B may cover fewer galaxies than the MOND/spin common set
    b_common = [n for n in common if n in (model_b_errors or {})]

    n   = len(common)
    n_b = len(b_common)
    w   = 14

    head = f"  {'Metric':<34} {'Per-galaxy k':>{w}}  {'MOND':>{w}}"
    sep  = "  " + "-"*(34 + 2*w + 4)
    if has_b and n_b:
        head += f"  {'Model B (global)':>{w}}"
        sep  += "-"*(w + 2)

    print("\n" + "="*(len(sep) - 2))
    print(f"COMPARISON SUMMARY  (N = {n} matched galaxies)")
    print("="*(len(sep) - 2))
    print(head)
    print(sep)

    def row(label, sv, mv, bv=None):
        s = f"  {label:<34} {sv:>{w}.2f}  {mv:>{w}.2f}"
        if has_b and n_b:
            s += f"  {bv:>{w}.2f}" if bv is not None else f"  {'n/a':>{w}}"
        print(s)

    row("Mean error (%)",   np.mean(spin_errs),   np.mean(mond_errs),
        np.mean(b_errs)   if b_errs else None)
    row("Median error (%)", np.median(spin_errs), np.median(mond_errs),
        np.median(b_errs) if b_errs else None)
    row("Std dev (%)",      np.std(spin_errs),    np.std(mond_errs),
        np.std(b_errs)    if b_errs else None)

    # Trimmed mean: exclude galaxies with MAE% > 100 % in EITHER model
    TRIM = 100.0
    trim_idx  = [i for i, (s, m) in enumerate(zip(spin_errs, mond_errs))
                 if s <= TRIM and m <= TRIM]
    trim_spin = [spin_errs[i] for i in trim_idx]
    trim_mond = [mond_errs[i] for i in trim_idx]
    n_trim    = len(trim_idx)
    if has_b and n_b:
        # re-align to b_common subset
        bcommon_set = set(b_common)
        trim_b_errs = [b_errs[b_common.index(common[i])]
                       for i in trim_idx if common[i] in bcommon_set]
        trim_b = np.mean(trim_b_errs) if trim_b_errs else None
    else:
        trim_b = None
    print(sep)
    row(f"Trimmed mean (MAE≤{TRIM:.0f}%, N={n_trim})",
        np.mean(trim_spin), np.mean(trim_mond), trim_b)

    print(sep)
    # chi2
    def row_chi2(label, sv, mv):
        s = f"  {label:<34} {sv:>{w}.3f}  {mv:>{w}.3f}"
        if has_b and n_b:
            s += f"  {'—':>{w}}"
        print(s)
    row_chi2("Mean chi2_red",   np.mean(spin_chi2s),   np.mean(mond_chi2s))
    row_chi2("Median chi2_red", np.median(spin_chi2s), np.median(mond_chi2s))

    print(sep)
    fp = f"  {'Free params':<34} {'1 per galaxy':>{w}}  {'1 per galaxy':>{w}}"
    if has_b and n_b:
        fp += f"  {'2 global (k0,β)':>{w}}"
    print(fp)
    print("="*(len(sep) - 2))
    if has_b and n_b and n_b < n:
        print(f"  * Model B column uses {n_b} galaxies with matched global-fit data.")


# =============================================================================
# TWO-POPULATION MODEL B ANALYSIS
# =============================================================================

def two_population_analysis(spin_results: Dict[str, Dict],
                             model_b_errors: Dict[str, float],
                             output_path: str):
    """
    Split galaxies into k≈0 and k>0 populations.
    Report Model B per-galaxy errors for each sub-population separately.
    """
    k_zero_names = {n for n, d in spin_results.items()
                    if d["log10_k"] <= K_ZERO_THRESHOLD_LOG10}
    k_pos_names  = {n for n, d in spin_results.items()
                    if d["log10_k"] > K_ZERO_THRESHOLD_LOG10}

    n_total = len(spin_results)
    n_zero  = len(k_zero_names)
    n_pos   = len(k_pos_names)

    all_b  = [model_b_errors[n] for n in spin_results   if n in model_b_errors]
    pos_b  = [model_b_errors[n] for n in k_pos_names    if n in model_b_errors]
    zero_b = [model_b_errors[n] for n in k_zero_names   if n in model_b_errors]

    lines = [
        "TWO-POPULATION MODEL B ANALYSIS",
        "=" * 60,
        f"k≈0 threshold  : log10(k) <= {K_ZERO_THRESHOLD_LOG10}",
        f"k≈0 galaxies   : {n_zero} / {n_total} ({100 * n_zero / n_total:.1f}%)",
        f"k>0 galaxies   : {n_pos}  / {n_total} ({100 * n_pos  / n_total:.1f}%)",
        "",
    ]

    if all_b:
        lines += [
            "Model B — ALL galaxies:",
            f"  Mean   error : {np.mean(all_b):.2f}% ± {np.std(all_b):.2f}%",
            f"  Median error : {np.median(all_b):.2f}%",
            "",
        ]
    if pos_b:
        lines += [
            f"Model B — k>0 only (N = {len(pos_b)}):",
            f"  Mean   error : {np.mean(pos_b):.2f}% ± {np.std(pos_b):.2f}%",
            f"  Median error : {np.median(pos_b):.2f}%",
            "",
        ]
    if zero_b:
        lines += [
            f"Model B — k≈0 only (N = {len(zero_b)}):",
            f"  Mean   error : {np.mean(zero_b):.2f}% ± {np.std(zero_b):.2f}%",
            f"  Median error : {np.median(zero_b):.2f}%",
        ]

    text = "\n".join(lines)
    print("\n" + text)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(text + "\n")
    print(f"\n  ✓ Two-population stats → {output_path}")


# =============================================================================
# PLOTS
# =============================================================================

def plot_error_histogram(mond_results: List[Dict],
                          spin_results: Dict[str, Dict],
                          output_path: str):
    """Overlaid per-galaxy error distribution: spin-mediated model vs MOND."""
    common = {r["galaxy_name"] for r in mond_results} & set(spin_results.keys())
    mond_errs = [r["mean_error_pct"] for r in mond_results
                 if r["galaxy_name"] in common]
    spin_errs = [spin_results[n]["mean_error_pct"] for n in common]

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.arange(0, 115, 5)

    ax.hist(spin_errs, bins=bins, alpha=0.55, color="tab:red",
            edgecolor="white", linewidth=0.5,
            label=(f"Spin-mediated  "
                   f"mean={np.mean(spin_errs):.1f}%  "
                   f"median={np.median(spin_errs):.1f}%"))
    ax.hist(mond_errs, bins=bins, alpha=0.55, color="tab:blue",
            edgecolor="white", linewidth=0.5,
            label=(f"MOND RAR       "
                   f"mean={np.mean(mond_errs):.1f}%  "
                   f"median={np.median(mond_errs):.1f}%"))

    for val, color, ls in [
        (np.mean(spin_errs),   "tab:red",  "--"),
        (np.median(spin_errs), "tab:red",  ":"),
        (np.mean(mond_errs),   "tab:blue", "--"),
        (np.median(mond_errs), "tab:blue", ":"),
    ]:
        ax.axvline(val, color=color, linestyle=ls, linewidth=1.3, alpha=0.85)

    ax.text(0.97, 0.97, "— mean   ···  median",
            transform=ax.transAxes, fontsize=8,
            ha="right", va="top", color="gray")

    ax.set_xlabel("Mean absolute relative error (%)", fontsize=12)
    ax.set_ylabel("Number of galaxies", fontsize=12)
    ax.set_title(
        f"Error distribution: spin-mediated model vs MOND RAR\n"
        f"({len(common)} SPARC galaxies, 1 free parameter per galaxy each)",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.set_xlim(0, 110)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Error histogram    → {output_path}")


def plot_rotation_curves(galaxy_names: List[str],
                          data_dir: str,
                          mond_map: Dict[str, Dict],
                          spin_map: Dict[str, Dict],
                          output_path: str):
    """
    Multi-panel rotation curve overlay for selected galaxies.
    Shows V_obs, V_MOND (fitted Υ), V_spin (fitted k), and V_bar (Υ=1).

    Fill PLOT_GALAXIES at the top of this file after reviewing the results CSV.
    """
    if not galaxy_names:
        print("  PLOT_GALAXIES is empty — skipping rotation curve plot.")
        print("  Fill in galaxy names at the top of the script and re-run.")
        return

    n     = len(galaxy_names)
    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5 * ncols, 4 * nrows),
                              squeeze=False)

    for idx, gname in enumerate(galaxy_names):
        ax   = axes[idx // ncols][idx % ncols]
        fp   = os.path.join(data_dir, f"{gname}_rotmod.dat")
        data = RotationCurveLoader.load(fp)

        if data is None:
            ax.set_title(f"{gname}\n(data not found)", fontsize=9)
            continue

        radii = np.array(data["radii"])

        # V_obs
        ax.errorbar(radii, data["v_obs"], yerr=data["v_err"],
                    fmt="ko", markersize=3, linewidth=0.8, capsize=2,
                    label=r"$V_{\rm obs}$", zorder=5)

        # V_bar at Υ=1 (reference, same normalization as spin model)
        v_bar = [math.sqrt(data["v_gas"][i]**2 + data["v_disk"][i]**2
                           + data["v_bul"][i]**2)
                 for i in range(len(radii))]
        ax.plot(radii, v_bar, color="gray", linestyle=":", linewidth=1.2,
                label=r"$V_{\rm bar}$ ($\Upsilon=1$)", zorder=2)

        # V_MOND
        if gname in mond_map:
            mr     = mond_map[gname]
            ud     = mr["upsilon_disk"]
            ub     = mr["upsilon_bul"]
            v_mond = [mond_velocity(radii[i],
                                    data["v_gas"][i], data["v_disk"][i],
                                    data["v_bul"][i], ud, ub)
                      for i in range(len(radii))]
            ax.plot(radii, v_mond, color="tab:blue", linestyle="-", linewidth=1.8,
                    label=(rf"MOND  $\Upsilon_d={ud:.2f}$"
                           rf"  err={mr['mean_error_pct']:.1f}%"),
                    zorder=4)

        # V_spin
        if gname in spin_map:
            sr     = spin_map[gname]
            v_spin = [spin_velocity(radii[i],
                                    data["v_gas"][i], data["v_disk"][i],
                                    data["v_bul"][i],
                                    sr["k"], sr["J"], sr["r0_kpc"])
                      for i in range(len(radii))]
            ax.plot(radii, v_spin, color="tab:red", linestyle="--", linewidth=1.8,
                    label=rf"Spin-med.  err={sr['mean_error_pct']:.1f}%",
                    zorder=3)

        ax.set_xlabel("Radius (kpc)", fontsize=9)
        ax.set_ylabel(r"$V$ (km s$^{-1}$)", fontsize=9)
        ax.set_title(gname, fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="lower right")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    # Hide unused panels
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.suptitle("Rotation curves: spin-mediated model vs MOND RAR",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Rotation curves    → {output_path}")


# =============================================================================
# SCATTER PLOT: MOND error vs spin-model error
# =============================================================================

def plot_mond_vs_spin_scatter(mond_results: List[Dict],
                               spin_results: Dict[str, Dict],
                               output_path: str):
    """
    Per-galaxy scatter: x = spin-model MAE%, y = MOND MAE%.
    Points coloured by k≈0 status (red = k≈0, blue = k>0).
    The cluster of mutual failures is visually compelling for the paper.
    """
    common = [r["galaxy_name"] for r in mond_results if r["galaxy_name"] in spin_results]
    if not common:
        print("  No matched galaxies for scatter plot — skipping.")
        return

    mond_map = {r["galaxy_name"]: r for r in mond_results}

    x_spin, y_mond, colors, labels = [], [], [], []
    for gname in common:
        sx = spin_results[gname]["mean_error_pct"]
        sy = mond_map[gname]["mean_error_pct"]
        is_kzero = spin_results[gname]["log10_k"] <= K_ZERO_THRESHOLD_LOG10
        x_spin.append(sx)
        y_mond.append(sy)
        colors.append("tab:red" if is_kzero else "tab:blue")
        labels.append(is_kzero)

    x_arr  = np.array(x_spin)
    y_arr  = np.array(y_mond)
    kzero  = np.array(labels)

    fig, ax = plt.subplots(figsize=(7, 6))

    ax.scatter(x_arr[~kzero], y_arr[~kzero],
               color="tab:blue", alpha=0.55, s=22, linewidths=0,
               label=f"k>0  (N={np.sum(~kzero)})")
    ax.scatter(x_arr[kzero],  y_arr[kzero],
               color="tab:red",  alpha=0.75, s=30, linewidths=0,
               label=f"k≈0  (N={np.sum(kzero)})")

    # Diagonal y=x reference
    lim = max(x_arr.max(), y_arr.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="gray", linestyle="--", linewidth=0.9,
            alpha=0.6, label="Equal errors")

    ax.set_xlabel("Spin-mediated model  MAE (%)", fontsize=12)
    ax.set_ylabel("MOND RAR  MAE (%)", fontsize=12)
    ax.set_title(
        "Per-galaxy errors: spin-mediated model vs MOND RAR\n"
        "(red = k≈0 population; dashed = equal-error line)",
        fontsize=10,
    )
    ax.legend(fontsize=9)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ MOND vs spin scatter → {output_path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    DATA_DIR       = "data/Rotmod_LTG"
    SPIN_CSV       = "results/batch_fitting_results_continuous_k.csv"
    GLOBAL_CSV     = "results/global_fit_results.csv"
    MOND_CSV       = "results/mond_comparison_results.csv"
    HISTOGRAM_PNG  = "results/mond_error_histogram.png"
    SCATTER_PNG    = "results/mond_vs_spin_scatter.png"
    ROT_CURVES_PNG = "results/mond_rotation_curves.png"
    TWO_POP_TXT    = "results/two_population_model_b.txt"

    # ------------------------------------------------------------------
    # PHASE 1: Per-galaxy MOND fit
    # ------------------------------------------------------------------
    print("\n" + "="*60)
    print("PHASE 1: PER-GALAXY MOND FIT")
    print("="*60)

    processor    = BatchProcessor(DATA_DIR)
    mond_rows    = processor.run()

    # Load spin results early so we can enrich the MOND CSV before saving.
    spin_results = load_spin_results(SPIN_CSV)

    # Enrich each MOND row with k population info for easy cross-analysis.
    for row in mond_rows:
        sr = spin_results.get(row["galaxy_name"])
        if sr:
            row["log10_k"]            = sr["log10_k"]
            row["k_zero"]             = sr["log10_k"] <= K_ZERO_THRESHOLD_LOG10
            row["spin_mean_error_pct"] = sr["mean_error_pct"]
        else:
            row["log10_k"]            = None
            row["k_zero"]             = None
            row["spin_mean_error_pct"] = None

    processor.save_csv(MOND_CSV)
    processor.print_summary()

    # ------------------------------------------------------------------
    # PHASE 2: Side-by-side comparison with spin model
    # ------------------------------------------------------------------
    print("\n" + "="*60)
    print("PHASE 2: COMPARISON WITH SPIN-MEDIATED MODEL")
    print("="*60)

    model_b_errors = load_model_b_errors(GLOBAL_CSV)

    if spin_results:
        print_comparison_summary(mond_rows, spin_results, model_b_errors or None)
    else:
        print("  Spin results unavailable — skipping comparison.")

    # ------------------------------------------------------------------
    # PHASE 3: Two-population Model B analysis
    # ------------------------------------------------------------------
    print("\n" + "="*60)
    print("PHASE 3: TWO-POPULATION MODEL B ANALYSIS")
    print("="*60)

    if spin_results and model_b_errors:
        two_population_analysis(spin_results, model_b_errors, TWO_POP_TXT)
    else:
        print("  Missing spin results or global fit CSV — skipping.")

    # ------------------------------------------------------------------
    # PHASE 4: Plots
    # ------------------------------------------------------------------
    print("\n" + "="*60)
    print("PHASE 4: PLOTS")
    print("="*60)

    if spin_results:
        plot_error_histogram(mond_rows, spin_results, HISTOGRAM_PNG)
        plot_mond_vs_spin_scatter(mond_rows, spin_results, SCATTER_PNG)

    mond_map = {r["galaxy_name"]: r for r in mond_rows}
    plot_rotation_curves(PLOT_GALAXIES, DATA_DIR, mond_map,
                          spin_results, ROT_CURVES_PNG)

    print("\n✓ g_mond_comparison.py complete.")
    print("="*60)
