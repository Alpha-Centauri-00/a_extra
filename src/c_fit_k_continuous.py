"""
Continuous k Parameter Optimization using J_visible
====================================================

UPDATED (Tasks 1 & 2):
  - Weighted chi-squared using SPARC error bars (col2 = v_err)
  - Global simultaneous fit for k0 and beta across all galaxies at once
  - Model A vs B comparison:
      Model A: single universal k (beta=0)
      Model B: power-law k(J) = k0 * J^beta

Output CSV: results/batch_fitting_results_continuous_k.csv
Global fit: results/global_fit_results.csv
"""

import os
import math
import csv
from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy.optimize import minimize
from tools.config import constants


# =============================================================================
# DATA LOADING
# =============================================================================

class RotationCurveLoader:
    """Load and parse rotation curve .dat files (SPARC format)."""

    @staticmethod
    def load(filepath: str) -> Optional[Dict]:
        result = {
            "filepath": filepath,
            "distance_mpc": None,
            "radii": [],
            "v_obs": [],
            "v_err": [],       # NEW: col2 = velocity error bar
            "v_gas": [],
            "v_disk": [],
            "v_bul": [],
        }

        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if line.startswith("# Distance"):
                    parts = line.split("=")
                    try:
                        result["distance_mpc"] = float(parts[1].strip().split()[0])
                    except Exception:
                        pass
                    break

            for line in lines:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue

                try:
                    parts = line.split()
                    if len(parts) < 3:
                        continue

                    r      = float(parts[0])
                    v_obs  = float(parts[1])
                    v_err  = float(parts[2]) if len(parts) > 2 else 0.0   # col2
                    v_gas  = float(parts[3]) if len(parts) > 3 else 0.0
                    v_disk = float(parts[4]) if len(parts) > 4 else 0.0
                    v_bul  = float(parts[5]) if len(parts) > 5 else 0.0

                    # Guard against zero / negative errors (use floor = 1 km/s)
                    v_err = max(v_err, 1.0)

                    result["radii"].append(r)
                    result["v_obs"].append(v_obs)
                    result["v_err"].append(v_err)
                    result["v_gas"].append(v_gas)
                    result["v_disk"].append(v_disk)
                    result["v_bul"].append(v_bul)

                except (ValueError, IndexError):
                    continue

        except FileNotFoundError:
            return None

        return result if result["radii"] else None


# =============================================================================
# MODEL
# =============================================================================

class RotationCurveModel:
    """Spin-mediated acceleration model."""

    @staticmethod
    def v_vis(v_gas: float, v_disk: float, v_bul: float) -> float:
        return math.sqrt(v_gas**2 + v_disk**2 + v_bul**2)

    @staticmethod
    def v_model(v_vis_kms: float, J: float, r_m: float, k: float,
                r0_m: float = constants.R0_M) -> float:
        """
        v_model = sqrt(v_vis^2 + r * k * J / (r + r0)^2)
        Returns velocity in km/s.
        """
        a_extra = k * J / ((r_m + r0_m) ** 2)
        v_vis_ms = v_vis_kms * 1000.0
        v_sq = v_vis_ms**2 + r_m * a_extra
        if v_sq < 0:
            return v_vis_kms
        return math.sqrt(v_sq) / 1000.0


# =============================================================================
# WEIGHTED CHI-SQUARED ERROR (replaces unweighted MAE)
# =============================================================================

def weighted_chi2(k: float, J: float, galaxy_data: Dict,
                  r0_m: float = constants.R0_M) -> float:
    """
    Reduced chi-squared using SPARC error bars:

        chi2 = (1/N) * sum[ (v_model - v_obs)^2 / v_err^2 ]

    Using reduced (divided by N) keeps the scale comparable across galaxies
    of different sizes and is minimised by the same k as the full chi2.
    """
    total = 0.0
    n = 0

    for i, r_kpc in enumerate(galaxy_data["radii"]):
        r_m = r_kpc * constants.KPC_TO_M
        v_vis = RotationCurveModel.v_vis(
            galaxy_data["v_gas"][i],
            galaxy_data["v_disk"][i],
            galaxy_data["v_bul"][i],
        )
        if v_vis == 0:
            continue

        v_pred = RotationCurveModel.v_model(v_vis, J, r_m, k, r0_m)
        v_obs  = galaxy_data["v_obs"][i]
        v_err  = galaxy_data["v_err"][i]

        if v_obs > 0 and v_err > 0:
            total += ((v_pred - v_obs) / v_err) ** 2
            n += 1

    return total / n if n > 0 else 1e10


def mean_abs_relative_error(k: float, J: float, galaxy_data: Dict,
                             r0_m: float = constants.R0_M) -> float:
    """Unweighted MAE for reporting purposes only (not used in optimisation)."""
    errors = []
    for i, r_kpc in enumerate(galaxy_data["radii"]):
        r_m = r_kpc * constants.KPC_TO_M
        v_vis = RotationCurveModel.v_vis(
            galaxy_data["v_gas"][i],
            galaxy_data["v_disk"][i],
            galaxy_data["v_bul"][i],
        )
        if v_vis == 0:
            continue
        v_pred = RotationCurveModel.v_model(v_vis, J, r_m, k, r0_m)
        v_obs  = galaxy_data["v_obs"][i]
        if v_obs > 0:
            errors.append(abs(v_pred - v_obs) / v_obs)
    return float(np.mean(errors)) if errors else 1e10


# =============================================================================
# PER-GALAXY FITTER  (unchanged logic, switched to weighted chi2)
# =============================================================================

class ContinuousKFitter:
    """Fit k parameter per galaxy using continuous optimisation (weighted chi2)."""

    GRID = [1e-45, 1e-43, 1e-42, 1e-41, 1e-40, 1e-38, 1e-36]
    BOUNDS = [(-55, -25)]

    def __init__(self, r0_m: float = constants.R0_M):
        self.r0_m = r0_m

    def _objective(self, log_k_arr, J, galaxy_data):
        return weighted_chi2(10 ** log_k_arr[0], J, galaxy_data, self.r0_m)

    def fit(self, J: float, galaxy_data: Dict) -> Tuple[float, float, float]:
        """
        Returns:
            (k_optimal, chi2_reduced, mean_abs_relative_error_pct)
        """
        best_log_k = min(
            (math.log10(k) for k in self.GRID),
            key=lambda lk: weighted_chi2(10**lk, J, galaxy_data, self.r0_m),
        )

        result = minimize(
            self._objective,
            x0=np.array([best_log_k]),
            args=(J, galaxy_data),
            method="L-BFGS-B",
            bounds=self.BOUNDS,
            options={"ftol": 1e-12, "maxiter": 5000},
        )

        k_opt = 10 ** result.x[0]
        chi2  = result.fun
        mare  = mean_abs_relative_error(k_opt, J, galaxy_data, self.r0_m)

        return k_opt, chi2, mare


# =============================================================================
# GLOBAL SIMULTANEOUS FIT  (Task 1, part 2 + Task 2)
# =============================================================================

class GlobalFitter:
    """
    Fit k0 and beta simultaneously across ALL galaxies using:

        k(J) = k0 * J^beta      (Model B)

    or

        k(J) = k0               (Model A, beta forced = 0)

    Minimises the sum of per-galaxy reduced chi-squared values.
    """

    # Upper bound on log10(k0) raised to +30: when beta < 0 and J ~ 10^65,
    # the effective k = k0 * J^beta can still be physically reasonable even
    # with large k0, because J^beta suppresses it. The physical ceiling
    # log10_k > 0 is enforced inside _total_chi2, not here.
    BOUNDS_LOG_K0 = (-80, +30)
    BOUNDS_BETA   = (-5.0, 5.0)

    def __init__(self, r0_m: float = constants.R0_M):
        self.r0_m = r0_m

    # ------------------------------------------------------------------
    # Internal: total objective over all galaxies
    # ------------------------------------------------------------------
    def _total_chi2(self, params, J_list, data_list, fix_beta=False,
                    r0_m_list=None):
        if fix_beta:
            log_k0 = params[0]
            beta   = 0.0
        else:
            log_k0, beta = params[0], params[1]

        total = 0.0

        for i, (J, data) in enumerate(zip(J_list, data_list)):
            r0_m = r0_m_list[i] if r0_m_list is not None else self.r0_m
            # Compute k entirely in log-space to prevent overflow/underflow
            # when J ~ 10^62 and beta is non-zero.
            if fix_beta or beta == 0.0:
                log10_k = log_k0
            else:
                log10_J = math.log10(J) if J > 0 else 0.0
                log10_k = log_k0 + beta * log10_J
            # log10_k > 0 means k > 1 SI, which is unphysical and causes
            # float overflow (10^308 limit) when J ~ 10^66 and beta ~ 5.
            if not math.isfinite(log10_k) or log10_k > 0:
                total += 1e10
                continue
            k = 10 ** log10_k
            total += weighted_chi2(k, J, data, r0_m)

        return total / max(len(J_list), 1)

    # ------------------------------------------------------------------
    # Helper: warm start from per-galaxy sample
    # ------------------------------------------------------------------
    def _warm_start_log_k(self, J_list, data_list, r0_m_list=None):
        """
        Fit k on a random sample of 20 galaxies and return the median
        log10(k).  This gives the global optimizer a sensible starting
        point and prevents it from getting stuck at the search boundary.
        """
        n    = min(20, len(J_list))
        idxs = np.random.choice(len(J_list), n, replace=False)
        log_ks = []
        for i in idxs:
            try:
                r0_m  = r0_m_list[i] if r0_m_list is not None else self.r0_m
                k_opt, _, _ = ContinuousKFitter(r0_m=r0_m).fit(J_list[i], data_list[i])
                if k_opt > 0:
                    log_ks.append(math.log10(k_opt))
            except Exception:
                continue
        return float(np.median(log_ks)) if log_ks else -42.0

    # ------------------------------------------------------------------
    # Model A: single universal k  (beta = 0)
    # ------------------------------------------------------------------
    def fit_model_A(self, J_list: List[float],
                    data_list: List[Dict],
                    r0_m_list: Optional[List[float]] = None) -> Dict:
        """
        Model A: k(J) = k0   (one number for all galaxies).
        Uses a warm start from the median of per-galaxy fits.
        r0_m_list: per-galaxy softening radii (m). Falls back to self.r0_m if None.
        """
        print("\n  [Model A] Fitting universal k (beta = 0) …")
        print("    Computing warm start from per-galaxy sample …")

        lk0_start = self._warm_start_log_k(J_list, data_list, r0_m_list)
        print(f"    Warm start: log10(k0) = {lk0_start:.2f}")

        # Fine grid around warm start
        lo   = max(lk0_start - 6, self.BOUNDS_LOG_K0[0])
        hi   = min(lk0_start + 6, self.BOUNDS_LOG_K0[1])
        grid = np.linspace(lo, hi, 25)
        best_lk0 = min(grid,
                       key=lambda lk: self._total_chi2([lk], J_list,
                                                        data_list, True, r0_m_list))

        result = minimize(
            self._total_chi2,
            x0=np.array([best_lk0]),
            args=(J_list, data_list, True, r0_m_list),
            method="L-BFGS-B",
            bounds=[self.BOUNDS_LOG_K0],
            options={"ftol": 1e-14, "gtol": 1e-8, "maxiter": 10000},
        )

        k0_A   = 10 ** result.x[0]
        chi2_A = result.fun
        print(f"    k0 = {k0_A:.4e}  |  chi2_reduced = {chi2_A:.4f}")

        per_galaxy_errors = []
        for i, (J, data) in enumerate(zip(J_list, data_list)):
            r0_m = r0_m_list[i] if r0_m_list is not None else self.r0_m
            mare = mean_abs_relative_error(k0_A, J, data, r0_m)
            per_galaxy_errors.append(mare * 100)

        return {
            "model": "A",
            "description": "Universal k (beta=0)",
            "k0": k0_A,
            "log10_k0": math.log10(k0_A),
            "beta": 0.0,
            "chi2_reduced": chi2_A,
            "mean_error_pct": float(np.mean(per_galaxy_errors)),
            "median_error_pct": float(np.median(per_galaxy_errors)),
            "per_galaxy_errors": per_galaxy_errors,
        }

    # ------------------------------------------------------------------
    # Model B: power-law k(J) = k0 * J^beta
    # ------------------------------------------------------------------
    def fit_model_B(self, J_list: List[float],
                    data_list: List[Dict],
                    r0_m_list: Optional[List[float]] = None) -> Dict:
        """
        Model B: k(J) = k0 * J^beta   (two parameters)

        Strategy:
          1. Warm start: median k from per-galaxy sample, beta = -1.0
          2. 2D grid around warm start to find best neighbourhood
          3. L-BFGS-B refinement
        r0_m_list: per-galaxy softening radii (m). Falls back to self.r0_m if None.
        """
        print("\n  [Model B] Fitting power-law k(J) = k0 * J^beta …")
        print("    Computing warm start …")

        lk0_start = self._warm_start_log_k(J_list, data_list, r0_m_list)
        print(f"    Warm start: log10(k0) = {lk0_start:.2f}, beta = -1.0")

        # 2D grid around warm start
        lk0_grid  = np.linspace(max(lk0_start - 8, self.BOUNDS_LOG_K0[0]),
                                min(lk0_start + 8, self.BOUNDS_LOG_K0[1]), 10)
        beta_grid = np.linspace(-5.0, 0.5, 12)  # extended: -3.0 caused boundary hits
        best_val  = np.inf
        best_x0   = [lk0_start, -1.0]

        for lk0 in lk0_grid:
            for b in beta_grid:
                val = self._total_chi2([lk0, b], J_list, data_list, False, r0_m_list)
                if val < best_val:
                    best_val = val
                    best_x0  = [lk0, b]

        print(f"    Grid best: log10(k0)={best_x0[0]:.1f}, beta={best_x0[1]:.2f},"
              f" chi2={best_val:.4f}")

        result = minimize(
            self._total_chi2,
            x0=np.array(best_x0),
            args=(J_list, data_list, False, r0_m_list),
            method="L-BFGS-B",
            bounds=[self.BOUNDS_LOG_K0, self.BOUNDS_BETA],
            options={"ftol": 1e-14, "gtol": 1e-8, "maxiter": 10000},
        )

        log_k0_B, beta_B = result.x[0], result.x[1]
        k0_B             = 10 ** log_k0_B
        chi2_B           = result.fun

        print(f"    k0 = {k0_B:.4e}  |  beta = {beta_B:.4f}"
              f"  |  chi2_reduced = {chi2_B:.4f}")

        # Per-galaxy errors (log-space k to avoid overflow)
        per_galaxy_errors = []
        for i, (J, data) in enumerate(zip(J_list, data_list)):
            r0_m    = r0_m_list[i] if r0_m_list is not None else self.r0_m
            log10_J = math.log10(J) if J > 0 else 0.0
            log10_k = log_k0_B + (beta_B * log10_J if beta_B != 0.0 else 0.0)
            k       = max(10 ** log10_k, 1e-60)
            mare    = mean_abs_relative_error(k, J, data, r0_m)
            per_galaxy_errors.append(mare * 100)

        return {
            "model": "B",
            "description": "Power-law k(J) = k0 * J^beta",
            "k0": k0_B,
            "log10_k0": log_k0_B,
            "beta": beta_B,
            "chi2_reduced": chi2_B,
            "mean_error_pct": float(np.mean(per_galaxy_errors)),
            "median_error_pct": float(np.median(per_galaxy_errors)),
            "per_galaxy_errors": per_galaxy_errors,
        }

    # ------------------------------------------------------------------
    # Model comparison (AIC / likelihood ratio)
    # ------------------------------------------------------------------
    @staticmethod
    def compare_models(result_A: Dict, result_B: Dict, n_total_datapoints: int) -> Dict:
        """
        Delta-chi2 comparison and AIC.

        AIC = 2*k_params - 2*ln(L)
        For Gaussian errors, -2ln(L) ∝ N * chi2_reduced, so:
            AIC_A = 2*1 + N * chi2_A
            AIC_B = 2*2 + N * chi2_B
        where N = total number of data points across all galaxies.
        Delta_AIC = AIC_B - AIC_A
          < 0 → Model B preferred
          > 0 → Model A preferred (simpler)
        """
        N = n_total_datapoints  # total data points, not number of galaxies
        aic_A = 2 * 1 + N * result_A["chi2_reduced"]
        aic_B = 2 * 2 + N * result_B["chi2_reduced"]
        delta_aic = aic_B - aic_A

        delta_chi2 = result_A["chi2_reduced"] - result_B["chi2_reduced"]

        return {
            "delta_chi2_A_minus_B": delta_chi2,
            "AIC_A": aic_A,
            "AIC_B": aic_B,
            "delta_AIC_B_minus_A": delta_aic,
            "preferred": "B (power-law)" if delta_aic < 0 else "A (universal k)",
        }


# =============================================================================
# J_VISIBLE LOADER
# =============================================================================

class J_VisibleLoader:

    @staticmethod
    def load_csv(csv_path: str) -> Dict[str, Dict]:
        """Return {galaxy_name: {'J': float, 'r_disk_kpc': float}}."""
        data = {}
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        data[row['galaxy_name']] = {
                            'J':          float(row['J_new_visible']),
                            'r_disk_kpc': float(row['r_disk_kpc']),
                        }
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(f"Warning: J_visible CSV not found at {csv_path}")
        return data


# =============================================================================
# BATCH PROCESSOR  (per-galaxy fits)
# =============================================================================

class BatchProcessor:

    # r0 = R0_ALPHA * r_disk_kpc  (robustness scan optimum: alpha=2.0)
    R0_ALPHA = 2.0

    def __init__(self, j_visible_csv: str):
        self.j_visible = J_VisibleLoader.load_csv(j_visible_csv)
        self.results   = []

    def find_files(self, directory: str) -> List[Tuple[str, str]]:
        files = []
        if not os.path.exists(directory):
            print(f"Directory not found: {directory}")
            return files
        for filename in sorted(os.listdir(directory)):
            if filename.endswith("_rotmod.dat"):
                files.append((
                    os.path.join(directory, filename),
                    filename.replace("_rotmod.dat", ""),
                ))
        return files

    def process_galaxy(self, filepath: str, galaxy_name: str) -> Optional[Dict]:
        print(f"  {galaxy_name} …", end=" ")
        data = RotationCurveLoader.load(filepath)
        if data is None:
            print("LOAD FAILED"); return None

        if galaxy_name not in self.j_visible:
            print("J not found"); return None

        gal_info   = self.j_visible[galaxy_name]
        J_new      = gal_info['J']
        r_disk_kpc = gal_info['r_disk_kpc']
        r0_kpc     = self.R0_ALPHA * r_disk_kpc
        r0_m       = r0_kpc * constants.KPC_TO_M

        fitter = ContinuousKFitter(r0_m=r0_m)
        k_opt, chi2, mare = fitter.fit(J_new, data)

        result = {
            "galaxy_name":       galaxy_name,
            "distance_mpc":      data["distance_mpc"],
            "r_max_kpc":         data["radii"][-1],
            "v_max_km_s":        data["v_obs"][-1],
            "J_new":             J_new,
            "log10_J_new":       math.log10(J_new) if J_new > 0 else 0.0,
            "r0_kpc":            r0_kpc,
            "k_optimal":         k_opt,
            "log10_k_optimal":   math.log10(k_opt) if k_opt > 0 else 0.0,
            "chi2_reduced":      chi2,
            "mean_error":        mare,
            "mean_error_pct":    mare * 100,
        }
        self.results.append(result)
        print(f"r0={r0_kpc:.2f}kpc  k={k_opt:.3e}  chi2={chi2:.3f}  err={mare*100:.2f}%")
        return result

    def process_all(self, directory: str):
        files = self.find_files(directory)
        if not files:
            return []
        print(f"\n  {len(files)} galaxies found, "
              f"{len(self.j_visible)} have J_new\n")
        for fp, name in files:
            self.process_galaxy(fp, name)
        return self.results

    def save_csv(self, output_path: str):
        if not self.results:
            return
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                "galaxy_name", "distance_mpc", "r_max_kpc", "v_max_km_s",
                "J_new", "log10_J_new", "r0_kpc",
                "k_optimal", "log10_k_optimal",
                "chi2_reduced", "mean_error", "mean_error_pct",
            ])
            writer.writeheader()
            writer.writerows(self.results)
        print(f"\n  ✓ Per-galaxy results → {output_path}")

    def print_summary(self):
        if not self.results:
            return
        errs  = [r['mean_error_pct'] for r in self.results]
        chi2s = [r['chi2_reduced']    for r in self.results]
        print("\n" + "="*80)
        print("PER-GALAXY FIT SUMMARY  (weighted chi2 optimisation)")
        print("="*80)
        print(f"  Galaxies          : {len(self.results)}")
        print(f"  Mean   error      : {np.mean(errs):.2f}% ± {np.std(errs):.2f}%")
        print(f"  Median error      : {np.median(errs):.2f}%")
        print(f"  Best / worst      : {min(errs):.2f}% / {max(errs):.2f}%")
        print(f"  Mean   chi2_red   : {np.mean(chi2s):.3f}")
        print(f"  Median chi2_red   : {np.median(chi2s):.3f}")
        print("="*80)


# =============================================================================
# SAVE GLOBAL FIT RESULTS
# =============================================================================

def save_global_results(result_A: Dict, result_B: Dict,
                        comparison: Dict, output_path: str,
                        galaxy_names: List[str]):
    """Save global fit summary + per-galaxy predicted errors to CSV."""

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Summary file
    summary_path = output_path.replace(".csv", "_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("GLOBAL FIT RESULTS\n")
        f.write("="*60 + "\n\n")

        for res in [result_A, result_B]:
            f.write(f"Model {res['model']}: {res['description']}\n")
            f.write(f"  k0              = {res['k0']:.6e}\n")
            f.write(f"  log10(k0)       = {res['log10_k0']:.4f}\n")
            f.write(f"  beta            = {res['beta']:.4f}\n")
            f.write(f"  chi2_reduced    = {res['chi2_reduced']:.4f}\n")
            f.write(f"  Mean  error     = {res['mean_error_pct']:.2f}%\n")
            f.write(f"  Median error    = {res['median_error_pct']:.2f}%\n\n")

        f.write("MODEL COMPARISON\n")
        f.write(f"  Delta chi2 (A-B)  = {comparison['delta_chi2_A_minus_B']:+.4f}\n")
        f.write(f"  AIC_A             = {comparison['AIC_A']:.2f}\n")
        f.write(f"  AIC_B             = {comparison['AIC_B']:.2f}\n")
        f.write(f"  Delta AIC (B-A)   = {comparison['delta_AIC_B_minus_A']:+.2f}\n")
        f.write(f"  Preferred model   : {comparison['preferred']}\n")

    print(f"  ✓ Global summary  → {summary_path}")

    # Per-galaxy CSV
    rows = []
    for i, name in enumerate(galaxy_names):
        rows.append({
            "galaxy_name":       name,
            "error_pct_modelA":  f"{result_A['per_galaxy_errors'][i]:.4f}",
            "error_pct_modelB":  f"{result_B['per_galaxy_errors'][i]:.4f}",
        })

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f,
                    fieldnames=["galaxy_name", "error_pct_modelA", "error_pct_modelB"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"  ✓ Global per-galaxy → {output_path}")


# =============================================================================
# MAIN
# =============================================================================

# =============================================================================
# TASK 4: VARIABLE r0 ROBUSTNESS SCAN
# =============================================================================

class R0RobustnessTester:
    """
    Task 4: Test whether fixing r0 = 0.5 kpc matters.

    Strategies tested:
      (a) r0 = alpha * h         (scales with disk scale length)
      (b) r0 = alpha * r_half    (scales with half-mass radius from b_compute_J)
      (c) r0 free per galaxy     (full freedom — upper bound on improvement)

    For each strategy we report mean chi2 and mean error across all galaxies.
    The goal is to see if the global chi2 improves significantly over the
    fixed-r0 baseline.  If not, r0 = 0.5 kpc is robust.
    """

    # Multiplier candidates for strategies (a) and (b)
    ALPHA_GRID = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

    def __init__(self):
        self.fitter_fixed = ContinuousKFitter(r0_m=constants.R0_M)

    # ------------------------------------------------------------------
    def _fit_with_r0(self, J: float, data: Dict,
                     r0_m: float) -> Tuple[float, float]:
        """Fit k with a specific r0 and return (chi2_reduced, mare_pct)."""
        fitter = ContinuousKFitter(r0_m=r0_m)
        k_opt, chi2, mare = fitter.fit(J, data)
        return chi2, mare * 100

    # ------------------------------------------------------------------
    def scan_fixed_r0(self, J_list, data_list,
                      r0_kpc_values=None) -> List[Dict]:
        """
        Scan a grid of fixed r0 values and report aggregate chi2.
        Baseline for comparison with adaptive strategies.
        """
        if r0_kpc_values is None:
            r0_kpc_values = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]

        results = []
        for r0_kpc in r0_kpc_values:
            r0_m   = r0_kpc * constants.KPC_TO_M
            chi2s  = []
            mares  = []
            for J, data in zip(J_list, data_list):
                chi2, mare = self._fit_with_r0(J, data, r0_m)
                chi2s.append(chi2)
                mares.append(mare)

            results.append({
                "strategy":          f"fixed r0 = {r0_kpc} kpc",
                "r0_kpc":            r0_kpc,
                "mean_chi2":         float(np.mean(chi2s)),
                "median_chi2":       float(np.median(chi2s)),
                "mean_error_pct":    float(np.mean(mares)),
                "median_error_pct":  float(np.median(mares)),
            })
            print(f"  r0={r0_kpc:5.2f} kpc | mean chi2={np.mean(chi2s):.4f}"
                  f" | mean err={np.mean(mares):.2f}%")

        return results

    # ------------------------------------------------------------------
    def scan_h_scaling(self, J_list, data_list,
                       h_kpc_list: List[Optional[float]]) -> List[Dict]:
        """
        Strategy (a): r0 = alpha * h
        Requires h_kpc for each galaxy (from morphology CSV).
        Galaxies with h = None use the default 0.5 kpc.
        """
        results = []

        for alpha in self.ALPHA_GRID:
            chi2s = []
            mares = []
            n_used = 0

            for J, data, h in zip(J_list, data_list, h_kpc_list):
                if h is not None and h > 0:
                    r0_m = alpha * h * constants.KPC_TO_M
                    n_used += 1
                else:
                    r0_m = constants.R0_M   # fallback
                chi2, mare = self._fit_with_r0(J, data, r0_m)
                chi2s.append(chi2)
                mares.append(mare)

            results.append({
                "strategy":          f"r0 = {alpha} * h",
                "alpha":             alpha,
                "n_with_h":          n_used,
                "mean_chi2":         float(np.mean(chi2s)),
                "median_chi2":       float(np.median(chi2s)),
                "mean_error_pct":    float(np.mean(mares)),
                "median_error_pct":  float(np.median(mares)),
            })
            print(f"  alpha={alpha:.2f} * h  | mean chi2={np.mean(chi2s):.4f}"
                  f" | mean err={np.mean(mares):.2f}%  (N_h={n_used})")

        return results

    # ------------------------------------------------------------------
    def scan_r_half_scaling(self, J_list, data_list,
                            r_half_kpc_list: List[Optional[float]]) -> List[Dict]:
        """
        Strategy (b): r0 = alpha * r_half_mass
        r_half comes from b_compute_J_visible (r_disk_kpc column).
        """
        results = []

        for alpha in self.ALPHA_GRID:
            chi2s = []
            mares = []
            n_used = 0

            for J, data, r_half in zip(J_list, data_list, r_half_kpc_list):
                if r_half is not None and r_half > 0:
                    r0_m = alpha * r_half * constants.KPC_TO_M
                    n_used += 1
                else:
                    r0_m = constants.R0_M
                chi2, mare = self._fit_with_r0(J, data, r0_m)
                chi2s.append(chi2)
                mares.append(mare)

            results.append({
                "strategy":          f"r0 = {alpha} * r_half",
                "alpha":             alpha,
                "n_with_r_half":     n_used,
                "mean_chi2":         float(np.mean(chi2s)),
                "median_chi2":       float(np.median(chi2s)),
                "mean_error_pct":    float(np.mean(mares)),
                "median_error_pct":  float(np.median(mares)),
            })
            print(f"  alpha={alpha:.2f} * r_half | mean chi2={np.mean(chi2s):.4f}"
                  f" | mean err={np.mean(mares):.2f}%  (N_r_half={n_used})")

        return results

    # ------------------------------------------------------------------
    @staticmethod
    def print_r0_comparison(fixed_results, h_results, r_half_results):
        """Print summary table comparing all r0 strategies."""
        print("\n" + "="*80)
        print("TASK 4: r0 ROBUSTNESS COMPARISON")
        print("="*80)
        print(f"\n{'Strategy':<28} {'Mean chi2':>10} {'Mean err (%)':>14}")
        print("-"*54)

        baseline = None
        for r in fixed_results:
            marker = " ← baseline" if abs(r['r0_kpc'] - 0.5) < 0.01 else ""
            if marker:
                baseline = r['mean_chi2']
            print(f"  {r['strategy']:<26} {r['mean_chi2']:>10.4f}"
                  f"  {r['mean_error_pct']:>10.2f}%{marker}")

        print()
        best_h = min(h_results, key=lambda x: x['mean_chi2'])
        print(f"  Best r0 = alpha*h:   alpha={best_h['alpha']:.2f}"
              f"  chi2={best_h['mean_chi2']:.4f}"
              f"  err={best_h['mean_error_pct']:.2f}%")

        best_rh = min(r_half_results, key=lambda x: x['mean_chi2'])
        print(f"  Best r0 = alpha*r½:  alpha={best_rh['alpha']:.2f}"
              f"  chi2={best_rh['mean_chi2']:.4f}"
              f"  err={best_rh['mean_error_pct']:.2f}%")

        if baseline:
            delta_h  = best_h['mean_chi2']  - baseline
            delta_rh = best_rh['mean_chi2'] - baseline
            print(f"\n  Δchi2 (best h-scaling  vs fixed 0.5 kpc): {delta_h:+.4f}")
            print(f"  Δchi2 (best r½-scaling vs fixed 0.5 kpc): {delta_rh:+.4f}")
            print(f"\n  Interpretation: Δchi2 < 0.01 → r0 = 0.5 kpc is robust.")

        print("="*80)

    # ------------------------------------------------------------------
    def save_r0_scan_csv(self, all_results: List[Dict], output_path: str):
        if not all_results:
            return
        fieldnames = list(all_results[0].keys())
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames,
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_results)
        print(f"  ✓ r0 scan results → {output_path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    J_VISIBLE_CSV   = "results/J_visible_from_data.csv"
    DATA_DIR        = "data/Rotmod_LTG"
    MORPHOLOGY_CSV  = "results/galaxy_properties_h_bt.csv"   # from script a
    PER_GAL_CSV     = "results/batch_fitting_results_continuous_k.csv"
    GLOBAL_CSV      = "results/global_fit_results.csv"
    R0_SCAN_CSV     = "results/r0_robustness_scan.csv"

    # ------------------------------------------------------------------
    # PHASE 1: Per-galaxy k optimisation (weighted chi2)
    # ------------------------------------------------------------------
    print("\n" + "="*80)
    print("PHASE 1: PER-GALAXY k OPTIMISATION  (weighted chi2)")
    print("="*80)

    batch = BatchProcessor(J_VISIBLE_CSV)
    batch.process_all(DATA_DIR)
    batch.save_csv(PER_GAL_CSV)
    batch.print_summary()

    # ------------------------------------------------------------------
    # PHASE 2: Global simultaneous fit  (Tasks 1 + 2)
    # ------------------------------------------------------------------
    print("\n" + "="*80)
    print("PHASE 2: GLOBAL SIMULTANEOUS FIT  (Model A vs B)")
    print("="*80)

    J_list     = None
    data_list  = None
    names_list = None

    if len(batch.results) < 5:
        print("  Not enough galaxies for global fit — skipping.")
    else:
        J_list     = [r["J_new"]        for r in batch.results]
        names_list = [r["galaxy_name"]  for r in batch.results]
        r0_m_list  = [r["r0_kpc"] * constants.KPC_TO_M for r in batch.results]

        # Re-load data for each galaxy (needed by GlobalFitter + R0Tester)
        data_list = []
        for r in batch.results:
            fp = os.path.join(DATA_DIR, f"{r['galaxy_name']}_rotmod.dat")
            d  = RotationCurveLoader.load(fp)
            data_list.append(d)

        gfitter    = GlobalFitter()
        result_A   = gfitter.fit_model_A(J_list, data_list, r0_m_list)
        result_B   = gfitter.fit_model_B(J_list, data_list, r0_m_list)
        n_total_pts = sum(len(d["radii"]) for d in data_list if d is not None)
        comparison = GlobalFitter.compare_models(result_A, result_B, n_total_pts)

        print("\n" + "="*80)
        print("MODEL COMPARISON SUMMARY")
        print("="*80)
        print(f"  {'Metric':<30} {'Model A':>14} {'Model B':>14}")
        print("  " + "-"*58)
        print(f"  {'k0':<30} {result_A['k0']:>14.4e} {result_B['k0']:>14.4e}")
        print(f"  {'beta':<30} {result_A['beta']:>14.4f} {result_B['beta']:>14.4f}")
        print(f"  {'chi2_reduced':<30} {result_A['chi2_reduced']:>14.4f}"
              f" {result_B['chi2_reduced']:>14.4f}")
        print(f"  {'Mean error (%)':<30} {result_A['mean_error_pct']:>14.2f}"
              f" {result_B['mean_error_pct']:>14.2f}")
        print(f"  {'Median error (%)':<30} {result_A['median_error_pct']:>14.2f}"
              f" {result_B['median_error_pct']:>14.2f}")
        print(f"  {'AIC':<30} {comparison['AIC_A']:>14.2f}"
              f" {comparison['AIC_B']:>14.2f}")
        print(f"  {'Delta AIC (B - A)':<30}"
              f" {comparison['delta_AIC_B_minus_A']:>14.2f}")
        print(f"\n  PREFERRED: {comparison['preferred']}")
        print("="*80)

        save_global_results(result_A, result_B, comparison,
                            GLOBAL_CSV, names_list)

    # ------------------------------------------------------------------
    # PHASE 3: r0 robustness scan  (Task 4)
    # ------------------------------------------------------------------
    print("\n" + "="*80)
    print("PHASE 3: r0 ROBUSTNESS SCAN  (Task 4)")
    print("="*80)

    if J_list is None or data_list is None:
        print("  Skipping r0 scan (no global fit data available).")
    else:
        # Load h and r_half from CSVs (if available)
        h_kpc_map      = {}
        r_half_kpc_map = {}

        if os.path.exists(MORPHOLOGY_CSV):
            with open(MORPHOLOGY_CSV, 'r') as f:
                for row in csv.DictReader(f):
                    gname = row['galaxy_name']
                    try:
                        h_kpc_map[gname] = float(row['h_kpc']) if row['h_kpc'] else None
                    except (ValueError, KeyError):
                        pass

        if os.path.exists(J_VISIBLE_CSV):
            with open(J_VISIBLE_CSV, 'r') as f:
                for row in csv.DictReader(f):
                    gname = row['galaxy_name']
                    try:
                        r_half_kpc_map[gname] = float(row['r_disk_kpc']) \
                            if row['r_disk_kpc'] else None
                    except (ValueError, KeyError):
                        pass

        h_list      = [h_kpc_map.get(n)      for n in names_list]
        r_half_list = [r_half_kpc_map.get(n) for n in names_list]

        tester = R0RobustnessTester()

        print("\n[a] Fixed r0 grid:")
        fixed_results  = tester.scan_fixed_r0(J_list, data_list)

        print("\n[b] r0 = alpha * h (disk scale length):")
        h_results      = tester.scan_h_scaling(J_list, data_list, h_list)

        print("\n[c] r0 = alpha * r_half (half-mass radius):")
        r_half_results = tester.scan_r_half_scaling(J_list, data_list, r_half_list)

        R0RobustnessTester.print_r0_comparison(
            fixed_results, h_results, r_half_results)

        all_r0_results = fixed_results + h_results + r_half_results
        tester.save_r0_scan_csv(all_r0_results, R0_SCAN_CSV)

    print("\n✓ c_fit_k_continuous.py complete.")
    print("="*80)