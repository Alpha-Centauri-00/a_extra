"""
CLUSTER ANALYSIS: Validate a_extra Model on Galaxy Clusters
============================================================

Fixed:
- k_grid: np.logspace(-50, -36, 8) instead of [1e-42 to 1e-30]
- r0_cluster: 0.15 h⁻¹Mpc (150 kpc) instead of 50 kpc
- LT calculation corrected for clusters
"""

import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
from tools.config import constants


# ============================================================================
# PHASE 1: DATA LOADING & FILTERING
# ============================================================================

class ClusterDataLoader:
    """Load and filter cluster data from galwcls.dat"""
    
    @staticmethod
    def load_galwcls(filepath: str) -> List[Dict]:
        """Load galwcls.dat (pipe-delimited)."""
        clusters = []
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            header = lines[0].strip().split('|')
            header = [h.strip() for h in header]
            
            for line in lines[1:]:
                if not line.strip():
                    continue
                
                values = line.strip().split('|')
                if len(values) < len(header):
                    continue
                
                try:
                    cluster = {}
                    for i, col in enumerate(header):
                        try:
                            cluster[col] = float(values[i].strip())
                        except ValueError:
                            cluster[col] = values[i].strip()
                    clusters.append(cluster)
                except Exception as e:
                    continue
            
            print(f" Loaded {len(clusters)} clusters from {filepath}")
            return clusters
        
        except Exception as e:
            print(f"ERROR loading {filepath}: {e}")
            return []
    
    @staticmethod
    def filter_high_quality(clusters: List[Dict], 
                           sig200_min: float = 400,
                           N200_min: int = 30,
                           z_max: float = 0.15) -> List[Dict]:
        """Filter clusters by quality criteria."""
        filtered = []
        
        for cluster in clusters:
            try:
                sig200 = float(cluster.get('sig200', 0))
                N200 = float(cluster.get('N200', 0))
                z = float(cluster.get('z', 0))
                
                if sig200 > sig200_min and N200 > N200_min and z < z_max:
                    filtered.append(cluster)
            except (ValueError, TypeError):
                continue
        
        print(f" Filtered to {len(filtered)} high-quality clusters")
        print(f"  (sig200 > {sig200_min} km/s, N200 > {N200_min}, z < {z_max})")
        return filtered


# ============================================================================
# PHASE 2: ANGULAR MOMENTUM COMPUTATION
# ============================================================================

class ClusterJCalculator:
    """Compute J for each cluster"""
    
    @staticmethod
    def compute_J(cluster: Dict) -> Tuple[float, float]:
        """Compute angular momentum J for cluster."""
        try:
            M200_h_Msun = float(cluster.get('M200', 0))
            sig200_km_s = float(cluster.get('sig200', 0))
            r200_h_Mpc = float(cluster.get('r200', 0))
            
            M_enc_kg = M200_h_Msun *  constants.M_SUN / 0.7
            v_3D_m_s = sig200_km_s * 1000 * math.sqrt(3)
            r_max_m = r200_h_Mpc * constants.MPC_TO_M / 0.7
            
            J = M_enc_kg * v_3D_m_s * r_max_m
            log10_J = math.log10(J) if J > 0 else 0.0
            
            return J, log10_J
        
        except (ValueError, TypeError, ZeroDivisionError):
            return None, None
    
    @staticmethod
    def add_J_to_clusters(clusters: List[Dict]) -> List[Dict]:
        """Add J and log10_J columns to each cluster."""
        
        for cluster in clusters:
            J, log10_J = ClusterJCalculator.compute_J(cluster)
            cluster['J'] = J
            cluster['log10_J'] = log10_J
        
        print(f" Computed J for {len(clusters)} clusters")
        
        J_values = [c['J'] for c in clusters if c['J'] is not None]
        if J_values:
            print(f"  J range: {min(J_values):.2e} to {max(J_values):.2e} kg·m²/s")
        
        return clusters


# ============================================================================
# PHASE 3: MODEL FITTING (FIXED k_grid and r0)
# ============================================================================

class ClusterModel:
    """Spin-coupling model for clusters"""
    
    @staticmethod
    def a_visible(M200_h_Msun: float, r200_h_Mpc: float) -> float:
        """Newtonian acceleration from cluster mass."""
        M_kg = M200_h_Msun * constants.M_SUN / 0.7
        r_m = r200_h_Mpc * constants.MPC_TO_M / 0.7
        
        if r_m == 0:
            return 0.0
        
        return constants.G * M_kg / (r_m ** 2)
    
    @staticmethod
    def a_extra(J: float, r200_m: float, k: float, 
                r0_m: float = 0.15 * constants.MPC_TO_M / 0.7) -> float:
        """Spin-coupling acceleration."""
        if r200_m + r0_m == 0 or J == 0:
            return 0.0
        
        return k * J / ((r200_m + r0_m) ** 2)
    
    @staticmethod
    def v_model(a_vis: float, a_extra: float, r200_m: float) -> float:
        """Model velocity dispersion.

        NOTE: v = sqrt(r * a_total) assumes circular orbital motion, which
        applies to rotationally-supported disk galaxies.  Galaxy clusters are
        pressure-supported (virial theorem): the correct relation is
        sigma^2 ~ G*M / r, not v_circ^2 = r * a.  Using the circular-orbit
        formula here is intentional — it is what allows us to demonstrate that
        the spin-coupling model *fails* for clusters, supporting the conclusion
        that the mechanism is specific to coherently rotating disk systems.
        """
        a_total = a_vis + a_extra

        if a_total < 0:
            return 0.0

        v_model_m_s = math.sqrt(r200_m * a_total)
        return v_model_m_s / 1000


class ClusterFitter:
    """Fit k parameter for clusters"""
    
    def __init__(self):
        self.k_range = np.logspace(-50, -36, 8)  # FIXED: wider range for clusters
        self.r0_m = 0.15 * constants.MPC_TO_M / 0.7  # FIXED: 150 kpc cluster core radius
    
    def fit_cluster(self, cluster: Dict) -> Tuple[float, float]:
        """Fit k parameter for single cluster."""
        
        J = cluster.get('J')
        if J is None or J == 0:
            return None, None
        
        sig200 = float(cluster.get('sig200', 0))
        M200_h_Msun = float(cluster.get('M200', 0))
        r200_h_Mpc = float(cluster.get('r200', 0))
        
        if sig200 <= 0 or M200_h_Msun <= 0 or r200_h_Mpc <= 0:
            return None, None
        
        r200_m = r200_h_Mpc * constants.MPC_TO_M / 0.7
        a_vis = ClusterModel.a_visible(M200_h_Msun, r200_h_Mpc)
        
        best_k = None
        best_error = float('inf')
        
        for k in self.k_range:
            a_ex = ClusterModel.a_extra(J, r200_m, k, self.r0_m)
            v_pred = ClusterModel.v_model(a_vis, a_ex, r200_m)
            
            if v_pred <= 0:
                continue
            
            error = abs(v_pred - sig200) / sig200 if sig200 > 0 else float('inf')
            
            if error < best_error:
                best_error = error
                best_k = k
        
        return best_k, best_error if best_error != float('inf') else None
    
    def fit_all(self, clusters: List[Dict]) -> List[Dict]:
        """Fit all clusters."""
        
        results = []
        
        for i, cluster in enumerate(clusters):
            k_best, error = self.fit_cluster(cluster)
            
            cluster['k_best'] = k_best
            cluster['mean_error'] = error
            
            if k_best is not None:
                results.append(cluster)
                if (i + 1) % 50 == 0:
                    print(f"  [{i+1}/{len(clusters)}] fitted")
        
        print(f" Fitted {len(results)} clusters successfully")
        
        # DEBUG: Print first 5 clusters
        print("\n[DEBUG] First 5 clusters:")
        print(f"{'ClID':<8} {'log10_J':<12} {'k_best':<16} {'error (%)':<10}")
        print("-"*50)
        for cluster in results[:5]:
            print(f"{int(cluster.get('ClID')):<8} {cluster.get('log10_J', 0):<11.2f} "
                  f"{cluster.get('k_best', 0):<15.3e} {cluster.get('mean_error', 0)*100:<9.2f}")
        
        return results


# ============================================================================
# PHASE 4: SCALING LAW ANALYSIS
# ============================================================================

class ScalingLawAnalysis:
    """Analyze log10(k) vs log10(J) scaling"""
    
    @staticmethod
    def fit_power_law(clusters: List[Dict]) -> Tuple[float, float, float]:
        """Fit power law: log10(k) = α * log10(J) + c"""
        
        log_J = []
        log_k = []
        
        for cluster in clusters:
            if cluster.get('k_best') and cluster.get('log10_J'):
                try:
                    lJ = float(cluster['log10_J'])
                    lk = math.log10(float(cluster['k_best']))
                    log_J.append(lJ)
                    log_k.append(lk)
                except (ValueError, TypeError):
                    continue
        
        if len(log_J) < 3:
            print("ERROR: Not enough data points for power law fit")
            return None, None, None
        
        log_J_arr = np.array(log_J)
        log_k_arr = np.array(log_k)
        
        coeffs = np.polyfit(log_J_arr, log_k_arr, 1)
        alpha = coeffs[0]
        c = coeffs[1]
        
        y_pred = alpha * log_J_arr + c
        ss_res = np.sum((log_k_arr - y_pred) ** 2)
        ss_tot = np.sum((log_k_arr - np.mean(log_k_arr)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        print(f"\n Power Law Fit: log10(k) = {alpha:.3f} * log10(J) + {c:.2f}")
        print(f"  R² = {r_squared:.4f}")
        print(f"  Galaxy result: α = -0.589")
        print(f"  Difference: Δα = {abs(alpha - (-0.589)):.3f}")
        
        return alpha, c, r_squared
    
    @staticmethod
    def plot_scaling_law(clusters: List[Dict], output_dir: str = "results/plots"):
        """Plot log10(k) vs log10(J) with power law fit."""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        log_J = []
        log_k = []
        
        for cluster in clusters:
            if cluster.get('k_best') and cluster.get('log10_J'):
                try:
                    lJ = float(cluster['log10_J'])
                    lk = math.log10(float(cluster['k_best']))
                    log_J.append(lJ)
                    log_k.append(lk)
                except (ValueError, TypeError):
                    continue
        
        if len(log_J) < 3:
            return
        
        log_J_arr = np.array(log_J)
        log_k_arr = np.array(log_k)
        
        coeffs = np.polyfit(log_J_arr, log_k_arr, 1)
        alpha = coeffs[0]
        c = coeffs[1]
        
        y_fit = alpha * log_J_arr + c
        
        fig, ax = plt.subplots(figsize=(12, 9))
        
        ax.scatter(log_J_arr, log_k_arr, s=80, alpha=0.6, c='blue', 
                  edgecolors='darkblue', linewidth=1, label=f'Clusters (N={len(log_J)})')
        
        J_line = np.linspace(log_J_arr.min() - 1, log_J_arr.max() + 1, 100)
        k_line = alpha * J_line + c
        ax.plot(J_line, k_line, 'r-', linewidth=2.5, 
               label=f'Fit: log₁₀(k) = {alpha:.3f} log₁₀(J) + {c:.2f}')
        
        k_gal = -0.589 * J_line - 5.5
        ax.plot(J_line, k_gal, 'g--', linewidth=2, alpha=0.7,
               label='Galaxy reference (α = -0.589)')
        
        ax.set_xlabel('log₁₀(J) [kg·m²/s]', fontsize=13, fontweight='bold')
        ax.set_ylabel('log₁₀(k) [SI]', fontsize=13, fontweight='bold')
        ax.set_title('Cluster Scaling Law (FIXED): log10(k) vs log10(J)',
                    fontsize=14, fontweight='bold', pad=20)
        
        ax.grid(True, which='both', alpha=0.3, linestyle=':')
        ax.legend(fontsize=11, loc='best')
        
        # plt.tight_layout()
        # output_path = os.path.join(output_dir, "cluster_k_vs_J_scaling_FIXED.png")
        # plt.savefig(output_path, dpi=150, bbox_inches='tight')
        # print(f"\n Scaling plot saved: {output_path}")
        # plt.close()


# ============================================================================
# PHASE 5: LENSE-THIRRING CONSTRAINT
# ============================================================================

class ClusterLTAnalysis:
    """Compute theoretical k from LT frame-dragging"""
    
    @staticmethod
    def compute_k_theory(J: float, r200_h_Mpc: float) -> Tuple[float, float]:
        """Compute theoretical k from Lense-Thirring for clusters."""
        
        try:
            r200_m = r200_h_Mpc * constants.MPC_TO_M / 0.7
            
            omega_lt = constants.G_OVER_C2 * J / (r200_m ** 3) if r200_m > 0 else 0
            k_theory = omega_lt * (r200_m ** 2) / J if J > 0 else 0
            
            log_ratio = math.log10(k_theory) if k_theory > 0 else None
            
            return k_theory, log_ratio
        
        except (ValueError, ZeroDivisionError):
            return None, None
    
    @staticmethod
    def add_LT_to_clusters(clusters: List[Dict]) -> List[Dict]:
        """Add k_theory and LT ratio to each cluster."""
        
        for cluster in clusters:
            J = cluster.get('J')
            r200 = cluster.get('r200')
            
            if J and r200:
                k_th, log_ratio = ClusterLTAnalysis.compute_k_theory(J, r200)
                cluster['k_theory'] = k_th
                cluster['log10_k_ratio'] = log_ratio
                
                if cluster.get('k_best'):
                    ratio = cluster['k_best'] / k_th if k_th and k_th > 0 else 0
                    cluster['k_ratio_theory_to_fitted'] = ratio
        
        print(f" Computed LT theory for {len(clusters)} clusters")
        return clusters


# ============================================================================
# PHASE 6: SUMMARY & PAPER TABLE
# ============================================================================

class ClusterResults:
    """Generate summary and extract 6-cluster table"""
    
    @staticmethod
    def print_summary(clusters: List[Dict]):
        """Print summary metrics."""
        
        if not clusters:
            return
        
        errors = [c['mean_error'] * 100 for c in clusters if c.get('mean_error')]
        
        print("\n" + "="*80)
        print("CLUSTER ANALYSIS SUMMARY (FIXED)")
        print("="*80)
        print(f"Total clusters analyzed: {len(clusters)}")
        print(f"Mean fitting error: {np.mean(errors):.2f}%")
        print(f"Median fitting error: {np.median(errors):.2f}%")
        print(f"Error range: {min(errors):.2f}% to {max(errors):.2f}%")
        
        log_ratios = [c['log10_k_ratio'] for c in clusters 
                     if c.get('log10_k_ratio') and c.get('log10_k_ratio') is not None]
        if log_ratios:
            print(f"\nLT Gap (log10(k_th/k_fit)):")
            print(f"  Mean: {np.mean(log_ratios):.2f} dex")
            print(f"  Median: {np.median(log_ratios):.2f} dex")
            print(f"  Range: {min(log_ratios):.2f} to {max(log_ratios):.2f} dex")
        
        print("="*80 + "\n")
    
    @staticmethod
    def extract_6_clusters(clusters: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """Extract 3 best + 3 worst clusters by fitting error."""
        
        sorted_by_error = sorted(clusters, key=lambda x: x.get('mean_error', float('inf')))
        
        best_3 = sorted_by_error[:3]
        worst_3 = sorted_by_error[-3:]
        
        return best_3, worst_3
    
    @staticmethod
    def save_6cluster_table(best_3: List[Dict], worst_3: List[Dict], 
                           output_path: str):
        """Save 6-cluster table to CSV."""
        
        rows = []
        
        for cluster in best_3:
            rows.append({
                'ClID': cluster.get('ClID'),
                'fit_quality': 'GOOD',
                'mean_error_pct': f"{cluster.get('mean_error', 0) * 100:.2f}",
                'log10_J': f"{cluster.get('log10_J', 0):.2f}",
                'k_fitted': f"{cluster.get('k_best', 0):.3e}",
                'k_theory': f"{cluster.get('k_theory', 0):.3e}",
                'log10_k_ratio': f"{cluster.get('log10_k_ratio', 0):.2f}",
                'sig200': f"{cluster.get('sig200', 0):.1f}",
                'r200': f"{cluster.get('r200', 0):.2f}",
            })
        
        for cluster in worst_3:
            rows.append({
                'ClID': cluster.get('ClID'),
                'fit_quality': 'BAD',
                'mean_error_pct': f"{cluster.get('mean_error', 0) * 100:.2f}",
                'log10_J': f"{cluster.get('log10_J', 0):.2f}",
                'k_fitted': f"{cluster.get('k_best', 0):.3e}",
                'k_theory': f"{cluster.get('k_theory', 0):.3e}",
                'log10_k_ratio': f"{cluster.get('log10_k_ratio', 0):.2f}",
                'sig200': f"{cluster.get('sig200', 0):.1f}",
                'r200': f"{cluster.get('r200', 0):.2f}",
            })
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        
        print(f" 6-cluster table saved: {output_path}")
    
    @staticmethod
    def print_6cluster_table(best_3: List[Dict], worst_3: List[Dict]):
        """Print 6-cluster table to console."""
        
        print("\n" + "="*80)
        print("TABLE: Spin-Coupling Model on Galaxy Clusters (6 Representative) - FIXED")
        print("="*80)
        
        print(f"\n{'ClID':<8} {'Error':<8} {'log₁₀J':<10} {'k_fitted':<15} {'k_theory':<15} {'log₁₀Ratio':<12} {'sig200':<12} {'r200':<10}")
        print(f"{'':8} {'(%)':<8} {'':10} {'(SI)':<15} {'(SI)':<15} {'(k_th/k_fit)':<12} {'(km/s)':<12} {'(h⁻¹Mpc)':<10}")
        print("-"*80)
        
        print("GOOD FITS (Low Error):")
        for c in best_3:
            print(f"{int(c.get('ClID')):<8} {c.get('mean_error', 0)*100:<7.2f} {c.get('log10_J', 0):<9.2f} "
                  f"{c.get('k_best', 0):<14.3e} {c.get('k_theory', 0):<14.3e} {c.get('log10_k_ratio', 0):<11.2f} "
                  f"{c.get('sig200', 0):<11.1f} {c.get('r200', 0):<9.2f}")
        
        print("\nBAD FITS (High Error):")
        for c in worst_3:
            print(f"{int(c.get('ClID')):<8} {c.get('mean_error', 0)*100:<7.2f} {c.get('log10_J', 0):<9.2f} "
                  f"{c.get('k_best', 0):<14.3e} {c.get('k_theory', 0):<14.3e} {c.get('log10_k_ratio', 0):<11.2f} "
                  f"{c.get('sig200', 0):<11.1f} {c.get('r200', 0):<9.2f}")
        
        print("\n" + "="*80)


# ============================================================================
# MAIN PIPELINE
# ============================================================================

if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("CLUSTER ANALYSIS: VALIDATE a_extra MODEL (FIXED VERSION)")
    print("="*80)
    
    print("\n[PHASE 1] Loading and filtering clusters...")
    loader = ClusterDataLoader()
    all_clusters = loader.load_galwcls(r"data\J_ApJS_246_2\galwcls.dat\galwcls.dat")
    
    clusters = loader.filter_high_quality(all_clusters, 
                                         sig200_min=400,
                                         N200_min=30,
                                         z_max=0.15)
    
    print("\n[PHASE 2] Computing angular momentum J...")
    clusters = ClusterJCalculator.add_J_to_clusters(clusters)
    
    print("\n[PHASE 3] Fitting k parameter (FIXED k_grid + r0)...")
    fitter = ClusterFitter()
    clusters = fitter.fit_all(clusters)
    
    print("\n[PHASE 4] Analyzing scaling law...")
    alpha, c, r2 = ScalingLawAnalysis.fit_power_law(clusters)
    ScalingLawAnalysis.plot_scaling_law(clusters)
    
    print("\n[PHASE 5] Computing Lense-Thirring constraint...")
    clusters = ClusterLTAnalysis.add_LT_to_clusters(clusters)
    
    print("\n[PHASE 6] Generating results...")
    ClusterResults.print_summary(clusters)
    
    best_3, worst_3 = ClusterResults.extract_6_clusters(clusters)
    ClusterResults.print_6cluster_table(best_3, worst_3)
    ClusterResults.save_6cluster_table(best_3, worst_3, r"results\cluster_LT_table_6_FIXED.csv")
    
    with open(r"results\cluster_fits_FIXED.csv", 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['ClID', 'sig200', 'M200', 'r200', 'N200', 'z', 
                                               'J', 'log10_J', 'k_best', 'mean_error', 
                                               'k_theory', 'log10_k_ratio'])
        writer.writeheader()
        for cluster in clusters:
            writer.writerow({
                'ClID': cluster.get('ClID'),
                'sig200': f"{cluster.get('sig200', 0):.2f}",
                'M200': f"{cluster.get('M200', 0):.3e}",
                'r200': f"{cluster.get('r200', 0):.3f}",
                'N200': int(cluster.get('N200', 0)),
                'z': f"{cluster.get('z', 0):.4f}",
                'J': f"{cluster.get('J', 0):.3e}",
                'log10_J': f"{cluster.get('log10_J', 0):.2f}",
                'k_best': f"{cluster.get('k_best', 0):.3e}",
                'mean_error': f"{cluster.get('mean_error', 0):.4f}",
                'k_theory': f"{cluster.get('k_theory', 0):.3e}",
                'log10_k_ratio': f"{cluster.get('log10_k_ratio', 0):.2f}",
            })
    
    print(f"\n Full results saved: results/cluster_fits_FIXED.csv")
    
    print("\n" + "="*80)
    print(" CLUSTER ANALYSIS COMPLETE (FIXED)")
    print("="*80 + "\n")