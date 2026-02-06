"""
Continuous k Parameter Optimization using J_visible
====================================================

Loads J_visible from J_visible_from_data.csv and rotation curve data.
Uses scipy.optimize.minimize for continuous k fitting (not grid search).
Compares old vs new results and tests k-J scaling law.

Output CSV: results/batch_fitting_results_continuous_k.csv
Comparison: results/optimization_comparison.csv
"""

import os
import math
import csv
from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy.optimize import minimize
from tools.config import constants


class RotationCurveLoader:
    """Load and parse rotation curve .dat files."""
    
    @staticmethod
    def load(filepath: str) -> Optional[Dict]:
        result = {
            "filepath": filepath,
            "distance_mpc": None,
            "radii": [],
            "v_obs": [],
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
                    result["distance_mpc"] = float(parts[1].strip().split()[0])
                    break
            
            for line in lines:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                
                try:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    
                    result["radii"].append(float(parts[0]))
                    result["v_obs"].append(float(parts[1]))
                    result["v_gas"].append(float(parts[3]) if len(parts) > 3 else 0.0)
                    result["v_disk"].append(float(parts[4]) if len(parts) > 4 else 0.0)
                    result["v_bul"].append(float(parts[5]) if len(parts) > 5 else 0.0)
                
                except (ValueError, IndexError):
                    continue
        
        except FileNotFoundError:
            return None
        
        return result if result["radii"] else None


class RotationCurveModel:
    """Spin-mediated acceleration model."""
    
    @staticmethod
    def estimate_v_vis(v_gas_km_s: float, v_disk_km_s: float, v_bul_km_s: float) -> float:
        """Estimate visible matter velocity."""
        return math.sqrt(v_gas_km_s**2 + v_disk_km_s**2 + v_bul_km_s**2)
    
    @staticmethod
    def a_extra(J: float, r_m: float, k: float, r0_m: float = 0.5 * constants.KPC_TO_M, power: float = 2.0) -> float:
        """Extra acceleration from spin coupling."""
        if r_m + r0_m == 0:
            return 0.0
        return k * J / ((r_m + r0_m) ** power)
    
    @staticmethod
    def v_model(v_vis_km_s: float, J: float, r_m: float, k: float, r0_m: float = 0.5 * constants.KPC_TO_M, power: float = 2.0) -> float:
        """Predicted velocity from model."""
        a_extra = RotationCurveModel.a_extra(J, r_m, k, r0_m, power)
        v_vis_m_s = v_vis_km_s * 1000
        v_model_sq = v_vis_m_s ** 2 + r_m * a_extra
        
        if v_model_sq < 0:
            return v_vis_km_s
        
        return (v_model_sq ** 0.5) / 1000


class ContinuousKFitter:
    """Fit k parameter using continuous optimization."""
    
    def __init__(self):
        self.r0_m = 0.5 * constants.KPC_TO_M
        self.power = 2.0
    
        # ---  error_function ---
    def error_function(self, log_k: float, J: float, galaxy_data: Dict) -> float:
        k = 10 ** log_k
        errors = []

        for i, r_kpc in enumerate(galaxy_data["radii"]):
            r_m = r_kpc * constants.KPC_TO_M

            v_vis = RotationCurveModel.estimate_v_vis(
                galaxy_data["v_gas"][i],
                galaxy_data["v_disk"][i],
                galaxy_data["v_bul"][i]
            )
            if v_vis == 0:
                continue

            v_pred = RotationCurveModel.v_model(
                v_vis, J, r_m, k, self.r0_m, self.power
            )
            v_obs = galaxy_data["v_obs"][i]

            if v_obs > 0:
                errors.append(abs(v_pred - v_obs) / v_obs)

        return np.mean(errors) if errors else 1e10

    
            # --- redo fit ---
    def fit(self, galaxy_name: str, J_new: float, galaxy_data: Dict):

        grid_k = [1e-45, 1e-43, 1e-42, 1e-41, 1e-40, 1e-38, 1e-36]
        best_k = min(
            grid_k,
            key=lambda k: self.error_function(math.log10(k), J_new, galaxy_data)
        )

        x0 = np.array([math.log10(best_k)])
        bounds = [(-55, -25)]

        result = minimize(
            self.error_function,
            x0,
            args=(J_new, galaxy_data),
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-10, "maxiter": 5000}
        )

        log_k_opt = result.x[0]
        k_opt = 10 ** log_k_opt

        return k_opt, result.fun



class J_VisibleLoader:
    """Load J_visible values from CSV."""
    
    @staticmethod
    def load_csv(csv_path: str) -> Dict[str, float]:
        """Load J_new values from J_visible CSV.
        
        Returns:
            {galaxy_name: J_new_value}
        """
        j_values = {}
        
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    galaxy_name = row['galaxy_name']
                    try:
                        j_new = float(row['J_new_visible'])
                        j_values[galaxy_name] = j_new
                    except (ValueError, KeyError):
                        continue
        
        except FileNotFoundError:
            print(f"Warning: J_visible CSV not found at {csv_path}")
        
        return j_values


class BatchProcessor:
    """Process multiple galaxies with continuous k optimization."""
    
    def __init__(self, j_visible_csv: str):
        self.fitter = ContinuousKFitter()
        self.j_loader = J_VisibleLoader()
        self.j_visible = self.j_loader.load_csv(j_visible_csv)
        self.results = []
        self.comparison = []
    
    def find_files(self, directory: str) -> List[Tuple[str, str]]:
        """Find all rotation curve files."""
        files = []
        
        if not os.path.exists(directory):
            print(f"Directory not found: {directory}")
            return files
        
        for filename in sorted(os.listdir(directory)):
            if filename.endswith("_rotmod.dat"):
                filepath = os.path.join(directory, filename)
                galaxy_name = filename.replace("_rotmod.dat", "")
                files.append((filepath, galaxy_name))
        
        return files
    
    def process_galaxy(self, filepath: str, galaxy_name: str) -> Optional[Dict]:
        """Process a single galaxy with continuous k optimization."""
        
        print(f"Processing {galaxy_name}...", end=" ")
        
        # Load galaxy data
        data = RotationCurveLoader.load(filepath)
        if data is None:
            print("Failed to load")
            return None
        
        # Get J_new (if available)
        if galaxy_name not in self.j_visible:
            print("J_new not found")
            return None
        
        J_new = self.j_visible[galaxy_name]
        
        # Fit k using continuous optimization
        k_optimal, error_optimal = self.fitter.fit(galaxy_name, J_new, data)
        
        result = {
            "galaxy_name": galaxy_name,
            "distance_mpc": data["distance_mpc"],
            "r_max_kpc": data["radii"][-1],
            "v_max_km_s": data["v_obs"][-1],
            "J_new": J_new,
            "log10_J_new": math.log10(J_new) if J_new > 0 else 0.0,
            "k_optimal": k_optimal,
            "log10_k_optimal": math.log10(k_optimal) if k_optimal > 0 else 0.0,
            "mean_error": error_optimal,
            "mean_error_pct": error_optimal * 100,
        }
        
        self.results.append(result)
        print(f"✓ k={k_optimal:.3e}, error={error_optimal*100:.2f}%")
        
        return result
    
    def process_all(self, directory: str) -> List[Dict]:
        """Process all galaxies in directory."""
        
        files = self.find_files(directory)
        
        if not files:
            print(f"No .dat files found in {directory}")
            return []
        
        print(f"\nFound {len(files)} galaxies to process")
        print(f"Found J_new for {len(self.j_visible)} galaxies\n")
        
        for filepath, galaxy_name in files:
            self.process_galaxy(filepath, galaxy_name)
        
        return self.results
    
    def save_csv(self, output_path: str):
        """Save results to CSV."""
        
        if not self.results:
            print("No results to save")
            return
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                "galaxy_name", "distance_mpc", "r_max_kpc", "v_max_km_s",
                "J_new", "log10_J_new", "k_optimal", "log10_k_optimal", 
                "mean_error", "mean_error_pct"
            ])
            writer.writeheader()
            
            for result in self.results:
                writer.writerow(result)
        
        print(f"\n✓ Results saved to {output_path}")
    
    def print_summary(self):
        """Print summary statistics."""
        
        if not self.results:
            return
        
        errors = [r['mean_error_pct'] for r in self.results]
        k_values = [r['k_optimal'] for r in self.results]
        j_values = [r['J_new'] for r in self.results]
        
        print("\n" + "="*100)
        print("CONTINUOUS k OPTIMIZATION SUMMARY")
        print("="*100)
        print(f"\nGalaxies processed: {len(self.results)}")
        print(f"Mean error: {np.mean(errors):.2f}% (±{np.std(errors):.2f}%)")
        print(f"Median error: {np.median(errors):.2f}%")
        print(f"Best fit: {min(errors):.2f}%")
        print(f"Worst fit: {max(errors):.2f}%")
        print(f"\nk range: {min(k_values):.3e} to {max(k_values):.3e}")
        print(f"J range: {min(j_values):.3e} to {max(j_values):.3e}")
        print("="*100)
    
    def test_scaling_law(self):
        """Test if k ∝ J^-β scaling law holds with new J values."""
        
        if len(self.results) < 5:
            print("Not enough galaxies for scaling law analysis")
            return
        
        log_J = np.array([r['log10_J_new'] for r in self.results])
        log_k = np.array([r['log10_k_optimal'] for r in self.results])
        
        # Check for variation in k values
        k_unique = len(np.unique(log_k))
        if k_unique < 3:
            print("\nWARNING: k values not varying enough for reliable scaling law fit")
            print(f"Only {k_unique} unique k values found across {len(self.results)} galaxies")
            print("This suggests optimizer did not converge properly.")
            return
        
        # Linear regression: log_k = a * log_J + b
        coeffs = np.polyfit(log_J, log_k, 1)
        beta = -coeffs[0]  # negative slope
        a = coeffs[1]
        
        # R-squared
        p = np.poly1d(coeffs)
        y_pred = p(log_J)
        ss_res = np.sum((log_k - y_pred) ** 2)
        ss_tot = np.sum((log_k - np.mean(log_k)) ** 2)
        
        if ss_tot == 0:
            print("\nERROR: Cannot compute R² (zero total sum of squares)")
            return
        
        r_squared = 1 - (ss_res / ss_tot)
        
        print("\n" + "="*100)
        print("SCALING LAW ANALYSIS: k(J)")
        print("="*100)
        print(f"\nFitted form: k = {10**a:.4f} × J^(-{beta:.3f})")
        print(f"R² = {r_squared:.4f}")
        print(f"Exponent β = {beta:.3f} (compare to old 0.589)")
        print("="*100)
        
        return {
            'exponent': beta,
            'prefactor': 10**a,
            'r_squared': r_squared
        }


if __name__ == "__main__":
    
    print("\n" + "="*100)
    print("PHASE B: CONTINUOUS k OPTIMIZATION WITH J_VISIBLE")
    print("="*100)
    
    j_visible_csv = "results/J_visible_from_data.csv"
    data_dir = "data/Rotmod_LTG"
    output_csv = "results/batch_fitting_results_continuous_k.csv"
    
    # Process all galaxies
    processor = BatchProcessor(j_visible_csv)
    processor.process_all(data_dir)
    
    # Save results
    processor.save_csv(output_csv)
    
    # Print summary
    processor.print_summary()
    
    # Test scaling law
    scaling_law = processor.test_scaling_law()
    
    print("\n✓ Continuous k optimization complete!")
    print("="*100)