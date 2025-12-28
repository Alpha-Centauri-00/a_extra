"""
Batch Galaxy Fitter: Automated k parameter fitting for multiple galaxies
=========================================================================

Loads rotation curve files, calculates J, fits k, and saves results.
All data read directly from data > Rotmod_LTG > .dat files.
"""

import os
import math
import csv
from typing import List, Tuple, Dict
import numpy as np

# Physical constants
M_sun = 1.989e30
kpc_to_m = 3.086e19
G = 6.674e-11


class RotationCurveLoader:
    """Load and parse rotation curve .dat files."""
    
    @staticmethod
    def load(filepath: str) -> Dict:
        """Load rotation curve file."""
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
    
    @staticmethod
    def calculate_J(data: Dict) -> Tuple[float, float, float]:
        """Calculate J from last data point."""
        if not data or len(data["radii"]) == 0:
            return None, None, None
        
        r_kpc = data["radii"][-1]
        v_km_s = data["v_obs"][-1]
        
        r_m = r_kpc * kpc_to_m
        v_m_s = v_km_s * 1000
        
        M_enclosed = (v_m_s ** 2) * r_m / G
        J = M_enclosed * v_m_s * r_m
        
        return J, math.log10(J), v_km_s


class RotationCurveModel:
    """Spin-mediated acceleration model."""
    
    @staticmethod
    def a_extra(J: float, r_m: float, k: float, r0_m: float = 0.5 * kpc_to_m, power: float = 2.0) -> float:
        """Extra acceleration from spin coupling."""
        if r_m + r0_m == 0:
            return 0.0
        return k * J / ((r_m + r0_m) ** power)
    
    @staticmethod
    def v_model(v_vis_km_s: float, J: float, r_m: float, k: float, r0_m: float = 0.5 * kpc_to_m, power: float = 2.0) -> float:
        """Predicted velocity from model."""
        a_extra = RotationCurveModel.a_extra(J, r_m, k, r0_m, power)
        v_vis_m_s = v_vis_km_s * 1000
        v_model_sq = v_vis_m_s ** 2 + r_m * a_extra
        
        if v_model_sq < 0:
            return v_vis_km_s
        
        return (v_model_sq ** 0.5) / 1000
    
    @staticmethod
    def estimate_v_vis(v_gas_km_s: float, v_disk_km_s: float, v_bul_km_s: float) -> float:
        """Estimate visible matter velocity."""
        return math.sqrt(v_gas_km_s**2 + v_disk_km_s**2 + v_bul_km_s**2)


class GalaxyFitter:
    """Fit k parameter for a galaxy."""
    
    def __init__(self):
        self.k_range = [1e-42, 1e-40, 1e-38, 1e-36, 1e-34, 1e-32, 1e-30]
        self.r0_m = 0.5 * kpc_to_m
        self.power = 2.0
    
    def fit(self, galaxy_name: str, data: Dict, J: float) -> Tuple[float, float]:
        """Fit k parameter. Returns (k_best, mean_error)."""
        
        best_k = None
        best_error = float('inf')
        
        for k in self.k_range:
            predictions = []
            errors = []
            
            for i, r_kpc in enumerate(data["radii"]):
                r_m = r_kpc * kpc_to_m
                
                v_vis = RotationCurveModel.estimate_v_vis(
                    data["v_gas"][i],
                    data["v_disk"][i],
                    data["v_bul"][i]
                )
                
                if v_vis == 0:
                    continue
                
                v_pred = RotationCurveModel.v_model(v_vis, J, r_m, k, self.r0_m, self.power)
                predictions.append(v_pred)
                
                v_obs = data["v_obs"][i]
                if v_obs > 0:
                    error = abs(v_pred - v_obs) / v_obs
                    errors.append(error)
            
            if errors:
                mean_error = sum(errors) / len(errors)
                
                if mean_error < best_error:
                    best_error = mean_error
                    best_k = k
        
        return best_k, best_error


class BatchProcessor:
    """Process multiple galaxies."""
    
    def __init__(self, directory: str):
        self.directory = directory
        self.results = []
        self.fitter = GalaxyFitter()
    
    def find_files(self) -> List[Tuple[str, str]]:
        """Find all rotation curve files."""
        files = []
        
        if not os.path.exists(self.directory):
            print(f"Directory not found: {self.directory}")
            return files
        
        for filename in sorted(os.listdir(self.directory)):
            if filename.endswith("_rotmod.dat"):
                filepath = os.path.join(self.directory, filename)
                galaxy_name = filename.replace("_rotmod.dat", "")
                files.append((filepath, galaxy_name))
        
        return files
    
    def process_galaxy(self, filepath: str, galaxy_name: str) -> Dict:
        """Process a single galaxy."""
        
        print(f"Processing {galaxy_name}...", end=" ")
        
        data = RotationCurveLoader.load(filepath)
        if data is None:
            print("Failed to load")
            return None
        
        J, log10_J, v_max = RotationCurveLoader.calculate_J(data)
        if J is None:
            print("Failed to calculate J")
            return None
        
        k_best, error = self.fitter.fit(galaxy_name, data, J)
        if k_best is None:
            print("Failed to fit k")
            return None
        
        result = {
            "galaxy_name": galaxy_name,
            "filepath": filepath,
            "distance_mpc": data["distance_mpc"],
            "r_max_kpc": data["radii"][-1],
            "v_max_km_s": v_max,
            "J": J,
            "log10_J": log10_J,
            "k_best": k_best,
            "mean_error": error,
        }
        
        self.results.append(result)
        print(f"✓ J={J:.2e}, k={k_best:.2e}, error={error*100:.2f}%")
        
        return result
    
    def process_all(self) -> List[Dict]:
        """Process all galaxies in directory."""
        
        files = self.find_files()
        
        if not files:
            print(f"No .dat files found in {self.directory}")
            return []
        
        print(f"\nFound {len(files)} galaxies to process\n")
        
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
                "J", "log10_J", "k_best", "mean_error"
            ])
            writer.writeheader()
            
            for result in self.results:
                writer.writerow({
                    "galaxy_name": result["galaxy_name"],
                    "distance_mpc": result["distance_mpc"],
                    "r_max_kpc": result["r_max_kpc"],
                    "v_max_km_s": result["v_max_km_s"],
                    "J": f"{result['J']:.3e}",
                    "log10_J": f"{result['log10_J']:.2f}",
                    "k_best": f"{result['k_best']:.3e}",
                    "mean_error": f"{result['mean_error']:.4f}",
                })
        
        print(f"\n✓ Results saved to {output_path}")
    
    def print_summary(self):
        """Print summary table."""
        
        if not self.results:
            return
        
        print("\n" + "="*110)
        print("BATCH FITTING SUMMARY")
        print("="*110)
        print()
        print(f"{'Galaxy':<20} {'Distance (Mpc)':<15} {'r_max (kpc)':<12} {'v_max (km/s)':<14} {'J':<16} {'k_best':<14} {'Error (%)':<10}")
        print("-"*110)
        
        for result in sorted(self.results, key=lambda x: x["J"]):
            print(f"{result['galaxy_name']:<20} {result['distance_mpc']:<15.2f} {result['r_max_kpc']:<12.2f} {result['v_max_km_s']:<14.1f} {result['J']:<16.3e} {result['k_best']:<14.3e} {result['mean_error']*100:<10.2f}")
        
        print()

if __name__ == "__main__":
    
    print("\n" + "="*110)
    print("BATCH GALAXY FITTER: AUTOMATED k PARAMETER FITTING")
    print("="*110)
    
    processor = BatchProcessor("data/Rotmod_LTG")
    processor.process_all()
    processor.print_summary()
    processor.save_csv("results/batch_fitting_results.csv")
    
    print("\n" + "="*110)
    print("✓ Batch processing complete!")
    print("="*110)
