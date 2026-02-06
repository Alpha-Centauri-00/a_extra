"""
Extract disk scale length (h) and bulge-to-total (B/T) from .dens files
========================================================================

Loads .dens files from data/BulgeDiskDec_LTG/
For each galaxy:
  - Fit exponential to SBdisk(r) to extract h
  - Integrate SBdisk and SBbulge profiles to compute B/T
  - Output: CSV with (galaxy_name, h_kpc, B_T, SB_disk_center, SB_bulge_center)
"""

import os
import csv
import math
import numpy as np
from typing import Dict, Optional, Tuple
from scipy.optimize import curve_fit
# from scipy.integrate import trapz
from scipy.integrate import trapezoid



class DensFileLoader:
    """Load and parse .dens files."""
    
    @staticmethod
    def load(filepath: str) -> Optional[Dict]:
        """Load .dens file with 3 columns: Rad[kpc], SBdisk[Lsun/pc^2], SBbulge[Lsun/pc^2]
        
        Handles multiple spaces and tabs as delimiters.
        """
        result = {
            "radii": [],
            "sb_disk": [],
            "sb_bulge": [],
        }
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                
                # Skip header and empty lines
                if not line or line.startswith("#"):
                    continue
                
                # Split on whitespace (handles multiple spaces/tabs)
                parts = line.split()
                
                try:
                    if len(parts) < 3:
                        continue
                    
                    rad = float(parts[0])
                    sb_disk = float(parts[1])
                    sb_bulge = float(parts[2])
                    
                    result["radii"].append(rad)
                    result["sb_disk"].append(sb_disk)
                    result["sb_bulge"].append(sb_bulge)
                
                except (ValueError, IndexError):
                    continue
        
        except FileNotFoundError:
            return None
        
        return result if result["radii"] else None


class PhotometricDecomposition:
    """Extract h and B/T from surface brightness profiles."""
    
    @staticmethod
    def fit_disk_scale_length(radii: list, sb_disk: list) -> Optional[float]:
        """
        Fit exponential profile to SBdisk(r): SB = SB0 * exp(-r/h)
        
        Fits only the outer part where disk dominates (to avoid bulge contribution).
        Uses points where SB_disk > 0.
        
        Returns:
            h in kpc
        """
        
        if not radii or not sb_disk or len(radii) < 3:
            return None
        
        # Filter: only use points with SB_disk > 0
        valid_indices = [i for i, sb in enumerate(sb_disk) if sb > 0]
        
        if len(valid_indices) < 3:
            return None
        
        r_fit = np.array([radii[i] for i in valid_indices])
        sb_fit = np.array([sb_disk[i] for i in valid_indices])
        
        # Use outer half of data (where exponential is cleaner)
        n = len(r_fit)
        r_fit = r_fit[n//2:]
        sb_fit = sb_fit[n//2:]
        
        if len(r_fit) < 3:
            # Fallback: use all valid points
            r_fit = np.array([radii[i] for i in valid_indices])
            sb_fit = np.array([sb_disk[i] for i in valid_indices])
        
        try:
            # Fit: log(SB) = log(SB0) - r/h
            # Linear fit to log(SB) vs r
            log_sb = np.log(sb_fit)
            coeffs = np.polyfit(r_fit, log_sb, 1)
            
            # slope = -1/h
            slope = coeffs[0]
            
            if slope >= 0:  # Should be negative for exponential decay
                return None
            
            h = -1.0 / slope
            
            # Sanity check: h should be reasonable (0.1 - 100 kpc)
            if h < 0.1 or h > 100:
                return None
            
            return h
        
        except Exception:
            return None
    
    @staticmethod
    def compute_bulge_to_total(radii: list, sb_disk: list, sb_bulge: list) -> Optional[float]:
        """
        Compute B/T = L_bulge / (L_bulge + L_disk)
        
        Integrates surface brightness profiles assuming axisymmetry:
        L ∝ ∫ SB(r) * 2πr dr
        
        Returns:
            B/T ratio (0 to 1)
        """
        
        if not radii or not sb_disk or not sb_bulge or len(radii) < 2:
            return None
        
        try:
            radii_arr = np.array(radii)
            sb_disk_arr = np.array(sb_disk)
            sb_bulge_arr = np.array(sb_bulge)
            
            # Luminosity element: L ∝ SB(r) * 2πr * dr
            # We integrate SB * r (the 2π cancels in the ratio)
            integrand_disk = sb_disk_arr * radii_arr
            integrand_bulge = sb_bulge_arr * radii_arr
            
            L_disk = trapezoid(integrand_disk, radii_arr)
            L_bulge = trapezoid(integrand_bulge, radii_arr)
            
            if L_disk + L_bulge == 0:
                return None
            
            B_T = L_bulge / (L_disk + L_bulge)
            
            # Sanity check
            if B_T < 0 or B_T > 1:
                return None
            
            return B_T
        
        except Exception:
            return None


class BatchDensProcessor:
    """Process all .dens files and extract parameters."""
    
    def __init__(self):
        self.results = []
    
    def process_directory(self, dens_dir: str) -> list:
        """Process all .dens files in directory."""
        
        if not os.path.exists(dens_dir):
            print(f"Directory not found: {dens_dir}")
            return []
        
        files = sorted([f for f in os.listdir(dens_dir) if f.endswith(".dens")])
        
        print("\n" + "="*120)
        print("EXTRACTING DISK SCALE LENGTH (h) AND BULGE-TO-TOTAL (B/T)")
        print("="*120)
        print(f"\nProcessing {len(files)} .dens files...\n")
        
        print(f"{'Galaxy':<20} {'h (kpc)':<12} {'B/T':<10} {'SB_disk(0)':<14} {'SB_bulge(0)':<14} {'Status':<15}")
        print("-"*120)
        
        successful = 0
        failed = 0
        
        for fname in files:
            galaxy_name = fname.replace(".dens", "")
            filepath = os.path.join(dens_dir, fname)
            
            data = DensFileLoader.load(filepath)
            
            if data is None:
                print(f"{galaxy_name:<20} {'':<12} {'':<10} {'':<14} {'':<14} {'LOAD FAILED':<15}")
                failed += 1
                continue
            
            # Extract h
            h = PhotometricDecomposition.fit_disk_scale_length(data["radii"], data["sb_disk"])
            
            # Extract B/T
            bt = PhotometricDecomposition.compute_bulge_to_total(data["radii"], data["sb_disk"], data["sb_bulge"])
            
            # Get central values
            sb_disk_center = data["sb_disk"][0] if data["sb_disk"] else None
            sb_bulge_center = data["sb_bulge"][0] if data["sb_bulge"] else None
            
            if h is None or bt is None:
                status = "CALC FAILED"
                failed += 1
            else:
                status = "✓"
                successful += 1
            
            result = {
                'galaxy_name': galaxy_name,
                'h_kpc': f"{h:.4f}" if h else "",
                'B_T': f"{bt:.4f}" if bt else "",
                'SB_disk_center': f"{sb_disk_center:.2f}" if sb_disk_center else "",
                'SB_bulge_center': f"{sb_bulge_center:.2f}" if sb_bulge_center else "",
                'h_numeric': h,
                'bt_numeric': bt,
            }
            
            self.results.append(result)
            
            print(f"{galaxy_name:<20} {(f'{h:.4f}' if h else ''):<12} {(f'{bt:.4f}' if bt else ''):<10} {(f'{sb_disk_center:.2f}' if sb_disk_center else ''):<14} {(f'{sb_bulge_center:.2f}' if sb_bulge_center else ''):<14} {status:<15}")
        
        print("\n" + "="*120)
        print(f"Results: {successful} successful, {failed} failed out of {len(files)} galaxies")
        print("="*120)
        
        return self.results
    
    def save_csv(self, output_path: str):
        """Save results to CSV."""
        
        if not self.results:
            print("No results to save")
            return
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', newline='') as f:
            fieldnames = ['galaxy_name', 'h_kpc', 'B_T', 'SB_disk_center', 'SB_bulge_center']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in self.results:
                row = {k: result[k] for k in fieldnames}
                writer.writerow(row)
        
        print(f"\n✓ Results saved to: {output_path}\n")


if __name__ == "__main__":
    
    dens_dir = "data/BulgeDiskDec_LTG"
    output_csv = "results/galaxy_properties_h_bt.csv"
    
    processor = BatchDensProcessor()
    processor.process_directory(dens_dir)
    processor.save_csv(output_csv)