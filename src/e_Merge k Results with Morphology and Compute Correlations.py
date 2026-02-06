"""
Merge fitting results with galaxy morphology and compute correlations
=====================================================================

Merges:
- batch_fitting_results_continuous_k.csv (k, error)
- galaxy_properties_h_bt.csv (h, B/T, SB_disk_center)

Computes Pearson correlations between k and morphological properties.
Generates diagnostic scatter plots.
"""

import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
from scipy.stats import pearsonr


class DataMerger:
    """Merge k results with morphology data."""
    
    @staticmethod
    def load_k_results(csv_path: str) -> Dict[str, Dict]:
        """Load k optimization results."""
        data = {}
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    gname = row['galaxy_name']
                    data[gname] = {
                        'k_optimal': float(row['k_optimal']),
                        'log10_k_optimal': float(row['log10_k_optimal']),
                        'mean_error': float(row['mean_error']),
                        'mean_error_pct': float(row['mean_error_pct']),
                        'J_new': float(row['J_new']),
                        'log10_J_new': float(row['log10_J_new']),
                    }
        except FileNotFoundError:
            print(f"File not found: {csv_path}")
        return data
    
    @staticmethod
    def load_morphology(csv_path: str) -> Dict[str, Dict]:
        """Load h and B/T data."""
        data = {}
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    gname = row['galaxy_name']
                    h = None
                    bt = None
                    sb_disk = None
                    
                    try:
                        if row['h_kpc']:
                            h = float(row['h_kpc'])
                    except (ValueError, KeyError):
                        pass
                    
                    try:
                        if row['B_T']:
                            bt = float(row['B_T'])
                    except (ValueError, KeyError):
                        pass
                    
                    try:
                        if row['SB_disk_center']:
                            sb_disk = float(row['SB_disk_center'])
                    except (ValueError, KeyError):
                        pass
                    
                    data[gname] = {
                        'h_kpc': h,
                        'B_T': bt,
                        'SB_disk_center': sb_disk,
                    }
        except FileNotFoundError:
            print(f"File not found: {csv_path}")
        return data
    
    @staticmethod
    def merge(k_results: Dict, morphology: Dict) -> List[Dict]:
        """Merge two datasets on galaxy_name."""
        merged = []
        
        for gname in sorted(k_results.keys()):
            if gname not in morphology:
                continue
            
            row = {
                'galaxy_name': gname,
                **k_results[gname],
                **morphology[gname],
            }
            merged.append(row)
        
        return merged
    
    @staticmethod
    def save_merged(data: List[Dict], output_path: str):
        """Save merged data to CSV."""
        if not data:
            print("No data to save")
            return
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        fieldnames = [
            'galaxy_name', 'k_optimal', 'log10_k_optimal', 'mean_error_pct',
            'J_new', 'log10_J_new', 'h_kpc', 'B_T', 'SB_disk_center'
        ]
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in data:
                writer.writerow({k: row.get(k, '') for k in fieldnames})
        
        print(f"✓ Merged data saved to: {output_path}")


class CorrelationAnalyzer:
    """Compute correlations and generate plots."""
    
    @staticmethod
    def compute_correlations(data: List[Dict]) -> Dict:
        """Compute Pearson correlations between k and morphology."""
        
        correlations = {}
        
        # Extract data
        log_k = np.array([r['log10_k_optimal'] for r in data])
        h = np.array([r['h_kpc'] if r['h_kpc'] is not None else np.nan for r in data])
        bt = np.array([r['B_T'] if r['B_T'] is not None else np.nan for r in data])
        sb = np.array([r['SB_disk_center'] if r['SB_disk_center'] is not None else np.nan for r in data])
        error = np.array([r['mean_error_pct'] for r in data])
        j_vis = np.array([r['log10_J_new'] for r in data])
        
        # Compute correlations (ignore NaN)
        pairs = [
            ('log_k vs log_J', log_k, j_vis),
            ('log_k vs h', log_k, h),
            ('log_k vs B/T', log_k, bt),
            ('log_k vs SB_disk', log_k, sb),
            ('error vs B/T', error, bt),
            ('error vs h', error, h),
        ]
        
        print("\n" + "="*100)
        print("CORRELATION ANALYSIS")
        print("="*100 + "\n")
        
        for name, x, y in pairs:
            # Filter NaN values
            valid = ~(np.isnan(x) | np.isnan(y))
            x_clean = x[valid]
            y_clean = y[valid]
            
            if len(x_clean) < 3:
                print(f"{name:<25} | N/A (insufficient data)")
                continue
            
            r, p_value = pearsonr(x_clean, y_clean)
            correlations[name] = {'r': r, 'p': p_value, 'n': len(x_clean)}
            
            sig = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            print(f"{name:<25} | r = {r:+.4f}  (p = {p_value:.4f}) {sig:3s}  [N={len(x_clean)}]")
        
        print("\n" + "="*100 + "\n")
        
        return correlations


class DiagnosticPlotter:
    """Generate diagnostic scatter plots."""
    
    @staticmethod
    def plot_k_vs_h_colored_by_bt(data: List[Dict], output_dir: str = "results/plots"):
        """Plot k vs h, colored by B/T."""
        os.makedirs(output_dir, exist_ok=True)
        
        log_k = [r['log10_k_optimal'] for r in data if r['h_kpc'] is not None]
        h = [r['h_kpc'] for r in data if r['h_kpc'] is not None]
        bt = [r['B_T'] if r['B_T'] is not None else 0 for r in data if r['h_kpc'] is not None]
        
        fig, ax = plt.subplots(figsize=(11, 8))
        
        scatter = ax.scatter(h, log_k, c=bt, s=100, alpha=0.6, cmap='RdYlBu_r', edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Disk Scale Length h (kpc)', fontsize=12, fontweight='bold')
        ax.set_ylabel('log₁₀(k_optimal)', fontsize=12, fontweight='bold')
        ax.set_title('Coupling Constant vs Disk Scale Length', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        cbar = plt.colorbar(scatter, ax=ax, label='Bulge-to-Total (B/T)')
        
        plt.tight_layout()
        output_path = os.path.join(output_dir, "k_vs_h_colored_by_bt.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {output_path}")
        plt.close()
    
    @staticmethod
    def plot_k_vs_bt_colored_by_error(data: List[Dict], output_dir: str = "results/plots"):
        """Plot k vs B/T, colored by fitting error."""
        os.makedirs(output_dir, exist_ok=True)
        
        log_k = [r['log10_k_optimal'] for r in data if r['B_T'] is not None]
        bt = [r['B_T'] for r in data if r['B_T'] is not None]
        error = [r['mean_error_pct'] for r in data if r['B_T'] is not None]
        
        fig, ax = plt.subplots(figsize=(11, 8))
        
        scatter = ax.scatter(bt, log_k, c=error, s=100, alpha=0.6, cmap='RdYlGn_r', edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Bulge-to-Total (B/T)', fontsize=12, fontweight='bold')
        ax.set_ylabel('log₁₀(k_optimal)', fontsize=12, fontweight='bold')
        ax.set_title('Coupling Constant vs Bulge Prominence', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        cbar = plt.colorbar(scatter, ax=ax, label='Fitting Error (%)')
        
        plt.tight_layout()
        output_path = os.path.join(output_dir, "k_vs_bt_colored_by_error.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {output_path}")
        plt.close()
    
    @staticmethod
    def plot_k_vs_sb_colored_by_error(data: List[Dict], output_dir: str = "results/plots"):
        """Plot k vs SB_disk_center, colored by fitting error."""
        os.makedirs(output_dir, exist_ok=True)
        
        log_k = [r['log10_k_optimal'] for r in data if r['SB_disk_center'] is not None]
        sb = [r['SB_disk_center'] for r in data if r['SB_disk_center'] is not None]
        error = [r['mean_error_pct'] for r in data if r['SB_disk_center'] is not None]
        
        fig, ax = plt.subplots(figsize=(11, 8))
        
        scatter = ax.scatter(sb, log_k, c=error, s=100, alpha=0.6, cmap='RdYlGn_r', edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Central Surface Brightness (L☉/pc²)', fontsize=12, fontweight='bold')
        ax.set_ylabel('log₁₀(k_optimal)', fontsize=12, fontweight='bold')
        ax.set_title('Coupling Constant vs Disk Surface Brightness', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        cbar = plt.colorbar(scatter, ax=ax, label='Fitting Error (%)')
        
        plt.tight_layout()
        output_path = os.path.join(output_dir, "k_vs_sb_colored_by_error.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {output_path}")
        plt.close()
    
    @staticmethod
    def plot_error_vs_bt(data: List[Dict], output_dir: str = "results/plots"):
        """Plot error vs B/T to test if bulges hurt the fit."""
        os.makedirs(output_dir, exist_ok=True)
        
        error = [r['mean_error_pct'] for r in data if r['B_T'] is not None]
        bt = [r['B_T'] for r in data if r['B_T'] is not None]
        
        fig, ax = plt.subplots(figsize=(11, 8))
        
        scatter = ax.scatter(bt, error, s=100, alpha=0.6, c='steelblue', edgecolors='black', linewidth=0.5)
        
        # Add trend line
        z = np.polyfit(bt, error, 1)
        p = np.poly1d(z)
        bt_line = np.linspace(min(bt), max(bt), 100)
        ax.plot(bt_line, p(bt_line), "r--", linewidth=2, label=f'Trend: y={z[0]:.2f}x+{z[1]:.2f}')
        
        ax.set_xlabel('Bulge-to-Total (B/T)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Fitting Error (%)', fontsize=12, fontweight='bold')
        ax.set_title('Model Performance vs Bulge Prominence', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10)
        
        plt.tight_layout()
        output_path = os.path.join(output_dir, "error_vs_bt.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {output_path}")
        plt.close()


if __name__ == "__main__":
    
    k_results_csv = "results/batch_fitting_results_continuous_k.csv"
    morphology_csv = "results/galaxy_properties_h_bt.csv"
    merged_csv = "results/merged_k_morphology.csv"
    plot_dir = "results/plots"
    
    print("\n" + "="*100)
    print("MERGING k RESULTS WITH MORPHOLOGY DATA")
    print("="*100 + "\n")
    
    # Load data
    k_results = DataMerger.load_k_results(k_results_csv)
    morphology = DataMerger.load_morphology(morphology_csv)
    
    print(f"Loaded {len(k_results)} k results")
    print(f"Loaded {len(morphology)} morphology entries")
    
    # Merge
    merged = DataMerger.merge(k_results, morphology)
    print(f"Merged {len(merged)} common galaxies\n")
    
    # Save merged
    DataMerger.save_merged(merged, merged_csv)
    
    # Compute correlations
    print("\n" + "="*100)
    print("PHASE C: CORRELATION ANALYSIS")
    print("="*100)
    
    correlations = CorrelationAnalyzer.compute_correlations(merged)
    
    # Generate plots
    print("\n" + "="*100)
    print("GENERATING DIAGNOSTIC PLOTS")
    print("="*100 + "\n")
    
    DiagnosticPlotter.plot_k_vs_h_colored_by_bt(merged, plot_dir)
    DiagnosticPlotter.plot_k_vs_bt_colored_by_error(merged, plot_dir)
    DiagnosticPlotter.plot_k_vs_sb_colored_by_error(merged, plot_dir)
    DiagnosticPlotter.plot_error_vs_bt(merged, plot_dir)
    
    print("\n" + "="*100)
    print("✓ PHASE C COMPLETE: Correlation analysis and diagnostic plots ready")
    print("="*100 + "\n")