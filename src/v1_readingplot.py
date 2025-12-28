"""
Phase 6 Plotter: Generate scaling law plots from batch fitting results
=====================================================================

Reads batch_fitting_results.csv and creates publication-quality plots.
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import math

class BatchResultsPlotter:
    """Plot results from batch galaxy fitting."""
    
    def __init__(self, csv_file: str):
        """Initialize with CSV file."""
        self.galaxies = []
        self.J_values = []
        self.k_values = []
        self.errors = []
        self.names = []
        
        self.load_csv(csv_file)
    
    def load_csv(self, filepath: str):
        """Load results from CSV file."""
        try:
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    self.names.append(row['galaxy_name'])
                    self.J_values.append(float(row['J']))
                    self.k_values.append(float(row['k_best']))
                    self.errors.append(float(row['mean_error']))
            
            print(f" Loaded {len(self.names)} galaxies from {filepath}")
        
        except FileNotFoundError:
            print(f"File not found: {filepath}")
    
    def fit_scaling_law(self) -> tuple:
        """Fit k = c * J^alpha in log space."""
        if len(self.J_values) < 2:
            return None, None, None
        
        log_J = np.log10(self.J_values)
        log_k = np.log10(self.k_values)
        
        coeffs = np.polyfit(log_J, log_k, 1)
        alpha = coeffs[0]
        log_c = coeffs[1]
        c = 10 ** log_c
        
        y_pred = alpha * log_J + log_c
        ss_res = np.sum((log_k - y_pred) ** 2)
        ss_tot = np.sum((log_k - np.mean(log_k)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return alpha, c, r2
    
    def plot_k_vs_J(self, save_path: str = None):
        """Main plot: k vs J."""
        if len(self.J_values) < 2:
            print("Need at least 2 data points")
            return
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Data points
        ax.loglog(self.J_values, self.k_values, 'o', markersize=8, 
                 color='steelblue', alpha=0.6, label=f'SPARC galaxies (N={len(self.J_values)})', zorder=3)
        
        # Fit line
        alpha, c, r2 = self.fit_scaling_law()
        
        if alpha is not None:
            J_range = np.logspace(np.log10(min(self.J_values)) - 0.5,
                                  np.log10(max(self.J_values)) + 0.5, 100)
            k_fit = c * (J_range ** alpha)
            ax.loglog(J_range, k_fit, '--', linewidth=2.5, color='red', 
                     label=f'Fit: k ∝ J^{alpha:.2f} (R² = {r2:.3f})', zorder=2)
            
            print(f"\n Scaling Law (N={len(self.J_values)} galaxies):")
            print(f"   k = {c:.3e} × J^{alpha:.3f}")
            print(f"   R² = {r2:.3f}")
        
        ax.set_xlabel(r'Angular Momentum $J$ (kg·m²/s)', fontsize=13, fontweight='bold')
        ax.set_ylabel(r'Coupling Constant $k$', fontsize=13, fontweight='bold')
        ax.set_title(f'Scaling Law: k vs J ({len(self.J_values)} SPARC Galaxies)\nDirect from rotation curve files', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, which='both', alpha=0.3, linestyle=':')
        ax.legend(fontsize=12, loc='best', framealpha=0.9)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f" Saved: {save_path}")
        
        plt.show()
    
    def plot_error_distribution(self, save_path: str = None):
        """Plot error distribution."""
        if not self.errors:
            return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        error_percent = [e * 100 for e in self.errors]
        
        ax.hist(error_percent, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
        ax.axvline(np.mean(error_percent), color='red', linestyle='--', linewidth=2, 
                  label=f'Mean = {np.mean(error_percent):.2f}%')
        ax.axvline(np.median(error_percent), color='green', linestyle='--', linewidth=2,
                  label=f'Median = {np.median(error_percent):.2f}%')
        
        ax.set_xlabel('Fitting Error (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Galaxies', fontsize=12, fontweight='bold')
        ax.set_title(f'Error Distribution ({len(self.errors)} galaxies)', fontsize=13, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f" Saved: {save_path}")
        
        plt.show()
    
    def plot_error_vs_J(self, save_path: str = None):
        """Plot error vs angular momentum."""
        if len(self.J_values) < 2:
            return
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        error_percent = [e * 100 for e in self.errors]
        
        ax.loglog(self.J_values, error_percent, 'o', markersize=7, 
                 color='darkgreen', alpha=0.6, label='Fitting error', zorder=3)
        
        mean_error = np.mean(error_percent)
        ax.axhline(mean_error, color='red', linestyle='--', linewidth=2,
                  label=f'Mean error = {mean_error:.2f}%')
        
        ax.set_xlabel(r'Angular Momentum $J$ (kg·m²/s)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Fitting Error (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Error vs Galaxy Size ({len(self.J_values)} galaxies)', fontsize=13, fontweight='bold')
        ax.grid(True, which='both', alpha=0.3, linestyle=':')
        ax.legend(fontsize=11)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f" Saved: {save_path}")
        
        plt.show()
    
    def print_summary(self):
        """Print summary statistics."""
        if not self.errors:
            return
        
        error_percent = [e * 100 for e in self.errors]
        
        print("\n" + "="*80)
        print("BATCH FITTING SUMMARY STATISTICS")
        print("="*80)
        print(f"Total galaxies: {len(self.names)}")
        print(f"J range: {min(self.J_values):.2e} to {max(self.J_values):.2e} kg·m²/s")
        print(f"k range: {min(self.k_values):.2e} to {max(self.k_values):.2e}")
        print(f"\nError statistics:")
        print(f"  Mean:   {np.mean(error_percent):.2f}%")
        print(f"  Median: {np.median(error_percent):.2f}%")
        print(f"  Std:    {np.std(error_percent):.2f}%")
        print(f"  Min:    {np.min(error_percent):.2f}%")
        print(f"  Max:    {np.max(error_percent):.2f}%")
        
        alpha, c, r2 = self.fit_scaling_law()
        if alpha is not None:
            print(f"\nScaling law:")
            print(f"  k = {c:.3e} × J^{alpha:.3f}")
            print(f"  R² = {r2:.3f}")
        
        print("="*80 + "\n")

if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("PHASE 6: SCALING LAW PLOTS FROM BATCH FITTING RESULTS")
    print("="*80 + "\n")
    
    plotter = BatchResultsPlotter("results/batch_fitting_results.csv")
    plotter.print_summary()
    
    print("Generating plots...\n")
    
    plotter.plot_k_vs_J(save_path="results/phase6_k_vs_J_full_sample.png")
    plotter.plot_error_distribution(save_path="results/phase6_error_distribution.png")
    plotter.plot_error_vs_J(save_path="results/phase6_error_vs_J_full_sample.png")
    
    print("\n" + "="*80)
    print(" All plots generated!")
    print("="*80)