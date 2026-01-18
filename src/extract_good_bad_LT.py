"""
6-Galaxy Lense-Thirring Analysis for Paper Section 4.1
======================================================

Produces:
(a) Detailed comparison table (markdown/LaTeX) for 6 representative galaxies
(b) Plot A: log-log k_theory vs k_fitted with 6 points highlighted
(c) Plot B: log10(J) vs log10(k_ratio) showing under-prediction factor

Selected galaxies:
  GOOD FITS (3): F563-V1, UGC11914, NGC3953
  BAD FITS (3):  UGC11557, NGC4389, UGC02455
"""

import csv
import math
import os
import matplotlib.pyplot as plt
import numpy as np


def load_full_analysis(csv_path: str) -> dict:
    """Load full LT analysis CSV into dictionary keyed by galaxy name."""
    data = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['galaxy_name']
            data[name] = {
                'galaxy_name': name,
                'distance_mpc': float(row['distance_mpc']),
                'r_max_kpc': float(row['r_max_kpc']),
                'v_max_km_s': float(row['v_max_km_s']),
                'J': float(row['J']),
                'log10_J': float(row['log10_J']),
                'k_fitted': float(row['k_fitted']),
                'log10_k_fitted': float(row['log10_k_fitted']),
                'mean_error_pct': float(row['mean_error_pct']),
                'omega_lt_max_rad_s': float(row['omega_lt_max_rad_s']),
                'v_boost_max_km_s': float(row['v_boost_max_km_s']),
                'relative_boost_pct': float(row['relative_boost_pct']),
                'k_theory': float(row['k_theory']),
                'log10_k_theory': float(row['log10_k_theory']),
                'k_ratio_theory_to_fitted': float(row['k_ratio_theory_to_fitted']),
                'log10_k_ratio': float(row['log10_k_ratio']),
            }
    return data


def print_detailed_table(selected_galaxies: list, all_data: dict):
    """Print detailed comparison table for 6 galaxies."""
    
    print("\n" + "="*180)
    print("TABLE: Lense-Thirring Frame-Dragging vs. Empirical Spin Coupling")
    print("6 Representative Galaxies (3 Good Fits + 3 Bad Fits)")
    print("="*180)
    
    print("\nDETAILED METRICS TABLE")
    print("-"*180)
    print(f"{'Galaxy':<15} {'Error':<8} {'log₁₀J':<10} {'k_fitted':<15} {'k_theory':<15} {'log₁₀Ratio':<12} {'ΩLT(rmax)':<16} {'δv(rmax)':<12}")
    print(f"{'':15} {'(%)':<8} {'':10} {'(SI)':<15} {'(SI)':<15} {'(k_th/k_fit)':<12} {'(rad/s)':<16} {'(km/s)':<12}")
    print("-"*180)
    
    # Print data in order: good fits first, then bad fits
    good_fits = selected_galaxies[:3]
    bad_fits = selected_galaxies[3:]
    
    print("GOOD FITS (Low Error):")
    for name in good_fits:
        d = all_data[name]
        print(f"{d['galaxy_name']:<15} {d['mean_error_pct']:<7.2f} {d['log10_J']:<9.2f} "
              f"{d['k_fitted']:<14.3e} {d['k_theory']:<14.3e} {d['log10_k_ratio']:<11.2f} "
              f"{d['omega_lt_max_rad_s']:<15.3e} {d['v_boost_max_km_s']:<11.4f}")
    
    print("\nBAD FITS (High Error):")
    for name in bad_fits:
        d = all_data[name]
        print(f"{d['galaxy_name']:<15} {d['mean_error_pct']:<7.2f} {d['log10_J']:<9.2f} "
              f"{d['k_fitted']:<14.3e} {d['k_theory']:<14.3e} {d['log10_k_ratio']:<11.2f} "
              f"{d['omega_lt_max_rad_s']:<15.3e} {d['v_boost_max_km_s']:<11.4f}")
    
    print("\n" + "="*180)
    print("INTERPRETATION:")
    print("  k_theory: GR prediction from Lense-Thirring frame-dragging")
    print("  k_fitted: Empirical coupling constant from rotation curve fitting")
    print("  log₁₀(Ratio): Theory under-predicts by 4-6 orders of magnitude")
    print("  Pattern holds for both good and bad fits → universal scaling")
    print("="*180)


def save_table_csv(selected_galaxies: list, all_data: dict, output_path: str):
    """Save 6-galaxy table as CSV."""
    
    rows = []
    for name in selected_galaxies:
        d = all_data[name]
        rows.append({
            'galaxy_name': d['galaxy_name'],
            'fit_quality': 'GOOD' if selected_galaxies.index(name) < 3 else 'BAD',
            'mean_error_pct': f"{d['mean_error_pct']:.2f}",
            'log10_J': f"{d['log10_J']:.2f}",
            'k_fitted_SI': f"{d['k_fitted']:.3e}",
            'k_theory_SI': f"{d['k_theory']:.3e}",
            'log10_k_ratio': f"{d['log10_k_ratio']:.2f}",
            'omega_lt_max_rad_s': f"{d['omega_lt_max_rad_s']:.3e}",
            'v_boost_max_km_s': f"{d['v_boost_max_km_s']:.4f}",
            'relative_boost_pct': f"{d['relative_boost_pct']:.2f}",
        })
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"\n✓ Table saved to {output_path}")


def plot_a_k_comparison_6gal(selected_galaxies: list, all_data: dict, output_dir: str = "results/plots"):
    """
    Plot A: Log-log scatter of k_theory vs k_fitted
    - All 175 galaxies as gray background
    - 6 selected galaxies highlighted with different markers/colors
    - Good fits: blue circles
    - Bad fits: red squares
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # All galaxies as gray background
    all_k_theory = [all_data[g]['k_theory'] for g in all_data.keys()]
    all_k_fitted = [all_data[g]['k_fitted'] for g in all_data.keys()]
    
    ax.scatter(all_k_fitted, all_k_theory, s=30, alpha=0.15, c='gray', 
              label='Other galaxies (N=169)', zorder=1)
    
    # Good fits (blue circles)
    good_k_fitted = [all_data[g]['k_fitted'] for g in selected_galaxies[:3]]
    good_k_theory = [all_data[g]['k_theory'] for g in selected_galaxies[:3]]
    good_names = selected_galaxies[:3]
    
    ax.scatter(good_k_fitted, good_k_theory, s=300, marker='o', c='blue', 
              alpha=0.7, edgecolors='darkblue', linewidth=2, label='Good fits (error < 10%)', zorder=3)
    
    # Bad fits (red squares)
    bad_k_fitted = [all_data[g]['k_fitted'] for g in selected_galaxies[3:]]
    bad_k_theory = [all_data[g]['k_theory'] for g in selected_galaxies[3:]]
    bad_names = selected_galaxies[3:]
    
    ax.scatter(bad_k_fitted, bad_k_theory, s=300, marker='s', c='red', 
              alpha=0.7, edgecolors='darkred', linewidth=2, label='Bad fits (error > 80%)', zorder=3)
    
    # 1:1 reference line
    k_min = 1e-49
    k_max = 1e-41
    ax.plot([k_min, k_max], [k_min, k_max], 'r--', linewidth=2.5, alpha=0.8, 
           label='1:1 (perfect match)', zorder=2)
    
    # Annotations for good fits
    for i, name in enumerate(good_names):
        ax.annotate(name, (good_k_fitted[i], good_k_theory[i]), 
                   xytext=(10, 10), textcoords='offset points',
                   fontsize=10, fontweight='bold', color='darkblue',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='darkblue', lw=1.5))
    
    # Annotations for bad fits
    for i, name in enumerate(bad_names):
        ax.annotate(name, (bad_k_fitted[i], bad_k_theory[i]), 
                   xytext=(10, -15), textcoords='offset points',
                   fontsize=10, fontweight='bold', color='darkred',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='darkred', lw=1.5))
    
    # Log scales and labels
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('k_fitted (SI units)', fontsize=13, fontweight='bold')
    ax.set_ylabel('k_theory (SI units)', fontsize=13, fontweight='bold')
    ax.set_title('Plot A: Lense-Thirring Theory vs. Empirical Fit\n6 Representative Galaxies in Context',
                fontsize=14, fontweight='bold', pad=20)
    
    ax.grid(True, which='both', alpha=0.3, linestyle=':')
    ax.legend(fontsize=11, loc='upper left', framealpha=0.95)
    
    # Text box
    textstr = ("Red dashed line: Perfect match (theory = fitted)\n"
              "Actual: k_fitted >> k_theory by 4-6 orders\n"
              "Explanation: Disk integration + coherence factors")
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='bottom', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "plot_A_k_theory_vs_fitted_6gal.png")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Plot A saved: {output_path}")
    plt.close()


def plot_b_j_vs_ratio(selected_galaxies: list, all_data: dict, output_dir: str = "results/plots"):
    """
    Plot B: log10(J) vs log10(k_ratio)
    - All 175 galaxies as gray background
    - 6 selected galaxies highlighted
    - Shows under-prediction factor is independent of J
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # All galaxies background
    all_J = [all_data[g]['log10_J'] for g in all_data.keys()]
    all_ratio = [all_data[g]['log10_k_ratio'] for g in all_data.keys()]
    
    ax.scatter(all_J, all_ratio, s=30, alpha=0.15, c='gray', 
              label='Other galaxies (N=169)', zorder=1)
    
    # Good fits
    good_J = [all_data[g]['log10_J'] for g in selected_galaxies[:3]]
    good_ratio = [all_data[g]['log10_k_ratio'] for g in selected_galaxies[:3]]
    good_names = selected_galaxies[:3]
    
    ax.scatter(good_J, good_ratio, s=300, marker='o', c='blue', 
              alpha=0.7, edgecolors='darkblue', linewidth=2, label='Good fits', zorder=3)
    
    # Bad fits
    bad_J = [all_data[g]['log10_J'] for g in selected_galaxies[3:]]
    bad_ratio = [all_data[g]['log10_k_ratio'] for g in selected_galaxies[3:]]
    bad_names = selected_galaxies[3:]
    
    ax.scatter(bad_J, bad_ratio, s=300, marker='s', c='red', 
              alpha=0.7, edgecolors='darkred', linewidth=2, label='Bad fits', zorder=3)
    
    # Horizontal band showing typical ratio range
    ax.axhspan(-6, -4, alpha=0.1, color='green', label='Typical range: -6 to -4 dex')
    
    # Annotations for good fits
    for i, name in enumerate(good_names):
        ax.annotate(name, (good_J[i], good_ratio[i]), 
                   xytext=(10, 10), textcoords='offset points',
                   fontsize=10, fontweight='bold', color='darkblue',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='darkblue', lw=1.5))
    
    # Annotations for bad fits
    for i, name in enumerate(bad_names):
        ax.annotate(name, (bad_J[i], bad_ratio[i]), 
                   xytext=(10, -15), textcoords='offset points',
                   fontsize=10, fontweight='bold', color='darkred',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', color='darkred', lw=1.5))
    
    ax.set_xlabel('log₁₀(J) [kg·m²/s]', fontsize=13, fontweight='bold')
    ax.set_ylabel('log₁₀(k_theory / k_fitted)', fontsize=13, fontweight='bold')
    ax.set_title('Plot B: Under-Prediction Factor vs. Angular Momentum\nUniversal Scaling Across Galaxy Types',
                fontsize=14, fontweight='bold', pad=20)
    
    ax.grid(True, which='both', alpha=0.3, linestyle=':')
    ax.legend(fontsize=11, loc='best', framealpha=0.95)
    
    # Text box
    textstr = ("Key: Ratio is ~constant across all J\n"
              "No special J where theory recovers\n"
              "Consistent 4-6 dex gap (disk integration)")
    ax.text(0.02, 0.02, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='bottom', horizontalalignment='left',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "plot_B_J_vs_k_ratio_6gal.png")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Plot B saved: {output_path}")
    plt.close()


if __name__ == "__main__":
    
    # Define 6 selected galaxies
    selected_galaxies = ['F563-V1', 'UGC11914', 'NGC3953', 'UGC11557', 'NGC4389', 'UGC02455']
    
    # Paths
    analysis_csv = r"results\lense_thirring_analysis.csv"  # Use raw string for Windows paths
    output_dir = "results/plots"
    table_csv = "results/lt_paper_table_6galaxies.csv"
    
    print("\n" + "="*180)
    print("6-GALAXY LENSE-THIRRING ANALYSIS FOR PAPER SECTION 4.1")
    print("="*180)
    
    # Load data
    print("\nLoading analysis data...")
    all_data = load_full_analysis(analysis_csv)
    print(f"✓ Loaded {len(all_data)} galaxies")
    
    # Print detailed table
    print_detailed_table(selected_galaxies, all_data)
    
    # Save table as CSV
    save_table_csv(selected_galaxies, all_data, table_csv)
    
    # Generate plots
    print("\nGenerating plots...")
    plot_a_k_comparison_6gal(selected_galaxies, all_data, output_dir)
    plot_b_j_vs_ratio(selected_galaxies, all_data, output_dir)
    
    print("\n" + "="*180)
    print("✓ COMPLETE: All outputs ready for paper Section 4.1")
    print("="*180)
    print(f"\nOutputs:")
    print(f"  - Detailed table (console + {table_csv})")
    print(f"  - Plot A: {output_dir}/plot_A_k_theory_vs_fitted_6gal.png")
    print(f"  - Plot B: {output_dir}/plot_B_J_vs_k_ratio_6gal.png")
    print("="*180 + "\n")