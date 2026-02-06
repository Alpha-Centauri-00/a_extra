"""
Calculate angular momentum J from visible matter (gas + disk + bulge)
===================================================================

Scans rotation curve .dat files (same format used by the batch fitter),
computes a visible-matter-based J for each galaxy using a half-mass-radius
heuristic, and writes a CSV with the old circular-J and the new J_visible.

Output CSV: results/J_visible_from_data.csv
"""

import os
import math
import csv
from typing import Dict, List, Optional, Tuple
from tools.config import constants




class RotationCurveLoader:
    """Small loader compatible with your existing .dat layout."""

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

            # try to find distance in header
            for line in lines:
                if line.startswith("# Distance"):
                    parts = line.split("=")
                    try:
                        result["distance_mpc"] = float(parts[1].strip().split()[0])
                    except Exception:
                        pass
                    break

            for line in lines:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                # expect at least radius and v_obs
                try:
                    r = float(parts[0])
                    v_obs = float(parts[1])
                    v_gas = float(parts[3]) if len(parts) > 3 else 0.0
                    v_disk = float(parts[4]) if len(parts) > 4 else 0.0
                    v_bul = float(parts[5]) if len(parts) > 5 else 0.0

                    result["radii"].append(r)
                    result["v_obs"].append(v_obs)
                    result["v_gas"].append(v_gas)
                    result["v_disk"].append(v_disk)
                    result["v_bul"].append(v_bul)
                except Exception:
                    continue

        except FileNotFoundError:
            return None

        return result if result["radii"] else None


def calculate_J_circular(data: Dict) -> Optional[float]:
    """Original circular J computed from last data point using v_obs."""
    try:
        r_kpc = data["radii"][-1]
        v_km_s = data["v_obs"][-1]
        r_m = r_kpc * constants.KPC_TO_M
        v_m_s = v_km_s * 1000.0
        M_enclosed = (v_m_s ** 2) * r_m / constants.G
        J = M_enclosed * v_m_s * r_m
        return J
    except Exception:
        return None


def calculate_J_visible(data: Dict) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Compute J using visible matter only.

    Heuristic:
    - compute v_visible(r) = sqrt(v_gas^2 + v_disk^2 + v_bul^2)
    - compute M_enclosed(r) = v_visible(r)^2 * r / G
    - choose r_disk as the radius where M_enclosed reaches 50% of its max
      (half-mass radius heuristic)
    - M_visible = M_enclosed(r_disk)
    - J_visible = M_visible * v_visible(r_disk) * r_disk
    
    Returns:
        (J_visible, r_disk_kpc, v_visible_at_disk)
    """
    try:
        radii = data["radii"]
        v_gas = data["v_gas"]
        v_disk = data["v_disk"]
        v_bul = data["v_bul"]

        if not radii:
            return None, None, None

        M_vals = []
        v_vis_vals = []
        
        for i, r_kpc in enumerate(radii):
            r_m = r_kpc * constants.KPC_TO_M
            
            v_gas_i = v_gas[i] if i < len(v_gas) else 0.0
            v_disk_i = v_disk[i] if i < len(v_disk) else 0.0
            v_bul_i = v_bul[i] if i < len(v_bul) else 0.0
            
            v_vis_km_s = math.sqrt(v_gas_i ** 2 + v_disk_i ** 2 + v_bul_i ** 2)
            v_vis_m_s = v_vis_km_s * 1000.0
            v_vis_vals.append(v_vis_km_s)
            
            # avoid zero-radius or zero-velocity
            if r_m <= 0 or v_vis_m_s <= 0:
                M_vals.append(0.0)
                continue
            
            M_enc = (v_vis_m_s ** 2) * r_m / constants.G
            M_vals.append(M_enc)

        # Check if we have any visible matter
        if not any(v_vis_vals) or max(v_vis_vals) == 0:
            return None, None, None

        if not any(M_vals):
            return None, None, None

        M_max = max(M_vals)
        half_M = 0.5 * M_max

        # Find first radius where M_enc >= half_M
        idx = next((i for i, m in enumerate(M_vals) if m >= half_M), None)
        
        if idx is None:
            # Fallback: use outermost point (maximum extent of data)
            idx = len(M_vals) - 1

        r_disk_kpc = radii[idx]
        r_disk_m = r_disk_kpc * constants.KPC_TO_M
        M_visible = M_vals[idx]
        v_vis_km_s = v_vis_vals[idx]
        v_vis_m_s = v_vis_km_s * 1000.0

        J_visible = M_visible * v_vis_m_s * r_disk_m
        
        return J_visible, r_disk_kpc, v_vis_km_s

    except Exception:
        return None, None, None


def process_directory(data_dir: str, output_csv: str):
    """Process all .dat files in directory and write results."""
    
    if not os.path.exists(data_dir):
        print(f"Data directory not found: {data_dir}")
        return

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    files = sorted([f for f in os.listdir(data_dir) if f.endswith(".dat") or f.endswith("_rotmod.dat")])

    print(f"\n{'='*120}")
    print("CALCULATING J_VISIBLE FROM ROTATION CURVES")
    print(f"{'='*120}\n")
    print(f"Processing {len(files)} galaxies...\n")

    rows = []
    successful = 0
    failed = 0

    print(f"{'Galaxy':<20} {'J_old (circular)':<20} {'J_new (visible)':<20} {'r_disk (kpc)':<15} {'Ratio':<12}")
    print("-" * 120)

    for fname in files:
        path = os.path.join(data_dir, fname)
        gname = fname.replace("_rotmod.dat", "").replace(".dat", "")
        
        data = RotationCurveLoader.load(path)
        if data is None:
            print(f"{gname:<20} {'FAILED':<20} {'(load error)':<20}")
            failed += 1
            continue

        J_old = calculate_J_circular(data)
        J_new, r_disk, v_vis = calculate_J_visible(data)

        if J_new is None:
            print(f"{gname:<20} {J_old:.2e} {'FAILED':<20} {'':<15} {'(calc error)':<12}")
            failed += 1
            continue

        ratio = None
        if J_old and J_new and J_old != 0:
            ratio = J_new / J_old

        rows.append({
            'galaxy_name': gname,
            'J_old_circular': f"{J_old:.6e}" if J_old else "",
            'J_new_visible': f"{J_new:.6e}" if J_new else "",
            'r_disk_kpc': f"{r_disk:.2f}" if r_disk else "",
            'v_visible_at_disk': f"{v_vis:.2f}" if v_vis else "",
            'ratio': f"{ratio:.6f}" if ratio else "",
        })

        print(f"{gname:<20} {J_old:.2e} {J_new:.2e} {r_disk:.2f} {ratio:.6f}")
        successful += 1

    print("\n" + "=" * 120)
    print(f"Results: {successful} successful, {failed} failed out of {len(files)} galaxies")
    print("=" * 120)

    # Write CSV
    with open(output_csv, 'w', newline='') as f:
        fieldnames = ['galaxy_name', 'J_old_circular', 'J_new_visible', 'r_disk_kpc', 'v_visible_at_disk', 'ratio']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"\n✓ Results saved to: {output_csv}\n")


if __name__ == '__main__':
    DATA_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data', 'Rotmod_LTG'))
    OUTPUT_CSV = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'results', 'J_visible_from_data.csv'))

    print(f"Scanning: {DATA_DIR}")
    process_directory(DATA_DIR, OUTPUT_CSV)