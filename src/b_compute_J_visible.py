"""
Calculate angular momentum J from visible matter (gas + disk + bulge)
===================================================================

UPDATED (Task 3):
  - Propagates J_vis uncertainty from SPARC inclination (i) and distance (D)
    error bars, which are the two dominant sources of scatter.

  Uncertainty propagation:
    J_vis  ∝  M_vis * v_vis * r_disk
    M_vis  ∝  v_vis^2 * r_disk / G
    v_vis  ∝  v_vis_observed  (corrected by 1/sin(i))
    r_disk ∝  angular_scale * D

  So:   J_vis ∝ v_vis^3 * r_disk^2 / G   (roughly)

  Partial derivatives (log-space):
    d(lnJ) / d(lni)  = -3  (v_vis ∝ 1/sin(i), dominant at large i)
    d(lnJ) / d(lnD)  = +2  (r_disk ∝ D)

  We propagate these analytically and report sigma_J as fractional uncertainty.

Output CSV: results/J_visible_from_data.csv  (adds J_err, J_lower, J_upper columns)
"""

import os
import math
import csv
import numpy as np
from typing import Dict, List, Optional, Tuple
from tools.config import constants


# =============================================================================
# SPARC INCLINATION / DISTANCE ERROR LOOKUP
# (from SPARC Table 1: Lelli et al. 2016, ApJS 227, 21)
# These are fiducial fractional uncertainties used when per-galaxy values
# are not available in the data files.
# =============================================================================

# Typical SPARC error budgets (fractional, 1-sigma):
DEFAULT_FRAC_ERR_DISTANCE    = 0.10   # ~10% distance uncertainty
DEFAULT_FRAC_ERR_INCLINATION = 0.05   # ~5% inclination uncertainty (in sin(i))

# For edge-on galaxies (i > 80 deg), sin(i) ~ 1 and inclination error is tiny
# For face-on galaxies (i < 30 deg), inclination error dominates strongly
# We read i from the header when available; fall back to 60 deg if not found.
DEFAULT_INCLINATION_DEG = 60.0


# =============================================================================
# ROTATION CURVE LOADER
# =============================================================================

class RotationCurveLoader:
    """Load and parse rotation curve .dat files."""

    @staticmethod
    def load(filepath: str) -> Optional[Dict]:
        result = {
            "filepath":       filepath,
            "distance_mpc":   None,
            "distance_err":   None,   # NEW: distance error in Mpc (if in header)
            "inclination_deg": None,  # NEW: inclination in degrees (if in header)
            "incl_err_deg":   None,   # NEW: inclination error in degrees
            "radii":          [],
            "v_obs":          [],
            "v_err":          [],
            "v_gas":          [],
            "v_disk":         [],
            "v_bul":          [],
        }

        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if not line.startswith("#"):
                    continue
                low = line.lower()

                # Distance: "# Distance = 7.72 Mpc"
                if "distance" in low and "=" in line:
                    try:
                        result["distance_mpc"] = float(line.split("=")[1].strip().split()[0])
                    except Exception:
                        pass

                # Distance error: "# e_Distance = 0.77 Mpc"  (some files have this)
                if "e_distance" in low and "=" in line:
                    try:
                        result["distance_err"] = float(line.split("=")[1].strip().split()[0])
                    except Exception:
                        pass

                # Inclination: "# Inc = 65 deg"
                if "inc" in low and "=" in line and "e_inc" not in low:
                    try:
                        val = float(line.split("=")[1].strip().split()[0])
                        result["inclination_deg"] = val
                    except Exception:
                        pass

                # Inclination error: "# e_Inc = 3 deg"
                if "e_inc" in low and "=" in line:
                    try:
                        result["incl_err_deg"] = float(line.split("=")[1].strip().split()[0])
                    except Exception:
                        pass

            # Parse data rows
            for line in lines:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                try:
                    r      = float(parts[0])
                    v_obs  = float(parts[1])
                    v_err  = float(parts[2]) if len(parts) > 2 else 0.0
                    v_gas  = float(parts[3]) if len(parts) > 3 else 0.0
                    v_disk = float(parts[4]) if len(parts) > 4 else 0.0
                    v_bul  = float(parts[5]) if len(parts) > 5 else 0.0

                    result["radii"].append(r)
                    result["v_obs"].append(v_obs)
                    result["v_err"].append(v_err)
                    result["v_gas"].append(v_gas)
                    result["v_disk"].append(v_disk)
                    result["v_bul"].append(v_bul)
                except Exception:
                    continue

        except FileNotFoundError:
            return None

        return result if result["radii"] else None


# =============================================================================
# J_VIS CALCULATOR
# =============================================================================

def calculate_J_circular(data: Dict) -> Optional[float]:
    """Original circular J from last point (v_obs)."""
    try:
        r_m   = data["radii"][-1] * constants.KPC_TO_M
        v_m_s = data["v_obs"][-1] * 1000.0
        M     = v_m_s**2 * r_m / constants.G
        return M * v_m_s * r_m
    except Exception:
        return None


def calculate_J_visible(
    data: Dict,
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    J_vis from baryonic components only.

    Returns (J_visible, r_disk_kpc, v_visible_at_disk)
    """
    try:
        radii  = data["radii"]
        v_gas  = data["v_gas"]
        v_disk = data["v_disk"]
        v_bul  = data["v_bul"]

        if not radii:
            return None, None, None

        M_vals    = []
        v_vis_vals = []

        for i, r_kpc in enumerate(radii):
            r_m = r_kpc * constants.KPC_TO_M
            vg  = v_gas[i]  if i < len(v_gas)  else 0.0
            vd  = v_disk[i] if i < len(v_disk) else 0.0
            vb  = v_bul[i]  if i < len(v_bul)  else 0.0

            v_vis_kms = math.sqrt(vg**2 + vd**2 + vb**2)
            v_vis_ms  = v_vis_kms * 1000.0
            v_vis_vals.append(v_vis_kms)

            if r_m <= 0 or v_vis_ms <= 0:
                M_vals.append(0.0)
                continue
            M_vals.append(v_vis_ms**2 * r_m / constants.G)

        if not any(v_vis_vals) or max(v_vis_vals) == 0:
            return None, None, None

        M_max  = max(M_vals)
        half_M = 0.5 * M_max
        idx    = next((i for i, m in enumerate(M_vals) if m >= half_M),
                      len(M_vals) - 1)

        r_disk_kpc = radii[idx]
        r_disk_m   = r_disk_kpc * constants.KPC_TO_M
        M_visible  = M_vals[idx]
        v_vis_kms  = v_vis_vals[idx]
        v_vis_ms   = v_vis_kms * 1000.0

        J_visible = M_visible * v_vis_ms * r_disk_m
        return J_visible, r_disk_kpc, v_vis_kms

    except Exception:
        return None, None, None


# =============================================================================
# UNCERTAINTY PROPAGATION  (Task 3)
# =============================================================================

def propagate_J_uncertainty(
    J_visible: float,
    data: Dict,
    r_disk_kpc: float,
) -> Tuple[float, float, float, float, float]:
    """
    Propagate distance and inclination uncertainties into sigma_J.

    Physical model:
    ---------------
    The velocity components in the .dat files are already inclination-corrected
    (SPARC divides by sin(i)).  The radius r_disk is physical (kpc) and
    scales linearly with distance D.

    Scaling relations (log-space):
        J_vis ≈ M_vis * v_vis * r_disk
              ≈ [v_vis^2 * r_disk / G] * v_vis * r_disk
              = v_vis^3 * r_disk^2 / G

    With v_vis ∝ 1/sin(i)  and  r_disk ∝ D:

        ln J  ≈ 3*ln(v_vis) + 2*ln(r_disk) - ln(G)
              = 3*ln(1/sin(i)) + 2*ln(D) + const

    Partial derivatives:
        ∂(ln J)/∂(ln sin(i))  = -3
        ∂(ln J)/∂(ln D)       = +2

    So:
        (σ_J / J)^2 = 9*(σ_sin_i / sin(i))^2 + 4*(σ_D / D)^2

    where σ_sin_i / sin(i) ≈ σ_i * |cos(i)| / sin(i)  = σ_i * |cot(i)|

    Returns:
        (sigma_J_frac, sigma_J_abs, J_lower, J_upper, sigma_J_from_D_frac)
    """
    # Resolve inclination
    i_deg  = data.get("inclination_deg") or DEFAULT_INCLINATION_DEG
    di_deg = data.get("incl_err_deg")    or (DEFAULT_FRAC_ERR_INCLINATION
                                              * i_deg)

    i_rad  = math.radians(i_deg)
    di_rad = math.radians(di_deg)

    sin_i  = math.sin(i_rad)
    cos_i  = abs(math.cos(i_rad))

    # Guard against face-on singularity
    if sin_i < 0.05:
        sin_i = 0.05

    frac_err_sin_i = di_rad * cos_i / sin_i   # σ(sin i) / sin i

    # Resolve distance fractional error
    D    = data.get("distance_mpc")   or 1.0
    dD   = data.get("distance_err")   or (DEFAULT_FRAC_ERR_DISTANCE * D)
    frac_err_D = dD / D if D > 0 else DEFAULT_FRAC_ERR_DISTANCE

    # Quadrature propagation
    sigma_J_frac = math.sqrt(
        9.0 * frac_err_sin_i**2 +
        4.0 * frac_err_D**2
    )

    # Contribution breakdown
    sigma_from_incl_frac = 3.0 * frac_err_sin_i
    sigma_from_D_frac    = 2.0 * frac_err_D

    sigma_J_abs = sigma_J_frac * J_visible
    J_lower     = J_visible * (1.0 - sigma_J_frac)
    J_upper     = J_visible * (1.0 + sigma_J_frac)

    return (sigma_J_frac, sigma_J_abs,
            J_lower, J_upper,
            sigma_from_incl_frac, sigma_from_D_frac)


# =============================================================================
# BATCH PROCESSOR
# =============================================================================

def process_directory(data_dir: str, output_csv: str):

    if not os.path.exists(data_dir):
        print(f"Data directory not found: {data_dir}")
        return

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    files = sorted(
        f for f in os.listdir(data_dir)
        if f.endswith(".dat") or f.endswith("_rotmod.dat")
    )

    header = (f"{'Galaxy':<20} {'J_old (circ)':<18} {'J_vis':<18}"
              f" {'r_disk (kpc)':<13} {'sigma_J (%)':<13}"
              f" {'from_incl (%)':<14} {'from_D (%)':<12}"
              f" {'Ratio':<10}")
    divider = "-" * len(header)

    print(f"\n{'='*len(header)}")
    print("CALCULATING J_VISIBLE WITH UNCERTAINTY PROPAGATION")
    print(f"{'='*len(header)}\n")
    print(f"Processing {len(files)} galaxies …\n")
    print(header)
    print(divider)

    rows       = []
    successful = 0
    failed     = 0

    for fname in files:
        path  = os.path.join(data_dir, fname)
        gname = fname.replace("_rotmod.dat", "").replace(".dat", "")

        data = RotationCurveLoader.load(path)
        if data is None:
            print(f"  {gname:<20} LOAD FAILED")
            failed += 1
            continue

        J_old = calculate_J_circular(data)
        J_new, r_disk, v_vis = calculate_J_visible(data)

        if J_new is None:
            print(f"  {gname:<20} CALC FAILED")
            failed += 1
            continue

        # Uncertainty propagation
        (sigma_frac, sigma_abs,
         J_lower, J_upper,
         sig_incl_frac, sig_D_frac) = propagate_J_uncertainty(J_new, data, r_disk)

        ratio = (J_new / J_old) if J_old else None

        rows.append({
            'galaxy_name':          gname,
            'J_old_circular':       f"{J_old:.6e}" if J_old else "",
            'J_new_visible':        f"{J_new:.6e}",
            'r_disk_kpc':           f"{r_disk:.2f}",
            'v_visible_at_disk':    f"{v_vis:.2f}",
            'ratio':                f"{ratio:.6f}" if ratio else "",
            # NEW uncertainty columns
            'sigma_J_frac':         f"{sigma_frac:.4f}",
            'sigma_J_pct':          f"{sigma_frac*100:.2f}",
            'sigma_J_abs':          f"{sigma_abs:.4e}",
            'J_lower_1sigma':       f"{J_lower:.6e}",
            'J_upper_1sigma':       f"{J_upper:.6e}",
            'sigma_from_incl_pct':  f"{sig_incl_frac*100:.2f}",
            'sigma_from_D_pct':     f"{sig_D_frac*100:.2f}",
            'inclination_deg':      f"{data.get('inclination_deg') or DEFAULT_INCLINATION_DEG:.1f}",
            'distance_mpc':         f"{data.get('distance_mpc') or '?'}",
        })

        print(f"  {gname:<20} {J_old:.2e}  {J_new:.2e}"
              f"  {r_disk:>9.2f}     "
              f"  {sigma_frac*100:>9.1f}%"
              f"  {sig_incl_frac*100:>10.1f}%"
              f"  {sig_D_frac*100:>8.1f}%"
              f"  {ratio:.4f}" if ratio else "")
        successful += 1

    print(divider)
    print(f"\nResults: {successful} successful, {failed} failed "
          f"out of {len(files)} galaxies")

    # ----------------------------------------------------------------
    # Summary statistics on the uncertainty
    # ----------------------------------------------------------------
    if rows:
        sig_pcts = [float(r['sigma_J_pct']) for r in rows]
        print(f"\nJ_vis uncertainty summary (1-sigma):")
        print(f"  Mean  σ_J : {np.mean(sig_pcts):.1f}%")
        print(f"  Median σ_J: {np.median(sig_pcts):.1f}%")
        print(f"  Range     : {min(sig_pcts):.1f}% – {max(sig_pcts):.1f}%")
        print(f"\n  The scatter in log10(k) vs log10(J_vis) should be compared")
        print(f"  to these uncertainties to assess how much is intrinsic vs noise.")

    # Write CSV
    fieldnames = [
        'galaxy_name', 'J_old_circular', 'J_new_visible',
        'r_disk_kpc', 'v_visible_at_disk', 'ratio',
        'sigma_J_frac', 'sigma_J_pct', 'sigma_J_abs',
        'J_lower_1sigma', 'J_upper_1sigma',
        'sigma_from_incl_pct', 'sigma_from_D_pct',
        'inclination_deg', 'distance_mpc',
    ]
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n✓ Results saved to: {output_csv}\n")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    DATA_DIR   = os.path.normpath(
        os.path.join(os.path.dirname(__file__), '..', 'data', 'Rotmod_LTG'))
    OUTPUT_CSV = os.path.normpath(
        os.path.join(os.path.dirname(__file__), '..', 'results',
                     'J_visible_from_data.csv'))

    print(f"Scanning: {DATA_DIR}")
    process_directory(DATA_DIR, OUTPUT_CSV)