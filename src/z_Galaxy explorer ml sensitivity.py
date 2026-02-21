import pandas as pd
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import sys
from scipy.optimize import minimize
import math
from tools.config import constants


# --- Configuration ---
DATA_FOLDER    = r"data\Rotmod_LTG"
RESULTS_CSV    = r"results\batch_fitting_results_continuous_k.csv"
GLOBAL_SUMMARY = r"results\global_fit_results_summary.txt"

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Inter", "Segoe UI", "Roboto", "Arial"],
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})


# =============================================================================
# Helper: parse global_fit_results_summary.txt
# =============================================================================

def load_global_params(path):
    """Read k0/beta from the summary file written by c_fit_k_continuous.py."""
    out = {
        "model_A_k0": None, "model_A_chi2": None,
        "model_B_k0": None, "model_B_beta": None, "model_B_chi2": None,
        "preferred": None,  "delta_AIC": None,
    }
    if not os.path.exists(path):
        return out
    try:
        block = None
        for raw in open(path, encoding="utf-8", errors="ignore"):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Model A"):
                block = "A"; continue
            if line.startswith("Model B"):
                block = "B"; continue
            if "MODEL COMPARISON" in line:
                block = "CMP"; continue
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip().lower()
            val = val.strip()
            def fv():
                try:    return float(val)
                except: return None
            if block == "A":
                if key == "k0":             out["model_A_k0"]   = fv()
                elif key == "chi2_reduced": out["model_A_chi2"] = fv()
            elif block == "B":
                if key == "k0":             out["model_B_k0"]   = fv()
                elif key == "beta":         out["model_B_beta"] = fv()
                elif key == "chi2_reduced": out["model_B_chi2"] = fv()
            elif block == "CMP":
                if "preferred" in key:      out["preferred"]    = val
                elif "delta aic" in key:    out["delta_AIC"]    = fv()
    except Exception as e:
        print(f"[load_global_params] {e}")
    return out


# =============================================================================
# Main application
# =============================================================================

class GalaxyExplorer:
    def __init__(self, root, results_csv, data_folder, global_summary):

        self.galaxy_buttons = {}
        self.active_button  = None
        self.current_galaxy = None
        self.current_curve  = None

        self.root = root
        try:
            root.iconbitmap("data/Untitled.ico")
        except:
            pass

        self.root.title("SPARC Spin-Coupling Explorer - M/L Sensitivity Test")
        self.root.geometry("1500x900")
        self.data_folder = data_folder

        try:
            self.df_params   = pd.read_csv(results_csv)
            self.galaxy_list = sorted(self.df_params["galaxy_name"].unique().tolist())
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.root.destroy()
            return

        # Load global fit results
        self.gp = load_global_params(global_summary)

        # Toggle variables: both start OFF (= original per-galaxy behaviour)
        self.use_global_model = tk.BooleanVar(value=False)
        self.use_model_B      = tk.BooleanVar(value=False)

        # ── Sidebar ──────────────────────────────────────────────────
        self.sidebar = ctk.CTkFrame(root, width=280, fg_color="#242424")
        self.sidebar.pack(side="left", fill="y", padx=10, pady=10)

        ctk.CTkLabel(
            self.sidebar, text="Galaxies",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(pady=(10, 5))

        # Global params + toggle panel
        self._build_global_panel()

        self.list_frame = ctk.CTkScrollableFrame(self.sidebar)
        self.list_frame.pack(fill="both", expand=True, padx=5, pady=5)

        for name in self.galaxy_list:
            btn = ctk.CTkButton(
                self.list_frame, text=name, anchor="w",
                fg_color="transparent", hover_color="#3A3A3A",
                command=lambda n=name: self.on_galaxy_click(n)
            )
            btn.pack(fill="x", padx=5, pady=2)
            self.galaxy_buttons[name] = btn

        # ── Main content  (identical to original) ────────────────────
        self.main_frame = ctk.CTkFrame(root)
        self.main_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        self.plot_frame = ctk.CTkFrame(self.main_frame)
        self.plot_frame.pack(fill="both", expand=True)

        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.toolbar.pack_forget()

        btn_frame = ctk.CTkFrame(self.plot_frame, fg_color="#242424")
        btn_frame.pack(side="bottom", fill="x")

        controls = ctk.CTkFrame(btn_frame, fg_color="transparent")
        controls.pack(expand=True, pady=6)

        ctk.CTkButton(controls, text="Home", command=self.toolbar.home).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Pan",  command=self.toolbar.pan).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Zoom", command=self.toolbar.zoom).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Save", command=self.toolbar.save_figure).pack(side="left", padx=5)

        # M/L Slider (identical to original)
        self.ml_control_frame = ctk.CTkFrame(self.main_frame, fg_color="#1a1a1a")
        self.ml_control_frame.pack(fill="x", pady=(10, 0))

        ctk.CTkLabel(
            self.ml_control_frame,
            text="Stellar Mass-to-Light Ratio (M/L) Sensitivity Test",
            font=ctk.CTkFont(size=14, weight="bold")
        ).pack(pady=(10, 5))

        ctk.CTkLabel(
            self.ml_control_frame,
            text="Lower M/L = Less stellar mass -> Blue line drops -> Spin-coupling needs larger k",
            font=ctk.CTkFont(size=11), text_color="#AAAAAA"
        ).pack()

        slider_frame = ctk.CTkFrame(self.ml_control_frame, fg_color="transparent")
        slider_frame.pack(fill="x", padx=20, pady=10)

        ctk.CTkLabel(
            slider_frame, text="M/L Scale:",
            font=ctk.CTkFont(size=12, weight="bold")
        ).pack(side="left", padx=(0, 10))

        self.ml_slider = ctk.CTkSlider(
            slider_frame, from_=0.2, to=1.0,
            number_of_steps=40, command=self.on_ml_change
        )
        self.ml_slider.set(1.0)
        self.ml_slider.pack(side="left", fill="x", expand=True, padx=5)

        self.ml_value_label = ctk.CTkLabel(
            slider_frame, text="1.00",
            font=ctk.CTkFont(size=12, weight="bold"), width=50
        )
        self.ml_value_label.pack(side="left", padx=(10, 0))

        info_frame = ctk.CTkFrame(self.ml_control_frame, fg_color="transparent")
        info_frame.pack(fill="x", padx=20, pady=(0, 10))

        self.k_comparison_label = ctk.CTkLabel(
            info_frame, text="Select a galaxy to see M/L sensitivity",
            font=ctk.CTkFont(size=11), text_color="#CCCCCC"
        )
        self.k_comparison_label.pack()

        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.ax.text(0.5, 0.5, "Select a galaxy from the list",
                     ha="center", va="center", fontsize=14,
                     transform=self.ax.transAxes)
        self.apply_ctk_theme()
        self.canvas.draw()

    # =========================================================================
    # Global params panel + toggles  (new addition in sidebar)
    # =========================================================================

    def _build_global_panel(self):
        gp = self.gp

        outer = ctk.CTkFrame(self.sidebar, fg_color="#1e1e1e", corner_radius=8)
        outer.pack(fill="x", padx=6, pady=(0, 8))

        ctk.CTkLabel(
            outer, text="Global Fit  (Tasks 1 & 2)",
            font=ctk.CTkFont(size=12, weight="bold"), text_color="#AADDFF"
        ).pack(pady=(8, 4))

        # Model A
        rA = ctk.CTkFrame(outer, fg_color="#2a2a2a", corner_radius=6)
        rA.pack(fill="x", padx=8, pady=3)
        ctk.CTkLabel(rA, text="Model A  (universal k)",
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color="#88CCFF").pack(anchor="w", padx=8, pady=(4, 1))
        if gp["model_A_k0"] is not None:
            tA = f"k0 = {gp['model_A_k0']:.3e}   chi2 = {gp['model_A_chi2']:.2f}"
        else:
            tA = "not available  (re-run script c)"
        ctk.CTkLabel(rA, text=tA, font=ctk.CTkFont(size=9),
                     text_color="#CCCCCC").pack(anchor="w", padx=14, pady=(0, 5))

        # Model B
        rB = ctk.CTkFrame(outer, fg_color="#2a2a2a", corner_radius=6)
        rB.pack(fill="x", padx=8, pady=3)
        ctk.CTkLabel(rB, text="Model B  (power-law k*J^b)",
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color="#FFCC88").pack(anchor="w", padx=8, pady=(4, 1))
        if gp["model_B_k0"] is not None:
            tB = (f"k0 = {gp['model_B_k0']:.3e}\n"
                  f"b  = {gp['model_B_beta']:.4f}   chi2 = {gp['model_B_chi2']:.2f}")
        else:
            tB = "not available  (re-run script c)"
        ctk.CTkLabel(rB, text=tB, font=ctk.CTkFont(size=9),
                     text_color="#CCCCCC").pack(anchor="w", padx=14, pady=(0, 5))

        # Preferred
        if gp["preferred"]:
            pt = f"Preferred: {gp['preferred']}"
            if gp["delta_AIC"] is not None:
                pt += f"  (dAIC={gp['delta_AIC']:+.1f})"
            ctk.CTkLabel(outer, text=pt,
                         font=ctk.CTkFont(size=9, weight="bold"),
                         text_color="#AAFFAA").pack(pady=(2, 6))

        # Two toggle switches
        tf = ctk.CTkFrame(outer, fg_color="transparent")
        tf.pack(fill="x", padx=8, pady=(0, 6))

        ctk.CTkSwitch(
            tf,
            text="Use global k  (instead of per-galaxy)",
            variable=self.use_global_model,
            command=self._on_toggle,
            font=ctk.CTkFont(size=10), width=40
        ).pack(anchor="w", padx=4, pady=2)

        ctk.CTkSwitch(
            tf,
            text="Power-law k(J)  [Model B]",
            variable=self.use_model_B,
            command=self._on_toggle,
            font=ctk.CTkFont(size=10), width=40
        ).pack(anchor="w", padx=4, pady=(0, 4))

    def _on_toggle(self):
        if self.current_galaxy is not None:
            self.update_plot()

    # =========================================================================
    # Original event handlers (unchanged)
    # =========================================================================

    def on_galaxy_click(self, name):
        if self.active_button is not None:
            self.active_button.configure(fg_color="transparent")
        btn = self.galaxy_buttons[name]
        btn.configure(fg_color="#1F6AA5")
        self.active_button  = btn
        self.current_galaxy = name
        self.ml_slider.set(1.0)
        self.ml_value_label.configure(text="1.00")
        self.update_plot()

    def on_ml_change(self, value):
        self.ml_value_label.configure(text=f"{value:.2f}")
        if self.current_galaxy is not None:
            self.update_plot()

    def apply_ctk_theme(self):
        dark      = ctk.get_appearance_mode() == "Dark"
        self.bg   = "#242424" if dark else "#F9F9FA"
        self.fg   = "#EAEAEA" if dark else "#1A1A1A"
        self.grid = "#3A3A3A" if dark else "#D0D0D0"
        self.box  = "#2B2B2B" if dark else "#FFFFFF"
        self.fig.patch.set_facecolor(self.bg)
        self.ax.set_facecolor(self.bg)
        self.ax.tick_params(colors=self.fg)
        self.ax.title.set_color(self.fg)
        self.ax.xaxis.label.set_color(self.fg)
        self.ax.yaxis.label.set_color(self.fg)
        for spine in self.ax.spines.values():
            spine.set_color(self.grid)

    def on_closing(self):
        plt.close("all")
        self.root.quit()
        self.root.destroy()
        sys.exit(0)

    # =========================================================================
    # Original helpers (unchanged)
    # =========================================================================

    def load_dat_file(self, galaxy_name):
        file_path = os.path.join(self.data_folder, f"{galaxy_name}_rotmod.dat")
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)
        r, v_obs, v_gas, v_disk, v_bul = [], [], [], [], []
        with open(file_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                r.append(float(p[0]));     v_obs.append(float(p[1]))
                v_gas.append(float(p[3])); v_disk.append(float(p[4]))
                v_bul.append(float(p[5]))
        return {
            "r":     np.array(r),    "v_obs":  np.array(v_obs),
            "v_gas": np.array(v_gas),"v_disk": np.array(v_disk),
            "v_bul": np.array(v_bul),
        }

    def calculate_J_vis(self, r_kpc, v_vis_kms):
        r_m      = r_kpc     * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        M_vis    = (v_vis_ms**2 * r_m) / constants.G
        M_max    = M_vis.max()
        M_half   = M_max / 2.0
        idx      = np.argmax(M_vis >= M_half)
        if idx == 0 and M_vis[0] < M_half:
            idx = len(M_vis) - 1
        return float(M_vis[idx] * v_vis_ms[idx] * r_m[idx])

    def fit_k(self, r_kpc, v_obs, v_vis_kms, J_vis):
        r_m      = r_kpc     * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S

        def error_function(log_k):
            k = 10 ** log_k
            errors = []
            for i in range(len(r_m)):
                v_extra_sq = r_m[i] * k * J_vis / (r_m[i] + constants.R0_M)**2
                v_pred_sq  = v_vis_ms[i]**2 + v_extra_sq
                v_pred_kms = (np.sqrt(v_pred_sq) if v_pred_sq >= 0
                              else v_vis_ms[i]) / constants.KM_S_TO_M_S
                if v_obs[i] > 0:
                    errors.append(abs(v_pred_kms - v_obs[i]) / v_obs[i])
            return np.mean(errors) if errors else 1e10

        grid_k    = [1e-45, 1e-43, 1e-42, 1e-41, 1e-40, 1e-38, 1e-36]
        best_init = min(grid_k, key=lambda k: error_function(math.log10(k)))
        result    = minimize(
            error_function, x0=math.log10(best_init),
            method="L-BFGS-B", bounds=[(-55, -25)],
            options={"ftol": 1e-10, "maxiter": 5000}
        )
        return 10 ** result.x[0], result.fun

    def calculate_model(self, r_kpc, v_vis_kms, k, J):
        r_m        = r_kpc     * constants.KPC_TO_M
        v_vis_ms   = v_vis_kms * constants.KM_S_TO_M_S
        v_extra_sq = r_m * k * J / (r_m + constants.R0_M)**2
        return np.sqrt(np.maximum(v_vis_ms**2 + v_extra_sq, 0)) / constants.KM_S_TO_M_S

    # =========================================================================
    # update_plot  —  original + model-toggle block
    # =========================================================================

    def update_plot(self):
        try:
            name = self.current_galaxy
            if name is None:
                return

            ml_scale    = self.ml_slider.get()
            params_orig = self.df_params[self.df_params["galaxy_name"] == name].iloc[0]
            k_orig      = float(params_orig["k_optimal"])
            J_orig      = float(params_orig["J_new"])
            err_orig    = float(params_orig["mean_error_pct"])

            curve = self.load_dat_file(name)

            v_disk_scaled = np.sqrt(ml_scale) * curve["v_disk"]
            v_bul_scaled  = np.sqrt(ml_scale) * curve["v_bul"]
            v_vis         = np.sqrt(curve["v_gas"]**2
                                    + v_disk_scaled**2
                                    + v_bul_scaled**2)

            J_new = self.calculate_J_vis(curve["r"], v_vis)

            # ── Decide which k to use ─────────────────────────────────
            use_global = self.use_global_model.get()
            use_B      = self.use_model_B.get()
            gp         = self.gp
            model_tag  = None   # None  =>  original title/label format

            if use_global and use_B and gp["model_B_k0"] is not None and J_new > 0:
                k_new     = max(gp["model_B_k0"] * (J_new ** gp["model_B_beta"]), 1e-60)
                model_tag = f"Global B  k0*J^b  (b={gp['model_B_beta']:.3f})"

            elif use_global and not use_B and gp["model_A_k0"] is not None:
                k_new     = gp["model_A_k0"]
                model_tag = "Global A  (universal k0)"

            else:
                # Original behaviour: fit k for this galaxy at this M/L
                k_new, _ = self.fit_k(curve["r"], curve["v_obs"], v_vis, J_new)
            # ─────────────────────────────────────────────────────────

            v_model = self.calculate_model(curve["r"], v_vis, k_new, J_new)

            # ── Plot ──────────────────────────────────────────────────
            self.ax.clear()
            self.apply_ctk_theme()

            self.ax.scatter(
                curve["r"], curve["v_obs"],
                s=70, color=self.fg, edgecolors=self.grid,
                linewidths=0.8, alpha=0.9,
                label="Observed (SPARC)", zorder=3
            )

            # Task 6: label is "Baryonic Newtonian" (was "Visible Matter")
            self.ax.plot(
                curve["r"], v_vis,
                linestyle="--", color="#4C72B0",
                linewidth=6.5, alpha=0.9,
                label=f"Baryonic Newtonian (M/L x {ml_scale:.2f})",
                zorder=1
            )

            spin_label = ("Spin-Coupling Model" if model_tag is None
                          else f"Spin-Coupling  [{model_tag}]")
            self.ax.plot(
                curve["r"], v_model,
                color="#DD8452", linewidth=3, alpha=1.0,
                label=spin_label, zorder=2
            )

            all_v = np.concatenate([curve["v_obs"], v_vis, v_model])
            self.ax.set_ylim(max(0, all_v.min()*0.9), all_v.max()*1.1)

            # Error for display
            err_display = float(np.mean(
                np.abs((v_model - curve["v_obs"])
                       / np.where(curve["v_obs"] > 0, curve["v_obs"], 1.0))
            ) * 100)

            k_change_pct = ((k_new - k_orig) / k_orig * 100) if k_orig != 0 else 0.0

            if model_tag:
                self.ax.set_title(
                    f"{name}  |  {model_tag}  |  "
                    f"M/L:{ml_scale:.2f}  Error:{err_display:.2f}%  k={k_new:.2e}"
                )
            else:
                self.ax.set_title(
                    f"{name} | M/L Scale: {ml_scale:.2f} | "
                    f"Error: {err_display:.2f}% | k = {k_new:.2e}"
                )

            self.ax.set_xlabel("Radius (kpc)")
            self.ax.set_ylabel("Velocity (km/s)")
            self.ax.grid(True, linestyle=":", color=self.grid, alpha=0.6)
            self.ax.legend(
                frameon=True, fancybox=True,
                facecolor=self.box, edgecolor=self.grid,
                labelcolor=self.fg, framealpha=0.9
            )

            # Info box
            info = (
                f"M/L Scale: {ml_scale:.2f}\n\n"
                f"Original (per-galaxy, M/L=1.0):\n"
                f"  k = {k_orig:.2e}\n"
                f"  J = {J_orig:.2e}\n"
                f"  Error = {err_orig:.2f}%\n\n"
                f"Current:\n"
                f"  k = {k_new:.2e}\n"
                f"  J = {J_new:.2e}\n"
                f"  Error = {err_display:.2f}%\n\n"
                f"Dk = {k_change_pct:+.1f}%"
            )
            self.ax.text(
                0.02, 0.98, info,
                transform=self.ax.transAxes, fontsize=8.5, color=self.fg,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.4",
                          facecolor=self.box, edgecolor=self.grid, alpha=0.95)
            )

            # Bottom label
            if model_tag:
                self.k_comparison_label.configure(
                    text=f"{model_tag}  |  k={k_new:.2e}  |  Dk vs per-galaxy: {k_change_pct:+.1f}%"
                )
            elif ml_scale < 1.0:
                self.k_comparison_label.configure(
                    text=(f"Lower M/L ({ml_scale:.2f}x) -> Less stellar mass -> "
                          f"k changed by {k_change_pct:+.1f}% to compensate")
                )
            else:
                self.k_comparison_label.configure(
                    text=f"Original SPARC M/L (1.0x) -> k = {k_orig:.2e}"
                )

            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", f"Error updating plot:\n{str(e)}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    root = ctk.CTk()
    app  = GalaxyExplorer(root, RESULTS_CSV, DATA_FOLDER, GLOBAL_SUMMARY)
    root.mainloop()