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
DATA_FOLDER = r"data\Rotmod_LTG"
RESULTS_CSV = r"results\batch_fitting_results_continuous_k.csv"

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


class GalaxyExplorer:
    def __init__(self, root, results_csv, data_folder):
        
        self.galaxy_buttons = {}
        self.active_button = None
        self.current_galaxy = None
        self.current_curve = None

        self.root = root
        try:
            root.iconbitmap("data/Untitled.ico")
        except:
            pass
        
        self.root.title("SPARC Spin-Coupling Explorer - M/L Sensitivity Test")
        self.root.geometry("1400x900")
        self.data_folder = data_folder

        try:
            self.df_params = pd.read_csv(results_csv)
            self.galaxy_list = sorted(self.df_params["galaxy_name"].unique().tolist())
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.root.destroy()
            return

        # Sidebar
        self.sidebar = ctk.CTkFrame(root, width=260, fg_color="#242424")
        self.sidebar.pack(side="left", fill="y", padx=10, pady=10)

        ctk.CTkLabel(
            self.sidebar,
            text="Galaxies",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(pady=(10, 5))

        self.list_frame = ctk.CTkScrollableFrame(self.sidebar)
        self.list_frame.pack(fill="both", expand=True, padx=5, pady=5)

        for name in self.galaxy_list:
            btn = ctk.CTkButton(
                self.list_frame,
                text=name,
                anchor="w",
                fg_color="transparent",
                hover_color="#3A3A3A",
                command=lambda n=name: self.on_galaxy_click(n)
            )
            btn.pack(fill="x", padx=5, pady=2)
            self.galaxy_buttons[name] = btn

        # Main content area
        self.main_frame = ctk.CTkFrame(root)
        self.main_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        # Plot frame
        self.plot_frame = ctk.CTkFrame(self.main_frame)
        self.plot_frame.pack(fill="both", expand=True)

        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.toolbar.pack_forget()

        btn_frame = ctk.CTkFrame(self.plot_frame, fg_color="#242424")
        btn_frame.pack(side="bottom", fill="x")

        controls = ctk.CTkFrame(btn_frame, fg_color="transparent")
        controls.pack(expand=True, pady=6)

        ctk.CTkButton(controls, text="Home", command=self.toolbar.home).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Pan", command=self.toolbar.pan).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Zoom", command=self.toolbar.zoom).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Save", command=self.toolbar.save_figure).pack(side="left", padx=5)

        # M/L Slider Control Panel
        self.ml_control_frame = ctk.CTkFrame(self.main_frame, fg_color="#1a1a1a")
        self.ml_control_frame.pack(fill="x", pady=(10, 0))

        # Title
        ctk.CTkLabel(
            self.ml_control_frame,
            text="Stellar Mass-to-Light Ratio (M/L) Sensitivity Test",
            font=ctk.CTkFont(size=14, weight="bold")
        ).pack(pady=(10, 5))

        # Explanation
        ctk.CTkLabel(
            self.ml_control_frame,
            text="Lower M/L = Less stellar mass → Blue line drops → Spin-coupling needs larger k",
            font=ctk.CTkFont(size=11),
            text_color="#AAAAAA"
        ).pack()

        # Slider frame
        slider_frame = ctk.CTkFrame(self.ml_control_frame, fg_color="transparent")
        slider_frame.pack(fill="x", padx=20, pady=10)

        # Label on left
        ctk.CTkLabel(
            slider_frame,
            text="M/L Scale:",
            font=ctk.CTkFont(size=12, weight="bold")
        ).pack(side="left", padx=(0, 10))

        # Slider
        self.ml_slider = ctk.CTkSlider(
            slider_frame,
            from_=0.2,
            to=1.0,
            number_of_steps=40,
            command=self.on_ml_change
        )
        self.ml_slider.set(1.0)  # Default = original SPARC M/L
        self.ml_slider.pack(side="left", fill="x", expand=True, padx=5)

        # Value display
        self.ml_value_label = ctk.CTkLabel(
            slider_frame,
            text="1.00",
            font=ctk.CTkFont(size=12, weight="bold"),
            width=50
        )
        self.ml_value_label.pack(side="left", padx=(10, 0))

        # Info display frame
        info_frame = ctk.CTkFrame(self.ml_control_frame, fg_color="transparent")
        info_frame.pack(fill="x", padx=20, pady=(0, 10))

        self.k_comparison_label = ctk.CTkLabel(
            info_frame,
            text="Select a galaxy to see M/L sensitivity",
            font=ctk.CTkFont(size=11),
            text_color="#CCCCCC"
        )
        self.k_comparison_label.pack()

        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Initial plot
        self.ax.text(
            0.5, 0.5,
            "Select a galaxy from the list",
            ha="center",
            va="center",
            fontsize=14,
            transform=self.ax.transAxes
        )
        self.apply_ctk_theme()
        self.canvas.draw()

    def on_galaxy_click(self, name):
        # Reset previously active button
        if self.active_button is not None:
            self.active_button.configure(fg_color="transparent")

        # Activate new button
        btn = self.galaxy_buttons[name]
        btn.configure(fg_color="#1F6AA5")
        self.active_button = btn

        # Store galaxy name and reset slider
        self.current_galaxy = name
        self.ml_slider.set(1.0)
        self.ml_value_label.configure(text="1.00")
        
        # Update plot
        self.update_plot()

    def on_ml_change(self, value):
        """Called when M/L slider moves."""
        self.ml_value_label.configure(text=f"{value:.2f}")
        
        if self.current_galaxy is not None:
            self.update_plot()

    def apply_ctk_theme(self):
        dark = ctk.get_appearance_mode() == "Dark"
        self.bg = "#242424" if dark else "#F9F9FA"
        self.fg = "#EAEAEA" if dark else "#1A1A1A"
        self.grid = "#3A3A3A" if dark else "#D0D0D0"
        self.box = "#2B2B2B" if dark else "#FFFFFF"

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

    def load_dat_file(self, galaxy_name):
        """Load rotation curve data from .dat file."""
        file_path = os.path.join(self.data_folder, f"{galaxy_name}_rotmod.dat")
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        r, v_obs, v_gas, v_disk, v_bul = [], [], [], [], []

        with open(file_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                r.append(float(p[0]))
                v_obs.append(float(p[1]))
                v_gas.append(float(p[3]))
                v_disk.append(float(p[4]))
                v_bul.append(float(p[5]))

        return {
            "r": np.array(r),
            "v_obs": np.array(v_obs),
            "v_gas": np.array(v_gas),
            "v_disk": np.array(v_disk),
            "v_bul": np.array(v_bul),
        }

    def calculate_J_vis(self, r_kpc, v_vis_kms):
        """Calculate visible angular momentum using half-mass radius method."""
        r_m = r_kpc * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        
        # Calculate enclosed mass
        M_vis = (v_vis_ms**2 * r_m) / constants.G
        
        # Find half-mass radius
        M_max = M_vis[-1]
        M_half = M_max / 2.0
        
        # Find index where M >= M_half
        idx = np.argmax(M_vis >= M_half)
        if idx == 0 and M_vis[0] < M_half:
            idx = len(M_vis) - 1
        
        r_disk_m = r_m[idx]
        M_disk = M_vis[idx]
        v_disk_ms = v_vis_ms[idx]
        
        J_vis = M_disk * v_disk_ms * r_disk_m
        
        return J_vis

    def fit_k(self, r_kpc, v_obs, v_vis_kms, J_vis):
        """Fit k parameter using scipy optimization."""
        r_m = r_kpc * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        
        def error_function(log_k):
            k = 10 ** log_k
            errors = []
            
            for i in range(len(r_m)):
                # Calculate model velocity
                v_extra_sq = r_m[i] * k * J_vis / (r_m[i] + constants.R0_M)**2
                v_pred_sq = v_vis_ms[i]**2 + v_extra_sq
                
                if v_pred_sq < 0:
                    v_pred = v_vis_kms[i]
                else:
                    v_pred = np.sqrt(v_pred_sq) / constants.KM_S_TO_M_S
                
                if v_obs[i] > 0:
                    errors.append(abs(v_pred - v_obs[i]) / v_obs[i])
            
            return np.mean(errors) if errors else 1e10
        
        # Grid search for initial guess
        grid_k = [1e-45, 1e-43, 1e-42, 1e-41, 1e-40, 1e-38, 1e-36]
        best_k_init = min(grid_k, key=lambda k: error_function(math.log10(k)))
        
        # Optimize
        result = minimize(
            error_function,
            x0=math.log10(best_k_init),
            method="L-BFGS-B",
            bounds=[(-55, -25)],
            options={"ftol": 1e-10, "maxiter": 5000}
        )
        
        k_optimal = 10 ** result.x[0]
        error_optimal = result.fun
        
        return k_optimal, error_optimal

    def calculate_model(self, r_kpc, v_vis_kms, k, J):
        """Calculate model velocity curve."""
        r_m = r_kpc * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        v_extra_sq = r_m * k * J / (r_m + constants.R0_M)**2
        v_pred_ms = np.sqrt(np.maximum(v_vis_ms**2 + v_extra_sq, 0))
        return v_pred_ms / constants.KM_S_TO_M_S

    def update_plot(self):
        """Update plot with current M/L scaling."""
        try:
            name = self.current_galaxy
            if name is None:
                return
            
            # Get M/L scale factor
            ml_scale = self.ml_slider.get()
            
            # Get original parameters for comparison
            params_orig = self.df_params[self.df_params["galaxy_name"] == name].iloc[0]
            k_orig = params_orig["k_optimal"]
            J_orig = params_orig["J_new"]
            err_orig = params_orig["mean_error_pct"]

            # Load curve data
            curve = self.load_dat_file(name)
            
            # Apply M/L scaling to stellar components
            # v_stellar^2 ∝ M_stellar ∝ M/L
            # Therefore v_stellar_new = sqrt(ml_scale) * v_stellar_orig
            v_disk_scaled = np.sqrt(ml_scale) * curve["v_disk"]
            v_bul_scaled = np.sqrt(ml_scale) * curve["v_bul"]
            
            # Calculate new v_vis with scaled stellar components
            v_vis = np.sqrt(curve["v_gas"]**2 + v_disk_scaled**2 + v_bul_scaled**2)
            
            # Calculate new J_vis
            J_new = self.calculate_J_vis(curve["r"], v_vis)
            
            # Fit new k
            k_new, err_new = self.fit_k(curve["r"], curve["v_obs"], v_vis, J_new)
            
            # Calculate model curve
            v_model = self.calculate_model(curve["r"], v_vis, k_new, J_new)

            # Update plot
            self.ax.clear()
            self.apply_ctk_theme()

            self.ax.scatter(
                curve["r"], curve["v_obs"],
                s=70, color=self.fg,
                edgecolors=self.grid, linewidths=0.8,
                alpha=0.9, label="Observed (SPARC)", zorder=3
            )

            self.ax.plot(
                curve['r'], v_vis,
                linestyle="--",
                color="#4C72B0",
                linewidth=6.5,
                alpha=0.9,
                label=f"Visible Matter (M/L × {ml_scale:.2f})",
                zorder=1
            )

            self.ax.plot(
                curve['r'], v_model,
                color="#DD8452",
                linewidth=3,
                alpha=1.0,
                label="Spin-Coupling Model",
                zorder=2
            )

            # Auto-scale y-axis
            all_v = np.concatenate([curve["v_obs"], v_vis, v_model])
            self.ax.set_ylim(max(0, all_v.min()*0.9), all_v.max()*1.1)

            # Title with comparison
            k_change_pct = ((k_new - k_orig) / k_orig) * 100
            self.ax.set_title(
                f"{name} | M/L Scale: {ml_scale:.2f} | Error: {err_new*100:.2f}% | k = {k_new:.2e}"
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
                f"Original (M/L=1.0):\n"
                f"  k = {k_orig:.2e}\n"
                f"  J = {J_orig:.2e}\n"
                f"  Error = {err_orig:.2f}%\n\n"
                f"Scaled:\n"
                f"  k = {k_new:.2e}\n"
                f"  J = {J_new:.2e}\n"
                f"  Error = {err_new*100:.2f}%\n\n"
                f"Δk = {k_change_pct:+.1f}%"
            )

            self.info_box = self.ax.text(
                0.02, 0.98, info,
                transform=self.ax.transAxes,
                fontsize=8.5,
                color=self.fg,
                verticalalignment="top",
                bbox=dict(
                    boxstyle="round,pad=0.4",
                    facecolor=self.box,
                    edgecolor=self.grid,
                    alpha=0.95
                )
            )

            # Update bottom label
            if ml_scale < 1.0:
                comparison_text = (
                    f"Lower M/L ({ml_scale:.2f}×) → Less stellar mass → "
                    f"k increased by {k_change_pct:+.1f}% to compensate"
                )
            else:
                comparison_text = f"Original SPARC M/L (1.0×) → k = {k_orig:.2e}"
            
            self.k_comparison_label.configure(text=comparison_text)

            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", f"Error updating plot:\n{str(e)}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    root = ctk.CTk()
    app = GalaxyExplorer(root, RESULTS_CSV, DATA_FOLDER)
    root.mainloop()