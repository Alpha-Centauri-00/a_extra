import pandas as pd
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import sys
from tools.config import constants


# --- Configuration ---
DATA_FOLDER = r"data\Rotmod_LTG"
RESULTS_CSV = r"results\batch_fitting_results_continuous_k.csv"

ctk.set_appearance_mode("dark")  # "dark", "light", "system"
ctk.set_default_color_theme("blue")  # "blue", "dark-blue", "green"

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
        self.root = root
        root.iconbitmap("data/Untitled.ico")
        self.root.title("SPARC Spin-Coupling Explorer")
        self.root.geometry("1200x800")
        self.data_folder = data_folder

        try:
            self.df_params = pd.read_csv(results_csv)
            self.galaxy_list = sorted(self.df_params["galaxy_name"].unique().tolist())
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.root.destroy()
            return

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
            ctk.CTkButton(
                self.list_frame,
                text=name,
                anchor="w",
                command=lambda n=name: self.update_plot(n)
            ).pack(fill="x", padx=5, pady=2)

        self.plot_frame = ctk.CTkFrame(root)
        self.plot_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # toolbar_frame = tk.Frame(self.plot_frame)
        # toolbar_frame.pack(side="bottom", fill="x")
        # self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        # self.toolbar.update()
#######################################################
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.toolbar.pack_forget()  # hide the default toolbar

        btn_frame = ctk.CTkFrame(self.plot_frame,fg_color="#242424")
        btn_frame.pack(side="bottom", fill="x")


        controls = ctk.CTkFrame(btn_frame, fg_color="transparent")
        controls.pack(expand=True, pady=6)

        ctk.CTkButton(controls, text="Home",
                    command=self.toolbar.home).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Pan",
                    command=self.toolbar.pan).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Zoom",
                    command=self.toolbar.zoom).pack(side="left", padx=5)
        ctk.CTkButton(controls, text="Save",
                    command=self.toolbar.save_figure).pack(side="left", padx=5)


        


        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

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

    def calculate_model(self, r_kpc, v_vis_kms, k, J):
        r_m = r_kpc * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        v_extra_sq = r_m * k * J / (r_m + constants.R0_M) ** 2
        v_pred_ms = np.sqrt(np.maximum(v_vis_ms**2 + v_extra_sq, 0))
        return v_pred_ms / constants.KM_S_TO_M_S

    def update_plot(self, name):
        try:
            params = self.df_params[self.df_params["galaxy_name"] == name].iloc[0]
            k, J, err = params["k_optimal"], params["J_new"], params["mean_error_pct"]

            curve = self.load_dat_file(name)
            v_vis = np.sqrt(curve["v_gas"]**2 + curve["v_disk"]**2 + curve["v_bul"]**2)
            v_model = self.calculate_model(curve["r"], v_vis, k, J)

            self.ax.clear()
            self.apply_ctk_theme()

            self.ax.scatter(
                curve["r"], curve["v_obs"],
                s=70, color=self.fg,
                edgecolors=self.grid, linewidths=0.8,
                alpha=0.9, label="Observed", zorder=3
            )

            # self.ax.plot(
            #     curve["r"], v_vis,
            #     "--", color="#6BAED6",
            #     linewidth=2.2, alpha=0.85,
            #     label="Visible Matter"
            # )

            # self.ax.plot(
            #     curve["r"], v_model,
            #     color="#FB6A4A",
            #     linewidth=3, alpha=0.95,
            #     label="Spin-Coupling Model"
            # )

            self.ax.plot(
                curve['r'], v_vis,
                linestyle="--",
                color="#4C72B0",
                linewidth=5.5,        # thickness
                alpha=0.9,            # slightly more opaque
                label="Visible Matter (Newtonian)",
                zorder=1              # ensure it's below the red/orange line
            )

            # Spin-Coupling Model line (red/orange)
            self.ax.plot(
                curve['r'], v_model,
                color="#DD8452",
                linewidth=3,
                alpha=0.95,
                label="Spin-Coupling Model",
                zorder=2              # drawn above the blue line
            )

            all_v = np.concatenate([curve["v_obs"], v_vis, v_model])
            self.ax.set_ylim(max(0, all_v.min()*0.9), all_v.max()*1.1)

            self.ax.set_title(
                f"{name} | Error: {err:.2f}% | k={k:.2e} | J={J:.2e}"
            )
            self.ax.set_xlabel("Radius (kpc)")
            self.ax.set_ylabel("Velocity (km/s)")

            self.ax.grid(True, linestyle=":", color=self.grid, alpha=0.6)

            self.ax.legend(
                frameon=True, fancybox=True,
                facecolor=self.box, edgecolor=self.grid,
                labelcolor=self.fg, framealpha=0.9
            )

            # if hasattr(self, "info_box"):
            #     self.info_box.remove()

            info = (
                f"log₁₀(J) = {np.log10(J):.2f}\n"
                f"log₁₀(k) = {np.log10(k):.2f}\n\n"
                f"Gap = Missing Mass\n"
                f"Red fills the gap!"
            )

            self.info_box = self.ax.text(
                0.02, 0.98, info,
                transform=self.ax.transAxes,
                fontsize=9,
                color=self.fg,  # text color
                verticalalignment="top",
                bbox=dict(
                    boxstyle="round,pad=0.4",
                    facecolor=self.box,   # match CTk background
                    edgecolor=self.grid,  # subtle border
                    alpha=0.95
                )
            )


            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", str(e))


if __name__ == "__main__":
    root = ctk.CTk()
    app = GalaxyExplorer(root, RESULTS_CSV, DATA_FOLDER)
    root.mainloop()
