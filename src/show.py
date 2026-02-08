import pandas as pd
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
from tools.config import constants


# --- Configuration ---
DATA_FOLDER = r"data\Rotmod_LTG"
RESULTS_CSV = r"results\batch_fitting_results_continuous_k.csv"

class GalaxyExplorer:
    def __init__(self, root, results_csv, data_folder):
        self.root = root
        self.root.title("SPARC Spin-Coupling Explorer")
        self.data_folder = data_folder
        
        # 1. Load Results Data
        try:
            self.df_params = pd.read_csv(results_csv)
            # Get galaxy names from the CSV
            self.galaxy_list = sorted(self.df_params['galaxy_name'].unique().tolist())
            print(f"Loaded {len(self.galaxy_list)} galaxies from {results_csv}")
        except Exception as e:
            messagebox.showerror("Error", f"Could not load {results_csv}: {e}")
            self.root.destroy()
            return

        # 2. UI Layout
        self.sidebar = tk.Frame(root, width=250, bg="#f0f0f0")
        self.sidebar.pack(side="left", fill="y", padx=10, pady=10)
        
        tk.Label(self.sidebar, text="Galaxies", font=("Arial", 12, "bold"), bg="#f0f0f0").pack(pady=5)
        
        # Add scrollbar to listbox
        scrollbar = tk.Scrollbar(self.sidebar)
        scrollbar.pack(side="right", fill="y")
        
        self.listbox = tk.Listbox(self.sidebar, font=("Arial", 10), yscrollcommand=scrollbar.set)
        self.listbox.pack(fill="both", expand=True)
        scrollbar.config(command=self.listbox.yview)
        
        for name in self.galaxy_list:
            self.listbox.insert(tk.END, name)
        self.listbox.bind("<<ListboxSelect>>", self.update_plot)
        
        # Plot frame
        self.plot_frame = tk.Frame(root)
        self.plot_frame.pack(side="right", fill="both", expand=True)
        
        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Initial message
        self.ax.text(0.5, 0.5, 'Select a galaxy from the list', 
                    ha='center', va='center', fontsize=14, transform=self.ax.transAxes)
        self.canvas.draw()

    def load_dat_file(self, galaxy_name):
        """Parses the SPARC .dat file format."""
        file_path = os.path.join(self.data_folder, f"{galaxy_name}_rotmod.dat")
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing file: {file_path}")
        
        # Read all lines
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse data (skip comment lines starting with #)
        data_lines = [line.strip() for line in lines if not line.startswith('#') and line.strip()]
        
        # Parse into columns
        r_list = []
        v_obs_list = []
        v_gas_list = []
        v_disk_list = []
        v_bul_list = []
        
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 6:
                r_list.append(float(parts[0]))
                v_obs_list.append(float(parts[1]))
                v_gas_list.append(float(parts[3]))
                v_disk_list.append(float(parts[4]))
                v_bul_list.append(float(parts[5]))
        
        return {
            'r': np.array(r_list),
            'v_obs': np.array(v_obs_list),
            'v_gas': np.array(v_gas_list),
            'v_disk': np.array(v_disk_list),
            'v_bul': np.array(v_bul_list)
        }

    def calculate_model(self, r_kpc, v_vis_kms, k_val, J_val):
        """Calculates the predicted velocity based on the spin-coupling model."""
        r_m = r_kpc * constants.KPC_TO_M
        v_vis_ms = v_vis_kms * constants.KM_S_TO_M_S
        
        # Acceleration extra term: a_extra = k * J / (r + r0)^2
        # v_model^2 = v_vis^2 + r * a_extra
        v_extra_sq = r_m * k_val * J_val / (r_m + constants.R0_M)**2
        
        v_pred_ms = np.sqrt(np.maximum(v_vis_ms**2 + v_extra_sq, 0))
        return v_pred_ms / constants.KM_S_TO_M_S

    def update_plot(self, event):
        selection = self.listbox.curselection()
        if not selection: 
            return
        
        name = self.listbox.get(selection[0])
        
        try:
            # 1. Get fitted parameters from CSV
            params = self.df_params[self.df_params['galaxy_name'] == name].iloc[0]
            k = params['k_optimal']
            J = params['J_new']
            err_pct = params['mean_error_pct']
            
            # 2. Load raw rotation curve from .dat file
            curve = self.load_dat_file(name)
            
            # 3. Calculate visible matter velocity
            v_vis_kms = np.sqrt(curve['v_gas']**2 + curve['v_disk']**2 + curve['v_bul']**2)
            
            # 4. Calculate model prediction
            v_model_kms = self.calculate_model(curve['r'], v_vis_kms, k, J)
            
            # 5. Refresh plot
            self.ax.clear()
            
            # Plot components
            self.ax.scatter(curve['r'], curve['v_obs'], color='black', s=50, 
                          label='Observed (SPARC)', zorder=3)
            self.ax.plot(curve['r'], v_vis_kms, '--', color='blue', linewidth=2,
                       label='Visible Matter (Newtonian)', alpha=0.7)
            self.ax.plot(curve['r'], v_model_kms, color='red', linewidth=2.5, 
                       label=f'Spin-Coupling Model', zorder=2)
            
            # Styling
            self.ax.set_title(f"{name} | Error: {err_pct:.2f}% | k = {k:.2e} | J = {J:.2e}", 
                            fontsize=12, fontweight='bold')
            self.ax.set_xlabel("Radius (kpc)", fontsize=11, fontweight='bold')
            self.ax.set_ylabel("Velocity (km/s)", fontsize=11, fontweight='bold')
            self.ax.legend(fontsize=10, loc='best')
            self.ax.grid(True, linestyle=':', alpha=0.4)
            
            # Add info text box
            info_text = f"log₁₀(J) = {np.log10(J):.2f}\nlog₁₀(k) = {np.log10(k):.2f}"
            self.ax.text(0.02, 0.98, info_text, transform=self.ax.transAxes,
                       fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            self.canvas.draw()
            
        except FileNotFoundError as e:
            messagebox.showerror("File Not Found", str(e))
        except KeyError as e:
            messagebox.showerror("Data Error", f"Missing column in CSV: {e}")
        except Exception as e:
            messagebox.showerror("Plot Error", f"Error loading data for {name}:\n{e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    root = tk.Tk()
    app = GalaxyExplorer(root, RESULTS_CSV, DATA_FOLDER)
    root.mainloop()