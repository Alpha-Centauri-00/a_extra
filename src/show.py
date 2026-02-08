import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.optimize import minimize
from scipy.interpolate import interp1d

from tools.config import constants

# --- Physics Engine & Constants ---
# G = 6.674e-11
# C_LIGHT = 3.0e8
# KPC_TO_M = 3.086e19
# KM_S_TO_M_S = 1000.0
# R0_M = 0.5 * KPC_TO_M

def run_model_fit(r_kpc, v_obs, v_gas, v_disk, v_bulge):
    r_m = r_kpc * constants.KPC_TO_M
    v_vis = np.sqrt(v_gas**2 + v_disk**2 + v_bulge**2)
    v_vis_m_s = v_vis * constants.KM_S_TO_M_S
    
    # Calculate J_vis
    M_vis = (v_vis_m_s**2 * r_m) / constants.G
    M_half = M_vis[-1] / 2.0
    f_r = interp1d(M_vis, r_m, fill_value="extrapolate")
    r_disk_m = f_r(M_half)
    f_v = interp1d(r_m, v_vis_m_s, fill_value="extrapolate")
    J_vis = M_half * f_v(r_disk_m) * r_disk_m

    def objective(log_k):
        k_val = 10**log_k
        v_extra_sq = r_m * k_val * J_vis / (r_m + constants.R0_M)**2
        v_pred = np.sqrt(np.maximum(v_vis_m_s**2 + v_extra_sq, 0)) / constants.KM_S_TO_M_S
        return np.mean(np.abs((v_pred - v_obs) / v_obs))

    res = minimize(objective, x0=-40, bounds=[(-60, -20)])
    best_k = 10**res.x[0]
    
    v_extra_sq_final = r_m * best_k * J_vis / (r_m + constants.R0_M)**2
    v_pred_final = np.sqrt(np.maximum(v_vis_m_s**2 + v_extra_sq_final, 0)) / constants.KM_S_TO_M_S
    
    return v_vis, v_pred_final, best_k, res.fun, J_vis

# --- Data Repository ---
# Add your 175 galaxies here!
galaxy_data = {
    "CamB": {
        "r": np.array([0.16, 0.41, 0.57, 0.73, 0.90, 1.06, 1.22, 1.47, 1.79]),
        "v_obs": np.array([1.99, 4.84, 6.79, 8.87, 10.90, 12.90, 14.70, 16.80, 20.10]),
        "v_gas": np.array([1.86, 4.24, 5.61, 6.77, 7.77, 8.44, 8.64, 8.08, 6.91]),
        "v_disk": np.array([3.75, 9.47, 11.76, 13.72, 14.80, 15.24, 15.11, 15.90, 14.91]),
        "v_bulge": np.zeros(9)
    }
}

# --- GUI Application ---
class GalaxyApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Spin-Mediated Acceleration Explorer")
        
        # Sidebar for List
        self.sidebar = tk.Frame(root, width=200, bg="#f0f0f0")
        self.sidebar.pack(side="left", fill="y", padx=10, pady=10)
        
        tk.Label(self.sidebar, text="Select Galaxy", font=("Arial", 12, "bold")).pack()
        self.listbox = tk.Listbox(self.sidebar, font=("Arial", 10))
        self.listbox.pack(fill="both", expand=True)
        for name in galaxy_data.keys():
            self.listbox.insert(tk.END, name)
        self.listbox.bind("<<ListboxSelect>>", self.on_select)
        
        # Plot Area
        self.plot_frame = tk.Frame(root)
        self.plot_frame.pack(side="right", fill="both", expand=True)
        
        self.fig, self.ax = plt.subplots(figsize=(7, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def on_select(self, event):
        selection = self.listbox.curselection()
        if not selection: return
        name = self.listbox.get(selection[0])
        data = galaxy_data[name]
        
        v_vis, v_pred, k, err, J = run_model_fit(data["r"], data["v_obs"], data["v_gas"], data["v_disk"], data["v_bulge"])
        
        self.ax.clear()
        self.ax.scatter(data["r"], data["v_obs"], color='black', label='Observed (SPARC)')
        self.ax.plot(data["r"], v_vis, '--', color='blue', label='Baryonic (Newtonian)')
        self.ax.plot(data["r"], v_pred, color='red', linewidth=2, label=f'Model (k={k:.1e})')
        
        self.ax.set_title(f"Galaxy: {name} | Fit Error: {err*100:.2f}%")
        self.ax.set_xlabel("Radius (kpc)")
        self.ax.set_ylabel("Velocity (km/s)")
        self.ax.legend()
        self.ax.grid(alpha=0.3)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = GalaxyApp(root)
    root.mainloop()