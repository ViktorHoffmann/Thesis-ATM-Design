import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sympy import Symbol, solve

# === CONFIGURATION ===
def load_json(filename):
    try:
        with open(filename, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        print(f"{filename} not found or invalid. Exiting.")
        exit()

setup = load_json("setup.json")
result = load_json("result.json")

DEVELOPMENT_MODE = setup["devmode"]
avionics_power = setup["avionics_power"]

# --- NEW: read new result keys with sensible fallbacks ---
hybrid_capacity_sim = result.get("hybrid_capacity_sim", None)
hybrid_capacity_nu  = result.get("hybrid_capacity_nu",  None)
pcm_capacity        = result.get("pcm_capacity", None)

if hybrid_capacity_sim is None:
    print("hybrid_capacity_sim (or legacy hybrid_capacity_target) missing in result.json. Exiting.")
    exit()
if hybrid_capacity_nu is None:
    print("hybrid_capacity_nu missing in result.json. Proceeding without Nusselt-based hybrid solution.")

L_target = float(np.sqrt(result["hybrid_radiator_area"]))

# Set global matplotlib style
mpl.rcParams.update({
    "figure.figsize": (4.9, 3.5),
    "font.size": 11.0,
    "font.family": "serif",
    "font.serif": ["cmr10"],
    "axes.titlesize": "medium",
    "figure.titlesize": "medium",
    "text.usetex": not DEVELOPMENT_MODE,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}[=v2]"
})

# === CONSTANTS ===
rho_alu = 2700     # aluminium density [kg*m^-3]
rho_pcm = 788      # pcm density [kg*m^-3]
h       = 240998.9 # pcm latent heat [J*kg^-1]
F       = 0.1      # void fraction
t       = 0.001    # wall thickness [m]

def total_mass(L, H): # pcm mass including case and fins
    return (rho_alu * (L**2 * H - (L - 2*t)**2 * (H - 2*t))
            + (F * rho_alu + (1 - F) * rho_pcm) * (L - 2*t)**2 * (H - 2*t)) 

def total_heat(L, H): # pcm latent heat capacity
    #...#
    pcm_heat  = (1 - F) * rho_pcm * (L - 2*t)**2 * (H - 2*t) * h
    return pcm_heat

# === GRID CALCULATIONS ===
L_vals, H_vals = np.meshgrid(np.linspace(2*h, 0.1, 100), np.linspace(2*h, 0.05, 100))
m_ges = total_mass(L_vals, H_vals)
Q_ges = total_heat(L_vals, H_vals)

# === PLOTTING ===
def plot_contour(X, Y, Z, title, label, filename):
    fig = plt.figure(figsize=(8, 6))
    contour = plt.contour(X, Y, Z, levels=15, colors='black')
    plt.clabel(contour, inline=True, levels=contour.levels[::2], fontsize=10, fmt=label)
    plt.xlabel("Länge [m]")
    plt.ylabel("Höhe [m]")
    plt.grid(True)
    if DEVELOPMENT_MODE:
        plt.title(title)
    else:
        fig.savefig(filename, bbox_inches="tight")

plot_contour(L_vals, H_vals, m_ges, "PCM Masse nach Länge und Höhe", '%1.3f [kg]', "pcm_mass.pdf")
plot_contour(L_vals, H_vals, Q_ges, "PCM Wärmekapazität nach Länge und Höhe", '%1.0f [J]', "pcm_heat_capacity.pdf")

# === SOLUTION HELPERS ===
def solve_positive_real(expr, sym):
    sols = solve(expr, sym)
    sols = [s for s in sols if getattr(s, "is_real", False) and float(s) > 0]
    if not sols:
        raise ValueError(f"No positive real solution for {sym}.")
    return float(sols[0])

def solve_hybrid_H_for(capacity):
    H = Symbol('H')
    eq = total_heat(L_target, H) - capacity
    return solve_positive_real(eq, H)

def solve_mass(L, H):
    # symbolic wrapper kept for consistency; returns float
    return float(total_mass(L, H))

def solve_normal_L_for(capacity):
    L = Symbol('L')
    eq = total_heat(L, L) - capacity
    return solve_positive_real(eq, L)

# === HYBRID SOLUTIONS ===
# From simulation fit
hybrid_H_sim = solve_hybrid_H_for(hybrid_capacity_sim)
hybrid_m_sim = solve_mass(L_target, hybrid_H_sim)

# From Nusselt surplus (if present)
if hybrid_capacity_nu is not None:
    hybrid_H_nu = solve_hybrid_H_for(hybrid_capacity_nu)
    hybrid_m_nu = solve_mass(L_target, hybrid_H_nu)
else:
    hybrid_H_nu = None
    hybrid_m_nu = None

# === NORMAL (no radiator) SOLUTION ===
normal_L_solution = solve_normal_L_for(pcm_capacity)
normal_m_solution = solve_mass(normal_L_solution, normal_L_solution)

# === OUTPUT ===
result.update({
    "hybrid_H_sim": hybrid_H_sim,
    "hybrid_m_sim": hybrid_m_sim,
    "hybrid_H_nu": hybrid_H_nu,
    "hybrid_m_nu": hybrid_m_nu,
    "normal_L_solution": normal_L_solution,
    "normal_m_solution": normal_m_solution
})

with open("result.json", "w") as f:
    json.dump(result, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()