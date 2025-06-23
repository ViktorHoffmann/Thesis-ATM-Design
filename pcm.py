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
pcm_capacity_target = result["pcm_capacity_target"]
L_target = np.sqrt(result["hybrid_radiator_area"])

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
rho_alu = 2700  # [kg*m^-3]
cp_alu = 896    # [J*kg^-1*K^-1]
rho_pcm = 788   # [kg*m^-3]
cp_pcm = 2950   # [J*kg^-1*kg^-1]
dT = 0.001      # [K]
dH = 266000     # [J*kg^-1]
F = 0.1         # voidfraction
h = 0.001       # wall thickness

def total_mass(L, H):
    return rho_alu * (L**2 * H - (L - 2*h)**2 * (H - 2*h)) + (F * rho_alu + (1 - F) * rho_pcm) * (L - 2*h)**2 * (H - 2*h)

def total_heat(L, H):
    alu_heat = rho_alu * (L**2 * H - (L - 2*h)**2 * (H - 2*h)) * cp_alu * dT
    alu_void = F * rho_alu * (L - 2*h)**2 * (H - 2*h) * cp_alu * dT
    pcm_heat = (1 - F) * rho_pcm * (L - 2*h)**2 * (H - 2*h) * (cp_pcm * dT + dH)
    return alu_heat + alu_void + pcm_heat

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

# === SOLUTION CALCULATION ===
def solve_hybrid_H():
    H = Symbol('H')
    eq = total_heat(L_target, H) - pcm_capacity_target
    return solve(eq, H)[0]

def solve_mass(L, H):
    m = Symbol('m')
    eq = total_mass(L, H) - m
    return solve(eq, m)[0]

def solve_normal_L():
    L = Symbol('L')
    eq = total_heat(L, L) - avionics_power * 1200
    sols = solve(eq, L)
    return [s for s in sols if s.is_real and s > 0][0]

# Hybrid
hybrid_H_target = solve_hybrid_H()
hybrid_m_target = solve_mass(L_target, hybrid_H_target)

# Normal
normal_L_solution = solve_normal_L()
normal_m_solution = solve_mass(normal_L_solution, normal_L_solution)

# === OUTPUT ===
result.update({
    "hybrid_H_target": float(hybrid_H_target),
    "hybrid_m_target": float(hybrid_m_target),
    "normal_L_solution": float(normal_L_solution),
    "normal_m_solution": float(normal_m_solution)
})

with open("result.json", "w") as f:
    json.dump(result, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()