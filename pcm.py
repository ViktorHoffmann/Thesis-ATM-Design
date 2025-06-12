# This Program calculates the pcm mass for given dimensions
#
# setup.json: devmode ------------------------>
# result.json: pcm_capacity, radiator_area ---> pcm.py: H_target, m_target ---> result.json

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sympy.solvers import solve     # sympy needed for solving equations for target values
from sympy import Symbol



# === CONFIGURATION ===
# load data from setup.json and result.json
try:
    with open("setup.json", "r") as f:
        data = json.load(f)
except FileNotFoundError:
    print("setup.json not found. Exiting.")
    exit()
except json.JSONDecodeError:
    print("setup.json is empty of invalid. Exiting.")
    exit()

DEVELOPMENT_MODE = data["devmode"]

try:
    with open("result.json", "r") as f:
        data = json.load(f)
except FileNotFoundError:
    print("setup.json not found. Exiting.")
    exit()
except json.JSONDecodeError:
    print("setup.json is empty or invalid. Exiting.")
    exit()

pcm_capacity_target = data["pcm_capacity_target"]
L_target = np.sqrt(data["radiator_area_target"])

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



# === Equations ===
rho_alu = 2700      #kg/m^3
cp_alu = 896        #J/kgK
rho_pcm = 788       #kg/m^3
cp_pcm = 2950       #J/kgK
dT = 30             #K
dH = 266000         #J/kg
F = 0.1             #void fraction
h = 0.001           #wall thickness [m]

# Create a grid of x and y values
L, H = np.meshgrid(np.linspace(2*h, 0.1, 100), np.linspace(2*h, 0.05, 100))

# Total mass [kg]
m_ges=rho_alu*(np.square(L)*H-np.square(L-2*h)*(H-2*h))+(F*rho_alu+(1-F)*rho_pcm)*np.square(L-2*h)*(H-2*h)

# Toatal Heat [J]
Q_ges=rho_alu*(np.square(L)*H-np.square(L-2*h)*(H-2*h))*cp_alu*dT+F*rho_alu*(np.square(L-2*h)*(H-2*h))*cp_alu*dT+(1-F)*rho_pcm*np.square(L-2*h)*(H-2*h)*(cp_pcm*dT+dH)

# Specific Heat [J/kg]
q_ges = Q_ges/m_ges

# Calculate H for pcm_capacity_target and L_target here



# === Plotting ===
# Mass Plot
fig1 = plt.figure(figsize=(8, 6))
contour = plt.contour(L, H, m_ges, levels=20, colors='black')
#plt.colorbar(contour, label='Gewicht [kg]')    # colorful contours
plt.clabel(contour, inline=True, levels=contour.levels[::2], fontsize=10, fmt='%1.3f [kg]') # levels makes every second contour have a label
plt.xlabel("Länge [m]")
plt.ylabel("Höhe [m]")
plt.grid(True)

if DEVELOPMENT_MODE:
    plt.title("PCM Masse nach Länge und Höhe")

if not DEVELOPMENT_MODE:
    fig1.savefig("pcm_mass.pdf", bbox_inches="tight")

# Heat Plot
fig2 = plt.figure(figsize=(8, 6))
contour = plt.contour(L, H, Q_ges, levels=20, colors='black')
#plt.colorbar(contour, label='Wärme [J]')   # colorful contours
plt.clabel(contour, inline=True, levels=contour.levels[::2], fontsize=10, fmt='%1.0f [J]') # levels makes every second contour have a label
plt.xlabel("Länge [m]")
plt.ylabel("Höhe [m]")
plt.grid(True)

if DEVELOPMENT_MODE:
    plt.title("PCM Wärmekapazität nach Länge und Höhe")

if not DEVELOPMENT_MODE:
    fig2.savefig("pcm_heat_capacity.pdf", bbox_inches="tight")



# === Target result ===
# Solution calculation needs to be down here because sympy.solve needs symbolic variable that is in the way for plotting
H = Symbol('H')
H_solutions = solve(rho_alu*((L_target**2)*H-((L_target-2*h)**2)*(H-2*h))*cp_alu*dT+F*rho_alu*(((L_target-2*h)**2)*(H-2*h))*cp_alu*dT+(1-F)*rho_pcm*((L_target-2*h)**2)*(H-2*h)*(cp_pcm*dT+dH)-pcm_capacity_target, H)
H_target = H_solutions[0]

m = Symbol('m')
m_solutions = solve(rho_alu*((L_target**2)*H_target-((L_target-2*h)**2)*(H_target-2*h))+(F*rho_alu+(1-F)*rho_pcm)*((L_target-2*h)**2)*(H_target-2*h)-m, m)
m_target = m_solutions[0]

# Append new data into result.json
with open("result.json", "r") as f:
    data = json.load(f)
data["H_target"] = float(H_target)
data["m_target"] = float(m_target)
with open("result.json", "w") as f:
    json.dump(data, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()