# This Program calculates the steady-state equation for a grey body radiator

import numpy as np
import json
from scipy.constants import Stefan_Boltzmann, pi
import matplotlib.pyplot as plt
import matplotlib as mpl



# === configuration ===
# load data from setup.json
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
avionics_power = data["avionics_power"] # W
e = data["emittance"]
a = data["absorptance"]
solar_flux = data["solar_flux"] # W*m^-2
target_temperature = data["target_temperature"]

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
rho_alu = 2700 # kg/m^3
T_target = 273.15+target_temperature # target temperature

# Define temperature in Celsius for plotting, and convert to Kelvin
temp_C = np.linspace(0, 100, 100) # °C
A = np.linspace(0.01*0.01, 0.3*0.3, 100) # Area [m²]

A_grid, T_C_grid = np.meshgrid(A, temp_C)
T_K_grid = T_C_grid + 273.15  # Convert Celsius to Kelvin for calculation

# Stefan–Boltzmann radiation formula
phi_radiation = (A_grid * e * Stefan_Boltzmann * T_K_grid**4) # [W]



# === Target result ===
hybrid_radiator_area = (avionics_power) / ((e * Stefan_Boltzmann * T_target**4)-(solar_flux/2)*a) # [m^2] necessary radiator area to get rid of avionics- and solarheatflux
hybrid_radiator_power = hybrid_radiator_area*e*Stefan_Boltzmann*T_target**4 # [W] total out-heatflux of radiator

# Putting Data into result.json
with open("result.json", "w") as f:
    json.dump({}, f, indent=4)  # clear result.json as this is the first write in the programm chain

data = {
    "hybrid_radiator_area": hybrid_radiator_area,
    "hybrid_radiator_power": hybrid_radiator_power
}

with open("result.json", "w") as f:
    json.dump(data, f, indent=4)



# === Plotting ===
# Plot
fig1 = plt.figure(figsize=(8, 6))
contour = plt.contour(A_grid, T_C_grid, phi_radiation, levels=20, colors='black')
#plt.colorbar(contour, label='Wärme [W]')   # colorful contours
plt.clabel(contour, inline=True, levels=contour.levels[::2], fontsize=10, fmt='%1.1f [W]') # levels makes every second contour have a label
plt.xlabel("Fläche [m²]")
plt.ylabel("Temperatur [°C]")
plt.grid(True)

if DEVELOPMENT_MODE:
    plt.title("Radiator Wärmestrahlung nach Fläche und Temperatur")

# Save it if not in dev mode
if not DEVELOPMENT_MODE:
    fig1.savefig("radiator_leistung.pdf", bbox_inches="tight")

if DEVELOPMENT_MODE:
    plt.show()