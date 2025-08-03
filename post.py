import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import norm
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset

try:
    with open("setup.json", "r") as f:
        data = json.load(f)
except FileNotFoundError:
    print("setup.json not found. Exiting.")
    exit()
except json.JSONDecodeError:
    print("setup.json is empty or invalid. Exiting.")
    exit()

DEVELOPMENT_MODE = data["devmode"]

heatfluxCompare = "heatfluxCompare.csv"
try:
    raw = np.genfromtxt(heatfluxCompare, delimiter=';', names=True, encoding='utf-8')
except IOError:
    print("heatfluxCompare not found. Exiting.")
    exit()
except ValueError:
    print("Error parsing heatfluxCompare. Exiting.")
    exit()

required_fields = [
    'x', 'maxQ', 'maxQm10', 'maxQp10', 'maxQp20'
]
missing_fields = [f for f in required_fields if f not in raw.dtype.names]
if missing_fields:
    print(f"Missing fields in CSV: {missing_fields}. Exiting.")
    exit()

# === FLIGHT DATA ===
position = raw['x']  # [m]
maxQHeatflux = raw['maxQ']  # [W/m^2]
maxQm10Heatflux = raw['maxQm10']  # [W/m^2]
maxQp10Heatflux = raw['maxQp10']  # [W/m^2]
maxQp20Heatflux = raw['maxQp20']  # [W/m^2]

# === MATPLOTLIB STYLE ===
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

x_values = position

# hull heatflux plots
fig1, ax1 = plt.subplots(constrained_layout=True)

ax1.plot(x_values, maxQHeatflux, color='dimgray', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ}}$')
ax1.plot(x_values, maxQm10Heatflux, color='black', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ-10}}$')
ax1.plot(x_values, maxQp10Heatflux, color='blue', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ+10}}$')
ax1.plot(x_values, maxQp20Heatflux, color='red', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ+20}}$')

ax1.set_xlabel('Position [m]')
ax1.set_ylabel(r'Spezifischer Wärmestrom [$W/m^2$]')
ax1.legend()
if DEVELOPMENT_MODE:
    ax1.set_title("spezifischer Wärmestrom Hülle")
else:
    fig1.savefig("maxQ_compare_heatflux.pdf", bbox_inches="tight")

# eicosane plot
# Constants from UDF and Fluent setup
Cps = 2132.4      # J/kg·K (solid)
Cpl = 2350.05     # J/kg·K (liquid)
Ts = 309.0        # Solidus temperature [K]
Tl = 311.0        # Liquidus temperature [K]
L = 240998.9      # Latent heat [J/kg]
Tc = (Ts + Tl) / 2  # Center of phase change

# Temperature range
T = np.linspace(300, 320, 2000)

# Define Gaussian spike: ensure area = L
width = (Tl - Ts) / 4  # Std dev
gaussian = norm(loc=Tc, scale=width)
spike_area = L
# Estimate the height needed so total spike area integrates to L
delta_T = T[1] - T[0]
spike_range = (T >= Ts) & (T <= Tl)
pdf_vals = gaussian.pdf(T[spike_range])
normalization_factor = np.sum(pdf_vals * delta_T)
spike_height = spike_area / normalization_factor

# eicosan effective Cp values
Cp = []
for temp in T:
    if temp < Ts:
        Cp_val = Cps
    elif temp > Tl:
        Cp_val = Cpl
    else:
        linear_blend = ((1 - (temp - Ts)/(Tl - Ts)) * Cps + ((temp - Ts)/(Tl - Ts)) * Cpl)
        spike = spike_height * gaussian.pdf(temp)
        Cp_val = linear_blend + spike
    Cp.append(Cp_val)

# effective cp vs T Plot
fig2, ax2 = plt.subplots(constrained_layout=True)
ax2.plot(T, Cp, color='darkred', label='Effective Cp')
ax2.axvspan(Ts, Tl, color='orange', alpha=0.2, label='Phasenwechselbereich (309–311 K)')
ax2.set_xlabel('Temperatur [K]')
ax2.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax2.grid(True)
ax2.legend()

if DEVELOPMENT_MODE:
    ax2.set_title('Eicosan effektives Cp')
else:
    fig2.savefig("eicosane_cpvst_spike.pdf", bbox_inches="tight")

#
Cp_udf = []
for temp in T:
    if temp < Ts:
        Cp_val = Cps
    elif temp > Tl:
        Cp_val = Cpl
    else:
        Gama = (temp - Ts) / (Tl - Ts)
        numerator = (1 - Gama) * 910.0 * Cps + Gama * 769.0 * Cpl
        denominator = (1 - Gama) * 910.0 + Gama * 769.0
        Cp_val = numerator / denominator
    Cp_udf.append(Cp_val)

# specific cp vs T Plot
fig3, ax3 = plt.subplots(constrained_layout=True)
ax3.plot(T, Cp_udf, color='darkgreen', label='Specific Cp')
ax3.axvspan(Ts, Tl, color='orange', alpha=0.2, label='Phasenwechselbereich (309–311 K)')
ax3.set_xlabel('Temperatur [K]')
ax3.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax3.grid(True)
ax3.legend()

if DEVELOPMENT_MODE:
    ax3.set_title('Eicosan Cp aus UDF')
else:
    fig3.savefig("eicosane_cpvst_udf.pdf", bbox_inches="tight")

if DEVELOPMENT_MODE:
    plt.show()