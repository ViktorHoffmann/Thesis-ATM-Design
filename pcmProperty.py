import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl

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

# === Constants from UDF and Fluent setup ===
Cps = 2132.4      # J/kg·K (solid)
Cpl = 2350.05     # J/kg·K (liquid)
Ts = 309.0        # Solidus temperature [K]
Tl = 311.0        # Liquidus temperature [K]
L = 240998.9      # Latent heat [J/kg]
Tc = (Ts + Tl) / 2  # Center of phase change

# === Temperature range ===
T = np.linspace(300, 320, 2000)

# === Define parabolic spike: ensure area = L ===
# f(T) = a * (T - Ts) * (Tl - T) over [Ts, Tl]
a = 6 * L / (Tl - Ts)**3  # From definite integral over the interval

# === Eicosane effective Cp values with parabolic spike ===
Cp = []
for temp in T:
    if temp < Ts:
        Cp_val = Cps
    elif temp > Tl:
        Cp_val = Cpl
    else:
        linear_blend = ((1 - (temp - Ts)/(Tl - Ts)) * Cps + ((temp - Ts)/(Tl - Ts)) * Cpl)
        spike = a * (temp - Ts) * (Tl - temp)
        Cp_val = linear_blend + spike
    Cp.append(Cp_val)

# === Effective Cp vs T Plot ===
fig1, ax1 = plt.subplots(constrained_layout=True)
ax1.plot(T, Cp, color='darkred', label='Effective Cp (Parabel)')
ax1.axvspan(Ts, Tl, color='orange', alpha=0.2, label='Phasenwechselbereich (309–311 K)')
ax1.set_xlabel('Temperatur [K]')
ax1.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax1.grid(True)
ax1.legend()

if DEVELOPMENT_MODE:
    ax1.set_title('Eicosan effektives Cp (Parabel)')
else:
    fig1.savefig("eicosane_cpvst_spike_parabola.pdf", bbox_inches="tight")

# === Cp from UDF ===
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

# === UDF Cp vs T Plot ===
fig2, ax2 = plt.subplots(constrained_layout=True)
ax2.plot(T, Cp_udf, color='darkgreen', label='Specific Cp (UDF)')
ax2.axvspan(Ts, Tl, color='orange', alpha=0.2, label='Phasenwechselbereich (309–311 K)')
ax2.set_xlabel('Temperatur [K]')
ax2.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax2.grid(True)
ax2.legend()

if DEVELOPMENT_MODE:
    ax2.set_title('Eicosan Cp aus UDF')
else:
    fig2.savefig("eicosane_cpvst_udf.pdf", bbox_inches="tight")



# === Show plots if in development ===
if DEVELOPMENT_MODE:
    plt.show()