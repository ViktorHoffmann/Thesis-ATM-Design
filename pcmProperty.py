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

# === Konstanten (UDF/Fluent Setup) ===
Cps = 2132.4      # J/kg·K (fest)
Cpl = 2350.05     # J/kg·K (flüssig)
Ts = 309.0        # Solidus [K]
Tl = 311.0        # Liquidus [K]
L  = 240998.9     # Latentwärme [J/kg]
Tc = (Ts + Tl) / 2

# === Temperaturbereich ===
T = np.linspace(300, 320, 2000)

# === 1) Cp aus UDF (sensible Wärmekapazität, dichte-gewichtet in Mushy-Zone) ===
Cp_udf = []
for temp in T:
    if temp < Ts:
        Cp_val = Cps
    elif temp > Tl:
        Cp_val = Cpl
    else:
        Gama = (temp - Ts) / (Tl - Ts)
        numerator   = (1 - Gama) * 910.0 * Cps + Gama * 769.0 * Cpl
        denominator = (1 - Gama) * 910.0        + Gama * 769.0
        Cp_val = numerator / denominator
    Cp_udf.append(Cp_val)

# === 2) Effektives Cp wie in Fluent (Block statt Spike) ===
# cp_eff_fluent = cp_sens + L/(Tl-Ts) in [Ts, Tl]
Cp_eff = []
block = L / (Tl - Ts)
for temp, cp_sens in zip(T, Cp_udf):
    if Ts <= temp <= Tl:
        Cp_eff.append(cp_sens + block)
    else:
        Cp_eff.append(cp_sens)

# === Plot A: Effektives Cp (Fluent-Block, KEIN Spike) ===
fig1, ax1 = plt.subplots(constrained_layout=True)
ax1.plot(T, Cp_eff, color='darkred', label='Effektive\nspezifische\nWärmekapazität')
ax1.axvspan(Ts, Tl, color='orange', alpha=0.2, label=f'Schmelzbereich')
ax1.set_xlabel('Temperatur [$\\mathrm{K}$]')
ax1.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax1.set_yscale("log")
ax1.grid(True)
ax1.legend()
# Optional: stärkere Stauchung, falls der Block dominiert
# ax1.set_yscale("symlog", linthresh=1800)

if DEVELOPMENT_MODE:
    ax1.set_title('Eicosan: effektives $c_p$ ohne Spike (Block in Mushy-Zone)')
else:
    fig1.savefig("eicosane_cpvst_total.pdf", bbox_inches="tight")

# === Plot B: Sensibles Cp (UDF) ===
fig2, ax2 = plt.subplots(constrained_layout=True)
ax2.plot(T, Cp_udf, color='darkgreen', label='Sensible\nspezifische\nWärmekapazität')
ax2.axvspan(Ts, Tl, color='orange', alpha=0.2, label=f'Schmelzbereich')
ax2.set_xlabel('Temperatur [$\\mathrm{K}$]')
ax2.set_ylabel(r'Spezifische Wärmekapazität [$\mathrm{J/(kg \cdot K)}$]')
ax2.grid(True)
ax2.legend()

if DEVELOPMENT_MODE:
    ax2.set_title('Eicosan: sensibles $c_p$ (UDF)')
    plt.show()
else:
    fig2.savefig("eicosane_cpvst_sensible.pdf", bbox_inches="tight")
