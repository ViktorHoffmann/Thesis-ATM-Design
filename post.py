import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
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
    'x', 'maxQ', 'maxQm1', 'maxQp1'
]
missing_fields = [f for f in required_fields if f not in raw.dtype.names]
if missing_fields:
    print(f"Missing fields in CSV: {missing_fields}. Exiting.")
    exit()

# === FLIGHT DATA ===
position = raw['x']  # [m]
maxQHeatflux = raw['maxQ']  # [W/m^2]
maxQm1Heatflux = raw['maxQm1']  # [W/m^2]
maxQp1Heatflux = raw['maxQp1']  # [W/m^2]

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

fig1, ax1 = plt.subplots(constrained_layout=True)

ax1.plot(x_values, maxQHeatflux, color='dimgray', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ}}$')
ax1.plot(x_values, maxQm1Heatflux, color='black', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ-1}}$')
ax1.plot(x_values, maxQp1Heatflux, color='darkgray', linestyle='-', label=r'$\dot{q}_{\mathrm{maxQ+1}}$')

ax1.set_xlabel('Position [m]')
ax1.set_ylabel(r'Spezifischer Wärmestrom [$W/m^2$]')
ax1.legend()
if DEVELOPMENT_MODE:
    ax1.set_title("spezifischer Wärmestrom Hülle")
else:
    fig1.savefig("maxQ_compare_specific_heatflux.pdf", bbox_inches="tight")

if DEVELOPMENT_MODE:
    plt.show()