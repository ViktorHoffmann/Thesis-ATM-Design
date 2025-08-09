import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# import setup
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

# import heatflux data
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

position = raw['x']  # [m]
maxQHeatflux = raw['maxQ']  # [W/m^2]
maxQm10Heatflux = raw['maxQm10']  # [W/m^2]
maxQp10Heatflux = raw['maxQp10']  # [W/m^2]
maxQp20Heatflux = raw['maxQp20']  # [W/m^2]

# import yplus data
heatfluxCompare = "yplusCompare.csv"
try:
    raw = np.genfromtxt(heatfluxCompare, delimiter=';', names=True, encoding='utf-8')
except IOError:
    print("yplusCompare not found. Exiting.")
    exit()
except ValueError:
    print("Error parsing yplusCompare. Exiting.")
    exit()

required_fields = [
    'x', 'maxQ', 'maxQm10', 'maxQp10', 'maxQp20'
]
missing_fields = [f for f in required_fields if f not in raw.dtype.names]
if missing_fields:
    print(f"Missing fields in CSV: {missing_fields}. Exiting.")
    exit()

position = raw['x']  # [m]
maxQyplus = raw['maxQ']  # [W/m^2]
maxQm10yplus = raw['maxQm10']  # [W/m^2]
maxQp10yplus = raw['maxQp10']  # [W/m^2]
maxQp20yplus = raw['maxQp20']  # [W/m^2]

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

ax1.plot(x_values, maxQm10Heatflux, color='black', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q - 10 s}}$')
ax1.plot(x_values, maxQHeatflux, color='dimgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q}}$')
ax1.plot(x_values, maxQp10Heatflux, color='darkgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q + 10 s}}$')
ax1.plot(x_values, maxQp20Heatflux, color='lightgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q + 20 s}}$')

ax1.set_xlabel('Position [m]')
ax1.set_ylabel(r'Spezifischer Wärmestrom [$W/m^2$]')
ax1.legend()
if DEVELOPMENT_MODE:
    ax1.set_title("spezifischer Wärmestrom Hülle")
else:
    fig1.savefig("maxQ_compare_heatflux.pdf", bbox_inches="tight")

# hull yplus plots
fig2, ax2 = plt.subplots(constrained_layout=True)

ax2.plot(x_values, maxQm10yplus, color='black', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q - 10 s}}$')
ax2.plot(x_values, maxQyplus, color='dimgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q}}$')
ax2.plot(x_values, maxQp10yplus, color='darkgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q + 10 s}}$')
ax2.plot(x_values, maxQp20yplus, color='lightgrey', linestyle='-', label=r'$\dot{q}_{\mathrm{max Q + 20 s}}$')

ax2.set_xlabel('Position [m]')
ax2.set_ylabel('y+')
ax2.legend()
if DEVELOPMENT_MODE:
    ax2.set_title("y+ entlang der Hülle")
else:
    fig2.savefig("maxQ_compare_yplus.pdf", bbox_inches="tight")

if DEVELOPMENT_MODE:
    plt.show()