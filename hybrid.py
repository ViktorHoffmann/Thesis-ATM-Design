# This Program calculates the transient equation of a radiator-pcm hybrid for the required pcm capacity
#
# setup.json: devmode, specificheatflux_environment, flight_duration ---> hybrid.py: pcm_capacity ---> result.json

import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# === CONFIGURATION ===

# Load data from setup.json and result.json
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
avionics_power = data["avionics_power"]
trajectory_data_path = data["trajectory_data_path"]
a = data["absorptance"]
solar_flux = data["solar_flux"]
target_temperature = data["target_temperature"]

try:
    with open("result.json", "r") as f:
        result_data = json.load(f)
except FileNotFoundError:
    print("result.json not found. Exiting.")
    exit()
except json.JSONDecodeError:
    print("result.json is empty or invalid. Exiting.")
    exit()

hybrid_radiator_area = result_data["hybrid_radiator_area"]
hybrid_radiator_power = result_data["hybrid_radiator_power"]

# Load trajectory CSV data
try:
    raw = np.genfromtxt(trajectory_data_path, delimiter=';', names=True, encoding='utf-8')
except IOError:
    print("trajectory_data_path not found. Exiting.")
    exit()
except ValueError:
    print("Error parsing trajectory_data_path. Exiting.")
    exit()

required_fields = [
    'time', 'altitude', 'velocity', 'acceleration',
    'air_temperature', 'air_pressure', 'air_density', 'reynoldsnumber'
]
missing_fields = [f for f in required_fields if f not in raw.dtype.names]
if missing_fields:
    print(f"Missing fields in CSV: {missing_fields}. Exiting.")
    exit()

# === FLIGHT DATA ===
time = raw['time']  # [s]
velocity = raw['velocity']  # [m/s]
air_temperature = raw['air_temperature'] + 273.15  # [K]
air_pressure = raw['air_pressure'] * 100 # [Pa]

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

# === CONSTANTS ===
eta_0 = 18.27e-6  # [Pa*s]
T_0 = 291.15      # [K]
C = 120           # [K]
kappa = 1.4
R = 287           # [J/(kg·K)]
c_p = 1005        # [J/(kg·K)]
x = 1.25             # [m]
T_w = 273.15 + target_temperature # [K] PCM melting point (icosane: 38°C)

# === THERMOPHYSICAL FUNCTIONS ===
def T_m(T1, T2): return (T1 + T2) / 2
def eta(T): return eta_0 * ((T_0 + C) / (T + C)) * (T / T_0) ** (3/2)
def lam(T): return 2.64638e-3 + 7.326e-5 * T - 1.746e-8 * T**2
def rho(p, T): return p / (R * T)
def Pr(T): return (c_p * eta(T)) / lam(T)
def Ma(V, T): return V / np.sqrt(kappa * R * T)
def Re(V, p, T, x): return V * rho(p, T) * x / eta(T)
def r(T): return Pr(T) ** (1/3)
def T_r(V, T): return T * (1 + r(T) * (kappa + 1) / 2 * Ma(V, T))
def qdot_air(p, T, V, x, T_w):
    Re_x = Re(V, p, T, x)
    Pr_x = Pr(T)
    Nu_x = 0.0296 * Re_x**0.8 * Pr_x**(1/3)
    alpha = Nu_x * lam(T) / x
    return alpha * (T_r(V, T) - T_w)

# === HEATFLUX CALCULATION ===
Qdot_env = np.array([
    qdot_air(p, T_m(T_w, T), V, x, T_w)
    for p, T, V in zip(air_pressure, air_temperature, velocity)
]) * hybrid_radiator_area + (solar_flux/2 * a * hybrid_radiator_area)  # add solar flux

Qdot_in = Qdot_env + avionics_power  # [W]

# === FLUID DYNAMICS PLOTS ===
Re_plot = np.array([Re(V, p, T, x) for V, p, T in zip(velocity, air_pressure, air_temperature)])
Pr_plot = np.array([Pr(T) for T in air_temperature])

# === PLOTTING HEATFLUX ===
x_values = time

fig1, ax1 = plt.subplots(constrained_layout=True)

# Main heat flow curves
ax1.plot(x_values, Qdot_env, color='blue', linestyle=(0, (1, 5)), label=r'$\dot{Q}_{\mathrm{Umwelt}}$')
ax1.plot(x_values, [avionics_power] * len(x_values), color='orange', linestyle='--', label=r'$\dot{Q}_{\mathrm{Avionik}}$')
ax1.plot(x_values, [hybrid_radiator_power] * len(x_values), color='red', linestyle='--', label=r'$\dot{Q}_{\mathrm{Radiator}}$')
ax1.plot(x_values, Qdot_in, color='blue', linestyle='-', label=r'$\dot{Q}_{\mathrm{Rein}}$')

# PCM melting range
ax1.fill_between(
    x_values, hybrid_radiator_power, Qdot_in,
    where=(Qdot_in > hybrid_radiator_power),
    facecolor='lightgrey',
    hatch='//',
    edgecolor='black',
    alpha=0.5,
    label='PCM\nSchmelz-\nbereich'
)

# Inset zoom
x1, x2, y1, y2 = 0, 100, 0, 150
axins = ax1.inset_axes([0.2, 0.4, 0.5, 0.5])
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.plot(x_values, Qdot_env, color='blue', linestyle=(0, (1, 5)))
axins.plot(x_values, [avionics_power] * len(x_values), color='orange', linestyle='--')
axins.plot(x_values, [hybrid_radiator_power] * len(x_values), color='red', linestyle='--')
axins.plot(x_values, Qdot_in, color='blue', linestyle='-')
axins.fill_between(
    x_values, hybrid_radiator_power, Qdot_in,
    where=(Qdot_in > hybrid_radiator_power),
    facecolor='lightgrey',
    hatch='//',
    edgecolor='black',
    alpha=0.5
)
mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# Labels and output
ax1.set_xlabel('Zeit [s]')
ax1.set_ylabel(r'Wärmestrom [W]')
ax1.legend()
if DEVELOPMENT_MODE:
    ax1.set_title("PCM Wärmestrom")
else:
    fig1.savefig("pcm_radiator_hybrid_heatflux_during_flight.pdf", bbox_inches="tight")

# === PLOTTING FLUID NUMBERS ===
fig2, ax2 = plt.subplots(constrained_layout=True)
ax2.plot(x_values, Re_plot, color='black', label='Reynolds-Zahl')
ax2.set_xlabel('Zeit [s]')
ax2.set_ylabel('Reynolds-Zahl')
ax2.tick_params(axis='y')

# Prandtl number on secondary axis
ax2b = ax2.twinx()
ax2b.plot(x_values, Pr_plot, linestyle='--', color='black', label='Prandtl-Zahl')
ax2b.set_ylabel('Prandtl-Zahl')
ax2b.tick_params(axis='y')

# Combined legend
lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax2b.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper right')
ax2.set_xlim(0, 200)

if DEVELOPMENT_MODE:
    ax2.set_title("Reynolds- und Prandtl-Zahl während kritischer Phase im Flug")
    Re_plot_clean = np.ma.masked_invalid(Re_plot)
    # Get max value and index
    max_index = np.ma.argmax(Re_plot_clean)
    max_value = Re_plot[max_index]
    max_x = x_values[max_index]
    print(f"Maximum Reynolds number: {max_value} at x = {max_x}")
else:
    fig2.savefig("re_pr_during_flight.pdf", bbox_inches="tight")

# === CALCULATE PCM CAPACITY REQUIREMENT ===
mask = Qdot_in > hybrid_radiator_power
x_limits = x_values[mask]
y_diff_hybrid = Qdot_in[mask] - hybrid_radiator_power
y_diff_pcm = avionics_power
hybrid_capacity_target = np.trapezoid(y_diff_hybrid, x_limits)  # hybrid Energy integral [J]
pcm_capacity_target = avionics_power * 1400 # normal pcm energy requirement

# Append result
result_data["hybrid_capacity_target"] = hybrid_capacity_target
result_data["pcm_capacity_target"] = pcm_capacity_target
with open("result.json", "w") as f:
    json.dump(result_data, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()