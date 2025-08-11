# This Program calculates the transient equation of a radiator-pcm hybrid for the required pcm capacity
# it additionally plots Re, Pr and dynamic pressure for given flight data

import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.optimize import curve_fit

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
acceleration = raw['acceleration']
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
x = 1.07          # [m] radiator centerpoint (0.06 m from hull top)
T_w = 273.15 + target_temperature # [K] PCM melting point (icosane: 38°C)

# === THERMOPHYSICAL FUNCTIONS ===
def T_m(T1, T2): return (T1 + T2) / 2                                   # average temperature
def eta(T): return eta_0 * ((T_0 + C) / (T + C)) * (T / T_0) ** (3/2)   # dynamic viscosity with surherlands formula
def lam(T): return 2.64638e-3 + 7.326e-5 * T - 1.746e-8 * T**2          # thermal conductivity with polynomial fit
def rho(p, T): return p / (R * T)                                       # air density
def Pr(T): return (c_p * eta(T)) / lam(T)                               # prandtl number
def Ma(V, T): return V / np.sqrt(kappa * R * T)                         # mach number
def Re(V, p, T, x): return V * rho(p, T) * x / eta(T)                   # reynolds number
def r(T): return Pr(T) ** (1/3)                                         # recovery factor
def T_r(V, T): return T * (1 + r(T) * (kappa + 1) / 2 * Ma(V, T))       # recovery temperature
def qdot_air(p, T, V, x, T_w):                                          # nusselt relation for wall heatflux
    Re_x = Re(V, p, T, x)
    Pr_x = Pr(T)
    Nu_x = 0.0296 * Re_x**0.8 * Pr_x**(1/3) # turbulent
    alpha = Nu_x * lam(T) / x
    return alpha * (T_r(V, T) - T_w)
def pdyn(V, T, p): return 0.5 * rho(p, T) * V**2                        # dynamic pressure

# === HEATFLUX CALCULATION ===
Qdot_env = np.array([
    qdot_air(p, T_m(T_w, T), V, x, T_w)
    for p, T, V in zip(air_pressure, air_temperature, velocity)
]) * hybrid_radiator_area + (solar_flux/2 * a * hybrid_radiator_area)  # add solar flux

Qdot_in = Qdot_env + avionics_power  # [W]

# === FLUID DYNAMICS PLOTS ===
Re_plot = np.array([Re(V, p, T, x) for V, p, T in zip(velocity, air_pressure, air_temperature)])
Pr_plot = np.array([Pr(T) for T in air_temperature])
pdyn_plot = np.array([pdyn(V, T, p) for V, p, T in zip(velocity, air_pressure, air_temperature)])

# === PLOTTING HEATFLUX ===
x_values = time

gauss_x = np.array([18.691, 28.691, 38.691, 48.7])  # seconds
gauss_y = np.array([
    34567.9 * hybrid_radiator_area,
    67706.8 * hybrid_radiator_area,
    72349.2 * hybrid_radiator_area,
    21479.6 * hybrid_radiator_area
])  # W

# Gaussian function
def gaussian(x, a, b, c, d):
    return a * np.exp(-((x - b)**2) / (2 * c**2)) + d

# Initial guess: [amplitude, mean, stddev, offset]
guess = [np.max(gauss_y), gauss_x[np.argmax(gauss_y)], 10, np.min(gauss_y)]

# Fit Gaussian to 4 points
popt, _ = curve_fit(gaussian, gauss_x, gauss_y, p0=guess)

# Generate smooth curve
gauss_fitted_raw = gaussian(x_values, *popt)
gauss_fitted = np.clip(gauss_fitted_raw, 0, None)  # Force values >= 0

# Main heat flow curves
fig1a, ax1a = plt.subplots(constrained_layout=True)
ax1a.plot(x_values, Qdot_env, color='blue', linestyle=(0, (1, 5)), label=r'$\dot{Q}_{\mathrm{Umwelt}}$')
ax1a.plot(x_values, [avionics_power]*len(x_values), color='orange', linestyle='--', label=r'$\dot{Q}_{\mathrm{Avionik}}$')
ax1a.plot(x_values, [hybrid_radiator_power]*len(x_values), color='red', linestyle='--', label=r'$\dot{Q}_{\mathrm{Radiator}}$')
ax1a.plot(x_values, Qdot_in, color='blue', linestyle='-', label=r'$\dot{Q}_{\mathrm{Rein}}$')
ax1a.fill_between(
    x_values, hybrid_radiator_power, Qdot_in,
    where=(Qdot_in > hybrid_radiator_power),
    facecolor='lightgrey',
    hatch='//',
    edgecolor='black',
    alpha=0.5,
    label='PCM\nSchmelz-\nbereich'
)
ax1a.set_xlabel('Zeit [s]')
ax1a.set_ylabel(r'W\"armestrom [W]')
ax1a.legend()
x1, x2, y1, y2 = 0, 125, 0, 100
axins = ax1a.inset_axes([0.2, 0.4, 0.5, 0.5])
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
mark_inset(ax1a, axins, loc1=2, loc2=4, fc="none", ec="0.5")
if DEVELOPMENT_MODE:
    ax1a.set_title("PCM W\"armestrom (ohne Simulation)")
else:
    fig1a.savefig("pcm_radiator_hybrid_heatflux_nosim.pdf", bbox_inches="tight")


fig1, ax1 = plt.subplots(constrained_layout=True)

# Main heat flow curves with simulation data
ax1.plot(x_values, Qdot_env, color='blue', linestyle=(0, (1, 5)), label=r'$\dot{Q}_{\mathrm{Umwelt}}$')
ax1.plot(x_values, [avionics_power] * len(x_values), color='orange', linestyle='--', label=r'$\dot{Q}_{\mathrm{Avionik}}$')
ax1.plot(x_values, [hybrid_radiator_power] * len(x_values), color='red', linestyle='--', label=r'$\dot{Q}_{\mathrm{Radiator}}$')
ax1.plot(x_values, Qdot_in, color='blue', linestyle='-', label=r'$\dot{Q}_{\mathrm{Rein}}$')
ax1.plot(28.691, (67706.8*hybrid_radiator_area), color='red',marker='o', label=r'$\dot{Q}_{\mathrm{Sim}}$')
ax1.plot(18.691, (34567.9*hybrid_radiator_area), color='red',marker='o')
ax1.plot(38.691, (72349.2*hybrid_radiator_area), color='red',marker='o')
ax1.plot(48.7, (21479.6*hybrid_radiator_area), color='red',marker='o')
ax1.plot(x_values, gauss_fitted, color='green', linestyle='-.', label=r'$\dot{Q}_{\mathrm{Sim-Fit}}$')

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
x1, x2, y1, y2 = 0, 60, 0, 7500
axins = ax1.inset_axes([0.2, 0.3, 0.4, 0.6])
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.plot(x_values, Qdot_env, color='blue', linestyle=(0, (1, 5)))
axins.plot(x_values, [avionics_power] * len(x_values), color='orange', linestyle='--')
axins.plot(x_values, [hybrid_radiator_power] * len(x_values), color='red', linestyle='--')
axins.plot(x_values, Qdot_in, color='blue', linestyle='-')
axins.plot(28.691, (67706.8*hybrid_radiator_area), color='red',marker='o', label=r'$\dot{Q}_{\mathrm{Sim}}$')
axins.plot(18.691, (34567.9*hybrid_radiator_area), color='red',marker='o')
axins.plot(38.691, (72349.2*hybrid_radiator_area), color='red',marker='o')
axins.plot(48.7, (21479.6*hybrid_radiator_area), color='red',marker='o')
axins.plot(x_values, gauss_fitted, color='green', linestyle='-.', label=r'$\dot{Q}_{\mathrm{Sim-Fit}}$')

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
    print(f"Fitted Gaussian parameters:\na={popt[0]:.2f}, b={popt[1]:.2f}, c={popt[2]:.2f}, d={popt[3]:.2f}")
else:
    fig1.savefig("pcm_radiator_hybrid_heatflux_with_sim.pdf", bbox_inches="tight")

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
ax2.set_xlim(0, 150)

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

# === PLOTTING DYNAMIC PRESSURE ===
fig3, ax3 = plt.subplots(constrained_layout=True)
ax3.plot(x_values, pdyn_plot, color='black', label='Dynamischer Druck')
ax3.set_xlabel('Zeit [s]')
ax3.set_ylabel(r'Dynamischer Druck [Pa]')
ax3.tick_params(axis='y')
ax3.set_xlim(0, 150)

if DEVELOPMENT_MODE:
    ax3.set_title("Dynamischer Druck während kritischer Phase im Flug")
    dynp_plot_clean = np.ma.masked_invalid(pdyn_plot)
    # Get max value and index
    max_index = np.ma.argmax(dynp_plot_clean)
    max_value = pdyn_plot[max_index]
    max_x = x_values[max_index]
    print(f"Maximum Dynamic Pressure: {max_value} at x = {max_x}")
else:
    fig3.savefig("dynp_during_flight.pdf", bbox_inches="tight")

# === PLOT: ALTITUDE OVER TIME ===
fig4, ax4 = plt.subplots(constrained_layout=True)
ax4.plot(time, raw['altitude'], color='navy')
ax4.set_xlabel("Zeit [s]")
ax4.set_ylabel("Höhe [m]")
ax4.set_title("Flug Höhe" if DEVELOPMENT_MODE else "")
if not DEVELOPMENT_MODE:
    fig4.savefig("altitude_over_time.pdf", bbox_inches="tight")

# === PLOT: VELOCITY OVER TIME ===
fig5, ax5 = plt.subplots(constrained_layout=True)
ax5.plot(time, velocity, color='darkgreen')
ax5.set_xlabel("Zeit [s]")
ax5.set_ylabel("Geschwindigkeit [m/s]")
ax5.set_title("Flug Geschwindigkeit" if DEVELOPMENT_MODE else "")
if not DEVELOPMENT_MODE:
    fig5.savefig("velocity_over_time.pdf", bbox_inches="tight")

# === PLOT: PRESSURE OVER TIME ===
fig6, ax6 = plt.subplots(constrained_layout=True)
ax6.plot(time, air_pressure, color='purple')
ax6.set_xlabel("Zeit [s]")
ax6.set_ylabel("Luftdruck [Pa]")
ax6.set_title("Flug Luftdruck" if DEVELOPMENT_MODE else "")
if not DEVELOPMENT_MODE:
    fig6.savefig("pressure_over_time.pdf", bbox_inches="tight")

# === PLOT: TEMPERATURE OVER TIME ===
fig7, ax7 = plt.subplots(constrained_layout=True)
ax7.plot(time, air_temperature, color='firebrick')
ax7.set_xlabel("Zeit [s]")
ax7.set_ylabel("Temperatur [K]")
ax7.set_title("Flug Lufttemperatur" if DEVELOPMENT_MODE else "")
if not DEVELOPMENT_MODE:
    fig7.savefig("temperature_over_time.pdf", bbox_inches="tight")

# === PLOT: ACCELERATION OVER TIME ===
fig8, ax8 = plt.subplots(constrained_layout=True)
ax8.plot(time, acceleration, color='darkkhaki')
ax8.set_xlabel("Zeit [s]")
ax8.set_ylabel(r"Beschleunigung [$\mathrm{m/s^2}$]")
ax8.set_title("Flug Beschleunigung" if DEVELOPMENT_MODE else "")
if not DEVELOPMENT_MODE:
    fig8.savefig("acceleration_over_time.pdf", bbox_inches="tight")

# === APPROXIMATED ACCELERATION FROM UDF ===
acc_udf = np.zeros_like(time)
for i, t in enumerate(time):
    if t < 20:
        acc_udf[i] = 34.81
    elif t < 50:
        acc_udf[i] = 109.81
    elif t < 150:
        acc_udf[i] = 19.62
    else:
        acc_udf[i] = 9.81

fig9, ax9 = plt.subplots(constrained_layout=True)
ax9.plot(time, acc_udf, color='darkred', label='Approximierte Beschleunigung')
ax9.set_xlabel("Zeit [s]")
ax9.set_ylabel(r"Beschleunigung [$\mathrm{m/s^2}$]")
ax9.set_title("Approximierte Beschleunigung" if DEVELOPMENT_MODE else "")
ax9.legend()

if not DEVELOPMENT_MODE:
    fig9.savefig("approximate_acceleration_over_time.pdf", bbox_inches="tight")


# calculate simulation fit capacity result
mask = gauss_fitted > hybrid_radiator_power
x_sim = x_values[mask]
y_sim = gauss_fitted[mask] - hybrid_radiator_power
hybrid_capacity_sim = float(np.trapezoid(y_sim, x_sim)) if x_sim.size else 0.0  # [J]

# calculate nusselt capacity result
mask_nu = Qdot_in > hybrid_radiator_power
x_nu = x_values[mask_nu]
y_nu = Qdot_in[mask_nu] - hybrid_radiator_power
hybrid_capacity_nu = float(np.trapezoid(y_nu, x_nu)) if x_nu.size else 0.0      # [J]

# simple PCM energy budget (constant duration)
pcm_capacity = float(avionics_power * 1200)  # [J]

# Append result
result_data.update({
    "hybrid_capacity_sim": hybrid_capacity_sim,
    "hybrid_capacity_nu": hybrid_capacity_nu,
    "pcm_capacity": pcm_capacity,
})

with open("result.json", "w") as f:
    json.dump(result_data, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()