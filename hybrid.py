# This Program calculates the transient equation of a radiator-pcm hybrid for the required pcm capacity
#
# setup.json: devmode, specificheatflux_environment, flight_duration ---> hybrid.py: pcm_capacity ---> result.json

import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib as mpl



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
Qdot_avionics = data["Qdot_avionics"]
trajectory_data_path = data["trajectory_data_path"]
a = data["absorptance"]

try:
    with open("result.json", "r") as f:
        data = json.load(f)
except FileNotFoundError:
    print("result.json not found. Exiting.")
    exit()
except json.JSONDecodeError:
    print("result.json is empty of invalid. Exiting.")
    exit()

phi_radiation_target = data["phi_radiation_target"]
radiator_area_target = data["radiator_area_target"]

# Load CSV data
try:
    raw = np.genfromtxt(trajectory_data_path, delimiter=';', names=True, encoding='utf-8')
except IOError:
    print("trajectory_data_path not found. Exiting.")
    exit()
except ValueError:
    print("Error parsing trajectory_data_path. Exiting.")
    exit()

required_fields = ['time', 'altitude', 'velocity', 'acceleration', 'air_temperature', 'air_pressure', 'air_density', 'reynoldsnumber']
missing_fields = [f for f in required_fields if f not in raw.dtype.names]
if missing_fields:
    print(f"Missing fields in CSV: {missing_fields}. Exiting.")
    exit()

time = raw['time'] # [s]
altitude = raw['altitude'] # [m]
velocity = raw['velocity'] # [m/s]
acceleration = raw['acceleration'] # [m/s²]
air_temperature = raw['air_temperature'] + 273.15 # convert °C to Kelvin
air_pressure = raw['air_pressure'] # [Pa]
air_density = raw['air_density'] # bad data
reynoldsnumber = raw['reynoldsnumber'] # unused because not local number

# Set global matplotlib style for latex integration
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

# Mock data setup for rocket skin heatflux [W]
# def skewed_gaussian(x, mu, sig_left, sig_right, peak=1000, noise_level=0.02, additive_noise=5):
#     y = np.zeros_like(x)
#     left = x < mu
#     right = x >= mu
#     y[left] = np.exp(-((x[left] - mu) ** 2) / (2 * sig_left ** 2))
#     y[right] = np.exp(-((x[right] - mu) ** 2) / (2 * sig_right ** 2))
#     y = y / np.max(y) * peak
#     if noise_level > 0:
#         y *= np.random.normal(1.0, noise_level, size=x.shape)
#     if additive_noise > 0:
#         y += np.random.normal(0.0, additive_noise, size=x.shape)
#     return np.clip(y, 0, None)

# Qdot_env = skewed_gaussian(
#     x_values, mu=200, sig_left=40, sig_right=5,
#     peak=100, noise_level=0.002, additive_noise=2
# )



# === Equations ===
eta_0 = 18.27e-6    # [Pa*s]
T_0 = 291.15        # [K] reference temperature for sutherlands law
C = 120             # [K]
kappa = 1.4
R = 287             # [J/(kg*K)]
c_p = 1005          # [J/(kg*K)]
x = 1               # [m] for local evaluation
T_w = 273.15+38     # [K] assumed constant at pcm melting point (icosane: 38°C)

# average temperature
def T_m(T1,T2):
    return (T1+T2)/2

# Sutherlands law
def eta(T):
    return eta_0 * ((T_0 + C)/(T + C)) * (T/T_0)**(3/2) # [Pa*s]

# Mach number
def Ma(V, T):
    return V/((kappa * R * T)**(1/2))

# thermal conductivity of air
def lam(T):
    # T in Kelvin
    A = 2.64638e-3
    B = 7.326e-5
    C = -1.746e-8
    return A + B * T + C * T**2  # [W/m·K]

# air density with ideal gas law
def rho(p, T):
    return p/(R * T) # [kg/m³]

# Prandtel number
def Pr(T):
    return (c_p * eta(T))/(lam(T))

# Reynoldsnumber
def Re(V, p, T, x):
    return (V * rho(p, T) * x)/eta(T)

# recovery factor with approximation
def r(T):
    return Pr(T)**(1/3)

# recovery temperature
def T_r(V, T):
    return T * (1 + r(T) * (kappa + 1)/2 * Ma(V, T)) # [K]

# specific heatflux laminar boundarylayer
def qdot_air(p, T, V, x, T_w):
    Re_x = Re(V, p, T, x)
    Pr_x = Pr(T)
    Nu_x = 0.332 * Re_x**(1/2) * Pr_x**(1/3)
    alpha = Nu_x * lam(T) / x
    Tr = T_r(V, T)
    return alpha * (Tr - T_w)  # [W/m²]

Qdot_env = np.array([
    qdot_air(p, T_m(T_w, T), V, x, T_w)
    for p, T, V in zip(air_pressure, air_temperature, velocity)
]) * radiator_area_target + (500*a*radiator_area_target) # environment includes solar flux

Qdot_in = Qdot_env + Qdot_avionics  # [W]

# testing
Re_plot = np.array([
    Re(V, p, T, x)
    for V, p, T in zip(velocity, air_pressure, air_temperature)
    ])
# testing
Pr_plot = np.array([
    Pr(T)
    for T in air_temperature
])


# === PLOTTING ===
x_values = time

# Plot the system heatflows
fig1, ax1 = plt.subplots(constrained_layout=True)

ax1.plot(x_values, Qdot_env, color='blue', linewidth=1, linestyle=(0, (1, 5)),
        label=r'$\dot{Q}_{\mathrm{Umwelt}}$')
ax1.plot(x_values, [Qdot_avionics]*len(x_values), color='orange', linewidth=1, linestyle='--',
        label=r'$\dot{Q}_{\mathrm{Avionik}}$')
ax1.plot(x_values, [phi_radiation_target]*len(x_values), color='red', linewidth=1, linestyle='--',
        label=r'$\dot{Q}_{\mathrm{Radiator}}$')
ax1.plot(x_values, Qdot_in, color='blue', linewidth=1, linestyle='-',
        label=r'$\dot{Q}_{\mathrm{Rein}}$')

ax1.fill_between(
    x_values, phi_radiation_target, Qdot_in,
    where=(Qdot_in > phi_radiation_target),
    facecolor='lightgrey',
    hatch='//',
    edgecolor='black',
    alpha=0.5,
    label='PCM\nSchmelz-\nbereich'
)

# Inset axes (zoom region)
x1, x2, y1, y2 = 0, 70, 0, 190 # adjust as needed
axins = ax1.inset_axes([0.2, 0.47, 0.5, 0.5])
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Plot same data in inset
axins.plot(x_values, Qdot_env, color='blue', linewidth=1, linestyle=(0, (1, 5)))
axins.plot(x_values, [Qdot_avionics] * len(x_values), color='orange', linewidth=1, linestyle='--')
axins.plot(x_values, [phi_radiation_target] * len(x_values), color='red', linewidth=1, linestyle='--')
axins.plot(x_values, Qdot_in, color='blue', linewidth=1, linestyle='-')
axins.fill_between(
    x_values, phi_radiation_target, Qdot_in,
    where=(Qdot_in > phi_radiation_target),
    facecolor='lightgrey',
    hatch='//',
    edgecolor='black',
    alpha=0.5,
)

# Labeling
ax1.set_xlabel('Zeit [s]')
ax1.set_ylabel(r'Wärmestrom [W]')
ax1.legend()
if DEVELOPMENT_MODE:
    ax1.set_title("PCM Wärmestrom")
else:
    # Save LaTeX figure
    fig1.savefig("pcm_radiator_hybrid_heatflux_during_flight.pdf", bbox_inches="tight")

# Plot the Reynoldsnumber and Prandtlnumber
fig2, ax2 = plt.subplots(constrained_layout=True)

# Plot Reynolds number
ax2.plot(x_values, Re_plot, color='black', label='Reynolds-Zahl')
ax2.set_xlabel('Zeit [s]')
ax2.set_ylabel('Reynolds-Zahl')
ax2.tick_params(axis='y')
# ax2.axhline(5e5, color='green', linestyle=':', linewidth=1, label=r'$Re_\mathrm{min} = 5 \times 10^5$')
# ax2.axhline(1e7, color='green', linestyle=':', linewidth=1, label=r'$Re_\mathrm{max} = 1 \times 10^7$')

# Create a second y-axis for Prandtl number
ax2b = ax2.twinx()
ax2b.plot(x_values, Pr_plot, linestyle='--', color='black', label='Prandtl-Zahl')
ax2b.set_ylabel('Prandtl-Zahl')
ax2b.tick_params(axis='y')

# Combine legends from both y-axes
lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax2b.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper right')
ax2.set_xlim(0,200)
if DEVELOPMENT_MODE:
    ax2.set_title("Reynolds- und Prandtl-Zahl während kritischer Phase im Flug")
else:
    # Save LaTex figure
    fig2.savefig("re_pr_during_flight.pdf", bbox_inches="tight")



# === Target Result ===
# Area calculation for PCM capacity requirement
mask = Qdot_in > phi_radiation_target
x_limits = x_values[mask]
y_diff = Qdot_in[mask] - phi_radiation_target
pcm_capacity_target = np.trapezoid(y_diff, x_limits)

# Append new data into result.json
with open("result.json", "r") as f:
    data = json.load(f)
data["pcm_capacity_target"] = pcm_capacity_target
with open("result.json", "w") as f:
    json.dump(data, f, indent=4)

if DEVELOPMENT_MODE:
    plt.show()