# Das Programm ist nur semi-relevant, da es pcm masse funktion plottet was bei hybrid lösung egal ist

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# === Constants ===
rho_alu = 2700      # kg/m^3
cp_alu = 896        # J/kgK
rho_pcm = 788       # kg/m^3
cp_pcm = 2950       # J/kgK
dT = 30             # K
dH = 266000         # J/kg
F = 0.1             # void fraction
h = 0.001           # wall thickness [m]
Q_target = 18000    # J, desired total heat

# === Matplotlib style ===
mpl.rcParams.update({
    "figure.figsize": (6.5, 4),
    "font.size": 11.0,
    "font.family": "serif",
    "axes.titlesize": "medium",
    "figure.titlesize": "medium"
})

# === Prepare grid ===
L_vals = np.linspace(2*h, 0.1, 150)
H_vals = np.linspace(2*h, 0.1, 150)
L_grid, H_grid = np.meshgrid(L_vals, H_vals)

# === Calculate Q_ges and m_ges ===
outer_vol = L_grid**2 * H_grid
inner_vol = (L_grid - 2*h)**2 * (H_grid - 2*h)
alu_shell_vol = outer_vol - inner_vol

Q_ges = (
    rho_alu * alu_shell_vol * cp_alu * dT +
    F * rho_alu * inner_vol * cp_alu * dT +
    (1 - F) * rho_pcm * inner_vol * (cp_pcm * dT + dH)
)

m_ges = (
    rho_alu * alu_shell_vol +
    (F * rho_alu + (1 - F) * rho_pcm) * inner_vol
)

# === Find contour for Q_target ===
fig, ax1 = plt.subplots(figsize=(6.5, 4))
contours = ax1.contour(L_grid, H_grid, Q_ges, levels=[Q_target], colors='red')

# === Prepare secondary axis for mass ===
ax2 = ax1.twinx()
ax2.set_ylabel("Masse $m$ [kg]")

# Extract contour points and compute mass
if len(contours.allsegs[0]) > 0:
    for segment in contours.allsegs[0]:
        Ls = segment[:, 0]
        Hs = segment[:, 1]
        ms = []

        for Lp, Hp in zip(Ls, Hs):
            outer = Lp**2 * Hp
            inner = (Lp - 2*h)**2 * (Hp - 2*h)
            alu = outer - inner
            m = rho_alu * alu + (F * rho_alu + (1 - F) * rho_pcm) * inner
            ms.append(m)

        # Sort by L for clean plotting
        sorted_indices = np.argsort(Ls)
        Ls_sorted = Ls[sorted_indices]
        Hs_sorted = Hs[sorted_indices]
        ms_sorted = np.array(ms)[sorted_indices]

        ax1.plot(Ls_sorted, Hs_sorted, label="Q = Target", color="red")
        ax2.plot(Ls_sorted, ms_sorted, label="Masse", color="blue", linestyle="--")

else:
    print(f"No contour found for Q_target = {Q_target} J. Try changing the range or resolution.")

# === Labels and layout ===
ax1.set_xlabel("Länge $L$ [m]")
ax1.set_ylabel("Höhe $H$ [m]")
ax1.set_title(f"Systemmasse entlang $Q_\\mathrm{{ges}}$ = {Q_target:.0f} J")
ax1.grid(True)

# Optional: combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

plt.tight_layout()
plt.savefig("mass_vs_geometry_fixed_heat.pdf", bbox_inches="tight")
plt.show()