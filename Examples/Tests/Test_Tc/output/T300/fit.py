#!/usr/bin/env python3
"""Plot Mavg versus temperature and estimate Tc from thermal.norm.dat."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Input and output files
input_file = "thermal.norm.dat"
output_file = "magnetization_vs_temperature.png"


def critical_law(T, Tc, beta):
    """Critical behaviour of the order parameter below Tc."""
    return (1.0 - T / Tc)**beta


# Read the first two numerical columns (Temp. and Mavg).
# Header lines, blank lines and comment lines are ignored.
data = np.loadtxt("thermal.norm.dat")
T = data[:,0]
M = data[:,1]

beta_guess = 0.36
Tc_guess = 1000.0

popt, pcov, infodict, errmsg, ier = curve_fit(critical_law,T,M,p0=[Tc_guess,beta_guess],full_output=True)

Tc, beta = popt
errors = np.sqrt(np.diag(pcov))
Tc_error = errors[0]

E = (infodict['fvec']**2).sum()
r2 = 1-E/((M-M.mean())**2).sum()


print("Critical-law fit: m(T) = (1 - T/Tc)^beta")
print("------------------------------------------------")
print(f"Tc   = {Tc:.3f} +/- {Tc_error:.3f} K")
print(f"beta = {beta:.4f} +/- {errors[1]:.4f}")
print("R2=",r2)

# Plot simulation data and fitted critical law.
T_curve = np.linspace(np.min(T), Tc, 400)
M_curve = critical_law(T_curve, Tc, beta)

plt.figure(figsize=(7, 5))
plt.plot(T, M, "o", label="Simulation")
plt.plot(T_curve, M_curve, "-", label=fr"Fit: $T_c={Tc:.2f}$ K, $\beta={beta:.3f}$")
plt.axvline(Tc, linestyle="--", linewidth=1, label=fr"$T_c={Tc:.2f}$ K")
plt.xlabel("Temperature (K)")
plt.ylabel("M/Ms")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.show()

print(f"Plot saved as: {output_file}")

