#!/usr/bin/env python3

import numpy as np
import fileinput
import sys, string, os
import ctypes
import matplotlib.pyplot as plt

from pymatgen.analysis.eos import EOS,EOSBase

f = open('ene_vs_vol_fm.dat','r')
l = f.readlines()
f.close

x = []
y = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))

vol = np.array(x)
ene = np.array(y)

#eos = EOS(eos_name='vinet')
eosfm = EOS(eos_name='murnaghan')
eos_fitfm = eosfm.fit(vol, ene)

eos_fitfm.plot()
plt.savefig("eos_fm.png")
plt.show()
plt.close()



print("EOS FM state:")

#eos_fit.plot()
print("B0 (GPa) =",eos_fitfm.b0_GPa)
print("E0 (eV) =",eos_fitfm.e0)
print("v_FM (A^3)=",eos_fitfm.v0)
print("dB0=",eos_fitfm.b1)


##################################################
f = open('ene_vs_vol_para.dat','r')
l = f.readlines()
f.close

x = []
y = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))

volp = np.array(x)
enep = np.array(y)

eosp = EOS(eos_name='murnaghan')
eos_fitp = eosp.fit(volp, enep)
eos_fitp.plot()
plt.savefig("eos_pm.png")
plt.show()
plt.close()

print("------------------")

print("EOS Paramagnetic state:")


eos_fitp.plot()
print("B0p (GPa) =",eos_fitp.b0_GPa)
print("E0p (eV) =",eos_fitp.e0)
print("v_PM (A^3) =",eos_fitp.v0)
print("dB0p=",eos_fitp.b1)

print("----------------")

print("Spontaneous volume magnetostriction:")

print("ws = (v_FM-v_PM)/v_PM = ",((eos_fitfm.v0-eos_fitp.v0)/eos_fitp.v0))

output_file = "ws.png"
plt.figure(figsize=(7, 5))
plt.plot(vol, ene, "-b", label="FM")
plt.plot(volp, enep, "--r", label="PM")
plt.xlabel("Volume (A^3)")
plt.ylabel("Energy (eV)")
plt.axvline(eos_fitfm.v0, linestyle=":", color='b', linewidth=1, label=fr"$V_0,FM={eos_fitfm.v0:.3f}$ A^3")
plt.axvline(eos_fitp.v0, linestyle="-.", color='r', linewidth=1, label=fr"$V_0,PM={eos_fitp.v0:.3f}$ A^3")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.show()

print(f"Plot saved as: {output_file}")

