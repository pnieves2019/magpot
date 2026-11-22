#!/usr/bin/env python3

import numpy as np
import fileinput
import sys, string, os
import ctypes
import matplotlib.pyplot as plt

from pymatgen.analysis.eos import EOS,EOSBase




f = open('ene_vs_vol.dat','r')
l = f.readlines()
f.close

x = []
y = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))

vol = np.array(x)
ene = np.array(y)

eos = EOS(eos_name='vinet')
#eos = EOS(eos_name='murnaghan')
eos_fit = eos.fit(vol, ene)

eos_fit.plot()

ax = plt.gca()

ax.set_xlabel(r"Volume ($\mathrm{\AA}^3$)", fontsize=13)
ax.set_ylabel("Energy (eV)", fontsize=13)

plt.tight_layout()
plt.savefig("eos.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

print("EOS:")

#eos_fit.plot()
print("B0(GPa)=",eos_fit.b0_GPa)
print("E0(eV)=",eos_fit.e0)
print("v0(A^3)=",eos_fit.v0)
print("dB0=",eos_fit.b1)

