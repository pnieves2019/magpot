#!/usr/bin/env python3

import numpy as np
import fileinput
import sys, string, os
import ctypes
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




f = open('ene_vs_vol_fm_110.dat','r')
l = f.readlines()
f.close

x = []
y = []
z = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))
    z.append(float(i.split()[2]))

vol = np.array(x)
press = np.array(y) # bar
ene110 = np.array(z)

f = open('ene_vs_vol_fm_100.dat','r')
l2 = f.readlines()
f.close

z2 = []

for i in l2:
    z2.append(float(i.split()[2]))

ene100 = np.array(z2)


def f(x, a, b):
  return a*x+b


vol0=vol*10**-30 #m^3
dE=(ene110-ene100)*(1.60218e-19) # Joule 
K1=(4.0*dE)/vol0 #J/m^3
press0=press*1.0e5 # Pa

xdib = np.linspace(press0[0], press0[-1], 20)

p, pcov, infodict, errmsg, ier = curve_fit(f, press0, K1, p0=[1.0, 1.0],full_output=True)

perr = np.sqrt(np.diag(pcov))
E = (infodict['fvec']**2).sum()
r2 = 1-E/((K1-K1.mean())**2).sum()

print("Fit to K1(P)=K1(0)+a*P:")
print("R-squared =",r2)
print("K1(0)=", p[1]," J/m^3")
print("a=", p[0]," J/(m^3*Pa)")
print("(1/K1(0))*dK1/dP=", (p[0]/p[1])*1e9," 1/GPa")

plt.xlabel('P(GPa)'); plt.ylabel('K1(J/m^3)')
plt.plot(press0*1e-9, K1, 'ro', label="data")
plt.plot(xdib*1e-9, f(xdib, *p), label="fit")
plt.legend()
plt.savefig("fit_dK1dP.png")
plt.show()

f = open('K1vsP.dat','w')
f.write("P(kbar) K1(kJ/m^3)\n")
for i in range(len(press0)):
  f.write(str(press0[i]*1e-5*1e-3))
  f.write(" ")
  f.write(str(K1[i]*1e-3))
  f.write("\n")
f.close









