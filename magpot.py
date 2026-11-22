'''
MAGPOT version 1.0
'''

from numpy import *

J2eV = 6.241509e18       # Convert Joule to eV
kB = 1.38064852e-23      # Boltzmann constant in J/K

print("#"*30)
print("MAGPOT code (version 1.0)")
print("#"*30)
print("Authors: P.Nieves, I.Korniienko, S.Arapan, A.Fraile, R.Iglesias and D.Legut")

print("Reading data ...")
data = loadtxt('input.dat',dtype=str)
data = atleast_1d(data)

if data.size != 11:
    raise ValueError("input.dat must contain exactly 11 numerical values")

p=data[0]
print("Program output format =",p)
c=data[1]
print("Crystal system =",c)
m=data[2]
print("Magnetic order =",m)
r0 = float(data[3])
print("Nearest neighbor distance at equilibrium volume without magnetic interactions in units of A =", r0)
B0 = float(data[4])
print("Bulk modulus = ", B0, " GPa")
elas=3.0*B0*1e9 # c11+2*c12 = 3*B0 (in Pa)
Tc = float(data[5])
print("Curie or Neel temperature = ", Tc, " K")
ws = float(data[6])
print("Spontaneous volume magnetostriction = ",ws)
K1 = float(data[7])
print("First magnetocrystalline anisotropy K1 = ",K1," J/m^3")
dK1 = float(data[8])
print("(1/K1)*(dK1/dP) = ",dK1, " 1/GPa")
b1 = float(data[9])
print("First anisotropic magnetoelastic constant b1 = ",b1, " MPa")
b1 = b1*1e6 # convert to Pa
b2 = float(data[10])
print("Second anisotropic magnetoelastic constant b2 = ",b2, " MPa")
b2 = b2*1e6 # convert to Pa

if r0 <= 0.0 or B0 <= 0.0 or Tc <= 0.0:
    raise ValueError("r0, B0 and Tc must be positive")
if c not in ('sc', 'bcc', 'fcc') or m not in ('fm', 'afm'):
    raise ValueError("Use c = 'sc', 'bcc' or 'fcc' and m = 'fm' or 'afm'")
if p!='lmp':
    raise ValueError("Current version only supports LAMMPS (lmp) output format") 

# ----------------------------------------------------------------------
# Simple cubic (SC)
# ----------------------------------------------------------------------
if c == 'sc' and m == 'afm':
    d1 = r0
    d2 = r0*sqrt(2.0)
    dcut = (d1+d2)/2.0
    n = 1.0
    v0 = (r0*1e-10)**3

    L = ((v0*b2)/(2.0*n))*J2eV
    dL = ((v0*b1)/(n))*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = r0

    q = ((K1*5.0*v0)/(26.0*n))*J2eV
    dq = ((K1*15.0*v0)/(26.0*n))*J2eV*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = r0

    J = -0.5*kB*Tc*J2eV
    dJ = -((ws*elas*v0)/(3.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = r0

elif c == 'sc' and m == 'fm':
    d1 = r0
    d2 = r0*sqrt(2.0)
    dcut = (d1+d2)/2.0
    n = 1.0
    v0 = (r0*1e-10)**3

    L = -((v0*b2)/(2.0*n))*J2eV
    dL = -((v0*b1)/(n))*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = r0

    q = ((K1*v0)/(2.0*n))*J2eV
    dq = 3.0*q*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = r0

    J = 0.5*kB*Tc*J2eV
    dJ = ((ws*elas*v0)/(3.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = r0

# ----------------------------------------------------------------------
# Body-centered cubic (BCC)
# ----------------------------------------------------------------------
elif c == 'bcc' and m == 'fm':
    a0 = 2.0*r0/sqrt(3.0)
    d1 = r0
    d2 = a0
    dcut = (d1+d2)/2.0
    n = 2.0
    v0 = (a0*1e-10)**3

    L = -(3.0*v0*b1/(8.0*n))*J2eV
    dL = (3.0*v0*(b1-3.0*b2)/(8.0*n))*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = r0

    q = -(9.0*v0*K1/(16.0*n))*J2eV
    dq = 3.0*q*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = r0

    J = (3.0*kB*Tc/8.0)*J2eV
    dJ = (ws*elas*v0/(4.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = r0

elif c == 'bcc' and m == 'afm':
    a0 = 2.0*r0/sqrt(3.0)
    d1 = r0
    d2 = a0
    dcut = (d1+d2)/2.0
    n = 2.0
    v0 = (a0*1e-10)**3

    # BCC-AFM type G: all eight first neighbours are antiparallel.
    L = (3.0*v0*b1/(8.0*n))*J2eV
    dL = -(3.0*v0*(b1-3.0*b2)/(8.0*n))*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = r0

    # Quadrupole term obtained from the AFM angular expression.
    q = -(45.0*v0*K1/(208.0*n))*J2eV
    dq = 3.0*q*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = r0

    J = -(3.0*kB*Tc/8.0)*J2eV
    dJ = -(ws*elas*v0/(4.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = r0

# ----------------------------------------------------------------------
# Face-centered cubic (FCC)
# ----------------------------------------------------------------------
elif c == 'fcc' and m == 'fm':
    a0 = sqrt(2.0)*r0
    d1 = r0
    d2 = sqrt(2.0)*r0
    dcut = (d1+d2)/2.0
    n = 4.0
    v0 = (a0*1e-10)**3

    L = (v0/(2.0*n))*(b2/2.0-b1)*J2eV
    dL = (v0/n)*(b1-1.5*b2)*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = r0

    q = -(v0*K1/n)*J2eV
    dq = 3.0*q*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = r0

    J = (kB*Tc/4.0)*J2eV
    dJ = (ws*elas*v0/(6.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = r0

elif c == 'fcc' and m == 'afm':
    # FCC-AFM type II: use only the six antiparallel second neighbours.
    # These neighbours form an SC magnetic sublattice.
    a0 = sqrt(2.0)*r0
    d1 = a0
    d2 = sqrt(2.0)*d1
    dcut = (d1+d2)/2.0
    n = 4.0
    v0 = (a0*1e-10)**3

    L = ((v0*b2)/(2.0*n))*J2eV
    dL = ((v0*b1)/(n))*J2eV
    alpha1 = (exp(1)/8.0)*(2.0*L-dL)
    gamma1 = -(dL)/(2.0*L-dL)
    delta1 = d1

    q = ((K1*5.0*v0)/(26.0*n))*J2eV
    dq = ((K1*15.0*v0)/(26.0*n))*J2eV*(1.0-B0*dK1)
    alpha2 = (exp(1)/8.0)*(2.0*q-dq)
    gamma2 = -(dq)/(2.0*q-dq)
    delta2 = d1

    J = -0.5*kB*Tc*J2eV
    dJ = -((ws*elas*v0)/(3.0*n))*J2eV
    alpha3 = (exp(1)/8.0)*(2.0*J-dJ)
    gamma3 = -(dJ)/(2.0*J-dJ)
    delta3 = d1

print("-"*40)
print("RESULTS FOR BETHE-SLATER PARAMETERS:")
print("-"*40)

# For FCC-AFM type II, different atom types are needed for the four FCC
# magnetic sublattices. This prevents the first-neighbour pairs from interacting.
if c == 'fcc' and m == 'afm':
    print("\nWARNING: FCC-AFM (type-II antiferromagnetic order)")
    print("-"*55)
    print("The 12 nearest neighbours cannot be all antiparallel.")
    print("For type-II AFM ordering, 6 nearest neighbours are")
    print("parallel and 6 are antiparallel.")
    print("")
    print("MAGPOT considers only the 6 second-nearest neighbours,")
    print("which are all antiparallel.")
    print("These interactions form four interpenetrating SC")
    print("sublattices within the FCC conventional unit cell.")
    print("")
    print("Therefore, the SC-AFM expressions are used for J(r),")
    print("l(r) and q(r), considering all four SC sublattices.")
    print("Interactions between different SC sublattices are neglected.")
    print("-"*55)
    print("pair_coeff  * * spin/exchange exchange ", dcut, 0.0, 0.0, delta3, " offset no")
    for i in range(1, 5):
        print("pair_coeff ", i, i, " spin/exchange exchange ", dcut, alpha3, gamma3, delta3, " offset no")
    print("pair_coeff  * * spin/neel neel ", dcut, 0.0, 0.0, delta1, 0.0, 0.0, delta2)
    for i in range(1, 5):
        print("pair_coeff ", i, i, " spin/neel neel ", dcut, alpha1, gamma1, delta1, alpha2, gamma2, delta2)
else:
    print("pair_coeff  * * spin/exchange exchange ", dcut, alpha3, gamma3, delta3, " offset no")
    print("pair_coeff  * * spin/neel neel ", dcut, alpha1, gamma1, delta1, alpha2, gamma2, delta2)

    
print(" ")
