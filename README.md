![GitHub Logo](logo.png)

MagPot code v1.0

Authors: P.Nieves, I.Korniienko, S.Arapan, A.Fraile, R.Iglesias and D.Legut

-------------------------
WHAT IS MAGPOT?
-------------------------

```MagPot``` is a program to generate magnetoelastic interatomic potentials for molecular dynamics simulations capable to reproduce the desired Curie (or Neel) temperature, spontaneous volume magnetostriction (w_s), pressure dependence of magnetocrystalline anisotropy K1(P) and anisotropic magnetoelastic constants (b). Current version only generates the interatomic potential with format suitable for the molecular dynamics code ```LAMMPS```.

------------------
USAGE:
------------------

The ```MagPot``` code requires to have ```Python3```. It can be easily used by running the python script ```magpot.py``` as
```bash
python3 magpot.py > results.dat
```
where the data of the generated potential will be written in file ```results.dat```. The settings of this program are specified throught the file ```input.dat``` using the same format and units as the provided example. The file ```input.dat``` should be in the same folder where you are running the python script ```magpot.py```. See more details in the User Manual.

