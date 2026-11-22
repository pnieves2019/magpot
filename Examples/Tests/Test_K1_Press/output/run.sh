#!/bin/bash

inputlmpfm=in.Fe_110   # input file name for LAMMPS spin direction 110

inputlmpfm2=in.Fe_100   # input file name for LAMMPS spin direction 100

exelmp=lmp_serial  # name of the executable file for LAMMPS


p1=995
p2=997
p3=999
p4=1000
p5=1001
p6=1003
p7=1005

for i in ${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7}

do

  cp ${inputlmpfm} in.lmp_0


  sed -i "s/AA/$i/"  ./in.lmp_0

  ${exelmp} < in.lmp_0 > results_fm_${i}

done


cat energy_${p1}.dat energy_${p2}.dat energy_${p3}.dat energy_${p4}.dat energy_${p5}.dat energy_${p6}.dat energy_${p7}.dat > ene_vs_vol_fm_110.dat


rm ./energy_*


for i in ${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7}

do

  cp ${inputlmpfm2} in.lmp_0


  sed -i "s/AA/$i/"  ./in.lmp_0

  ${exelmp} < in.lmp_0 > results_fm2_${i}

done


cat energy_${p1}.dat energy_${p2}.dat energy_${p3}.dat energy_${p4}.dat energy_${p5}.dat energy_${p6}.dat energy_${p7}.dat > ene_vs_vol_fm_100.dat


rm ./energy_*



python3 fit.py > output.dat
