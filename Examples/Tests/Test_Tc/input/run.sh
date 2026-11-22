#! /bin/bash 
for tempe in 25 100 200 300 400 500 600 700 800 900
do
 mkdir T$tempe/ 2>/dev/null
 echo "Temp: " $tempe
 cp ./* T$tempe/
 cd T$tempe/
 
 if [ $tempe -eq 25 ]
 then
   sed -i "s/tempe/$tempe/g" in.Fe0
   lmp_serial < in.Fe0 
   cp ./restart_Fe.equil ../restart_Fe.equil0
 fi
 if [ $tempe -ne 25 ]
 then
   sed -i "s/tempe/$tempe/g" in.Fe
   cp ../restart_Fe.equil0 ./   
   lmp_serial < in.Fe 
   cp ./restart_Fe.equil ../restart_Fe.equil0
 fi

 cd ..
done


echo "# Temp.   Mavg    " > thermal.norm.dat
for Temp in 25 100 200 300 400 500 600 700 800 900
 do
          tail -1 T$Temp/output.dat | awk '{ printf "%6i  %f\n", temp, $7 }' temp=$Temp >> thermal.norm.dat
 done


python3 fit.py > results.dat
