#!/bin/bash

#bash file to iterate through the ising model code - not on the HPCC
#change the temperature in the input file

iFile="input.in"
mFile="mag_temp.txt"

if [ -f "$iFile" ] && [ -f "$mFile" ];
then 
   #change the temperature
   in1="5.3		!temperature";
   in2="5_3.txt		!output file name (to read into stats.f)"
   
   for (( c=0; c<=40; c=c+1 ))
   do 
      #turn c into a decimal to iterate temp by 0.1
      d=$(echo "scale=2; $c/10.0" | bc)
      #replace temp in input file
      out1=$d"		!temperature"
      out2="0"${d/./_}".txt		!output file name (to read into stats.f)"
      #echo "$out2"
      sed -i "s/$in1/$out1/g" $iFile
	  sed -i "s/$in2/$out2/g" $iFile
      ~/My\ Documents/CompPhys/Ising/ising < $iFile
      ~/My\ Documents/CompPhys/Ising/stats < $iFile
      sed -i "s/$out1/$in1/g" $iFile
	  sed -i "s/$out2/$in2/g" $iFile

      #read the temp, avg. m, and std. dev. from mFile
      #and compile into one file
      while read line
      do 
        echo "$line" >> mag_file.txt
	echo "$line"
      done < $mFile
   done
else 
   echo "Input file $iFile or $mFile does not exist, exiting"
   exit
fi
