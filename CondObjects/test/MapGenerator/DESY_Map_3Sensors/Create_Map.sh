#!/bin/bash

#EE part

rm -rf Codes/*
rm -rf Executables/*
rm -rf map_CERN_Hexaboard_September_*.txt 

for i in `seq 1 3`
do
  if [ "$i" -eq 1 ]
  then
      cp Mapping_Hexaboard_L1_UnFlipped_Pi_EE.cpp Codes/Mapping_Layer"$i".cc
  fi

  if [ "$i" -eq 2 ]
  then
      cp Mapping_Hexaboard_L2_Flipped_MPiB2_EE.cpp Codes/Mapping_Layer"$i".cc
  fi

  if [ "$i" -eq 3 ]
  then
      cp Mapping_Hexaboard_L3_UnFlipped_PiB2_EE.cpp Codes/Mapping_Layer"$i".cc
  fi

  sed -i -e 's/FOO/'$i'/g' Codes/Mapping_Layer"$i".cc

done

#Create the executables
for i in `seq 1 3`
do
  g++ Codes/Mapping_Layer"$i".cc -o Executables/Mapping_Layer"$i"
  ./Executables/Mapping_Layer"$i"
done


for counter in `seq 1 3`
do
cat  map_CERN_Hexaboard_September_L_"$counter".txt >> map_CERN_Hexaboard_DESY_3Sensors_3EELayers_V1.txt
done
