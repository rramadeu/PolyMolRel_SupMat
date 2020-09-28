#!/bin/bash

## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## Computes the relationship matrix considering the multiallelic methods (LO, RI, WE, MM, ML) for all the files from 3 using PolyRelatedness V1.8

#Clonning data structure for the output
cp -r multiallelicdata PolyRelatedness_output_LO
cp -r multiallelicdata PolyRelatedness_output_RI
cp -r multiallelicdata PolyRelatedness_output_WE
cp -r multiallelicdata PolyRelatedness_output_MM
cp -r multiallelicdata PolyRelatedness_output_ML

#Listing files to loop
find multiallelicdata -name \*.txt -print > PolyRelatednessInput.txt
find PolyRelatedness_output_LO -name \*.txt -print > PolyRelatednessOutput_LO.txt
find PolyRelatedness_output_RI -name \*.txt -print > PolyRelatednessOutput_RI.txt
find PolyRelatedness_output_WE -name \*.txt -print > PolyRelatednessOutput_WE.txt
find PolyRelatedness_output_MM -name \*.txt -print > PolyRelatednessOutput_MM.txt
find PolyRelatedness_output_ML -name \*.txt -print > PolyRelatednessOutput_ML.txt

#Saving file names in a array
readarray -t in < PolyRelatednessInput.txt
readarray -t out_LO < PolyRelatednessOutput_LO.txt
readarray -t out_RI < PolyRelatednessOutput_RI.txt
readarray -t out_WE < PolyRelatednessOutput_WE.txt
readarray -t out_MM < PolyRelatednessOutput_MM.txt
readarray -t out_ML < PolyRelatednessOutput_ML.txt

#Doing the loop
len=${#in[@]}
echo $len

for (( i=0; i<$len; i++ ))
do
	rm ${out_LO[$i]}
	rm ${out_RI[$i]}
	rm ${out_WE[$i]}
	rm ${out_MM[$i]}
	rm ${out_ML[$i]}
done


time for (( i=0; i<$len; i++ ))
do
	echo ${in[$i]} 
	
	echo ${out_LO[$i]}
	./PolyRelatedness_1.8/PolyRelatedness.out ${in[$i]} ${out_LO[$i]} e 6
	
	echo ${out_RI[$i]}
	./PolyRelatedness_1.8/PolyRelatedness.out ${in[$i]} ${out_RI[$i]} e 4
	
	echo ${out_WE[$i]}
	./PolyRelatedness_1.8/PolyRelatedness.out ${in[$i]} ${out_WE[$i]} e 7
	
	echo ${out_MM[$i]}
	./PolyRelatedness_1.8/PolyRelatedness.out ${in[$i]} ${out_MM[$i]} e 1
	
	echo ${out_ML[$i]}
	./PolyRelatedness_1.8/PolyRelatedness.out ${in[$i]} ${out_ML[$i]} e 2
done


rm PolyRelatednessInput.txt
rm PolyRelatednessOutput_LO.txt
rm PolyRelatednessOutput_RI.txt
rm PolyRelatednessOutput_WE.txt
rm PolyRelatednessOutput_MM.txt
rm PolyRelatednessOutput_ML.txt

