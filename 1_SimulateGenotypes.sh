#!/bin/bash

## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## Loops the PedigreeSim software through all the parameters files N times (here, N=3)

cd PedigreeSimInput
time for i in {1..3} #increase 3 to number of reps you want to analyse
do
	java -jar ../PedigreeSim/PedigreeSim.jar 2.par
	java -jar ../PedigreeSim/PedigreeSim.jar 4.par
	java -jar ../PedigreeSim/PedigreeSim.jar 6.par
	java -jar ../PedigreeSim/PedigreeSim.jar 8.par
	java -jar ../PedigreeSim/PedigreeSim.jar 4_natural.par
	java -jar ../PedigreeSim/PedigreeSim.jar 6_natural.par
	java -jar ../PedigreeSim/PedigreeSim.jar 8_natural.par
	rm *_out_genotypes* *_out_alleledose* *hsa *hsb
	mkdir ../output/$i
	mv *_out_* ../output/$i
done 
