#!/bin/bash

## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## Execute all the necessary script to perform the simulation and analysis

start=`date +%s`

bash 1_SimulateGenotypes.sh
Rscript 2_funderalleles2biallelic.R
Rscript 3_funderalleles2multiallelic.R
Rscript 4_ObservedRelatedness.R
bash 5_PolyRelatednes.sh
Rscript 6_AGHmatrix.R
Rscript 7_ProcessingResults.R

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

