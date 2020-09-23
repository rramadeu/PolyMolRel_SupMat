# Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops

Supplementary Material (scripts) to reproduce all the analysis of the manuscript:

**"Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"**

by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia

Journal/Issue/Year: TBA

All the scripts were evaluated under Linux Ubuntu 20.04. The pipeline of analysis uses mainly bash, R, and java.

To reproduce our results, you need to run the scripts in order. For a comprehensive understanding, I order the scripts by 1 to N, they should be run in order fashion. It follows a short description of the content and steps:

## Folders:
`PedigreeSimInput` folder with the necessary parameters for PedigreeSim software. The files starting with 0 are common for all the ploidies:
- `0.chrom`: general genomic information (number of chromosomes, length, centromere position)
- `0.map`: genetic map with marker names and position
- `0.ped`: the pedigree used in this study
- `X.gen`: we are not using this information during the analysis, you can ignore it, but it is necessary for PedigreeSim to run. `X` is the ploidy
- `X.par`: parameters for the simulation
- `X_natural.par`: parameters for the simulation considering natural pairing


**PedigreeSim** folder with the PedigreeSim software V2.0

## Scripts:
