# Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops

Supplementary Material (scripts) to reproduce all the analysis of the manuscript:

**"Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"**

by Rodrigo R Amadeu, Leticia de Castro Lara, Patricio R Munoz, Antonio Augusto Franco Garcia

Journal/Issue/Year: TBA

All the scripts were evaluated under Linux Ubuntu 20.04. The pipeline of analysis uses mainly R and bash scripts. For a comprehensive reproduction, I order the scripts by 1 to N, they should be run in order fashion. It follows a short description of the content and steps:

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
`.sh` are bash scripts and `.R` are R scripts
- `1_SimulateGenotypes.sh`: loops the PedigreeSim software through all the parameters files N times (here, N=3)
- `2_funderalleles2biallelic.R`: reads PedigreeSim output and sample biallelic markers (0 or 1) to all the loci. It saves the data as a matrix in Rdata format.
- `3_funderalleles2multiallelic.R`: reads PedigreeSim output and sample biallelic and multiallelic markers to all loci. It saves the data as a .txt following PolyRelatedness V1.8 format
- `4_PolyRelatednes.sh`: computes the relationship matrix considering the multiallelic methods (LO, RI, WE, MM, ML)  for all the files from 3 using PolyRelatedness V1.8
- `5_AGHmatrix.R`: computes the relationship matrix considering the biallelic methods (VR, PD, FA) for all the files from 2 using AGHmatrix V2.0
- `6_ObservedRelatedness.R`: computes the observed relationship matrix on the original data following delta computation (Equation 2 from the paper)
- `7_SummaryResults.R`: computes the metrics between observed and estimated relationships and save them in a table format.





