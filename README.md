# Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops

Supplementary Material (scripts) to reproduce all the analysis of the manuscript:

**"Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"**  
by Rodrigo R Amadeu, Leticia de Castro Lara, Patricio R Munoz, Antonio Augusto Franco Garcia  
G3: Genes, Genomes, Genetics, Early online October 13, 2020  
https://doi.org/10.1534/g3.120.401669  

## Before running:
All the scripts were evaluated under Linux Ubuntu 20.04 64 bits. Before run it, verify if `PedigreeSim` and `PolyRelatedness` software are fully working. For `PedigreeSim` you need to have Java Runtime Environment installed. `PedigreeSim` details at https://www.wur.nl/en/show/Software-PedigreeSim.htm. `PolyRelatedness` details at https://github.com/huangkang1987/polyrelatedness. In addition, you need to have `R` and `AGHmatrix` package installed. Details in https://cran.r-project.org/ and https://cran.r-project.org/package=AGHmatrix respectivally. 

This is a comprehensive set of scripts designed to be able to run in a personal computer. In practice for a higher number of replicates you would need to split the `for loops` into parallel computation. The details on this task depends upon your computational cluster specifications. 

Any question, please open an `Issue` here in this page or send me an e-mail (rramadeu at ufl dot edu). 

## Overall pipeline of analysis
It uses mainly R and bash scripts. `0_runAll.sh` is a wrap-up bash script to run all the steps. For a comprehensive reproduction, the scripts are in order from 1 to 7. It just loops for for **three** replicates (the original work has 100 replicates). You can extend the loops throughout the script for your aimed number. The whole pipeline considering three replicates took approximately 25 hours to run using a Intel® Core™ i7-8650U CPU @ 1.90GHz × 8 with 16 GB RAM. It follows a short description of the content and steps.

## Content
`PedigreeSimInput` folder with the necessary parameters for PedigreeSim software. The files starting with 0 are common for all the ploidies:
- `0.chrom`: general genomic information (number of chromosomes, length, centromere position)
- `0.map`: genetic map with marker names and position
- `0.ped`: the pedigree used in this study
- `X.gen`: we are not using this information during the analysis, you can ignore it, but it is necessary for PedigreeSim to run. `X` is the ploidy
- `X.par`: parameters for the simulation
- `X_natural.par`: parameters for the simulation considering natural pairing

`PedigreeSim` folder with the PedigreeSim software V2.0
`PolyRelatedness_1.8` folder with the PolyRelatedness software V1.8
`StatResults_Example.csv` output example with the summary statistics computed for each methodology for each population and replicate.

## Scripts:
`.sh` are bash scripts and `.R` are R scripts
- `1_SimulateGenotypes.sh`: loops the PedigreeSim software through all the parameters files N times (here, N=3)
- `2_funderalleles2biallelic.R`: reads PedigreeSim output and sample biallelic markers (0 or 1) to all the loci. It saves the data as a matrix in Rdata format.
- `3_funderalleles2multiallelic.R`: reads PedigreeSim output and sample biallelic and multiallelic markers to all loci. It saves the data as a .txt following PolyRelatedness V1.8 format
- `4_ObservedRelatedness.R`: computes the observed relationship matrix on the original data following delta computation (Equation 2 from the paper)
- `5_PolyRelatednes.sh`: computes the relationship matrix considering the multiallelic methods (LO, RI, WE, MM, ML)  for all the files from 3 using PolyRelatedness V1.8
- `6_AGHmatrix.R`: computes the relationship matrix considering the biallelic methods (VR, PD, FA) for all the files from 2 using AGHmatrix V2.0
- `7_ProcessingResults.R`: computes the metrics between observed and estimated relationships and save them in a table format.
