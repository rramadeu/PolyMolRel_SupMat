#!/bin/bash

## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## Computes the relationship matrix considering the biallelic methods (VR, PD, FA) for all the files from 2 using AGHmatrix V2.0

library(AGHmatrix)
system("mkdir biallelic_results")
system("mkdir biallelic_results/2ploidy")
system("mkdir biallelic_results/4ploidy")
system("mkdir biallelic_results/6ploidy")
system("mkdir biallelic_results/8ploidy")
system("mkdir biallelic_results/4ploidy_natural")
system("mkdir biallelic_results/6ploidy_natural")
system("mkdir biallelic_results/8ploidy_natural")

loci <- c(5,10,15,20,50,100,500,1000)
reptotal <- 3

for(loci.n in 1:length(loci)){
    Time = proc.time()
    mat.array.VanRaden.2 <- mat.array.Yang.2 <- mat.array.Slater.2 <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.4 <- mat.array.Yang.4 <- mat.array.Slater.4 <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.6 <- mat.array.Yang.6 <- mat.array.Slater.6 <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.8 <- mat.array.Yang.8 <- mat.array.Slater.8 <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.4.natural <- mat.array.Yang.4.natural <- mat.array.Slater.4.natural <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.6.natural <- mat.array.Yang.6.natural <- mat.array.Slater.6.natural <- array(NA,dim=c(183,183,reptotal))
    mat.array.VanRaden.8.natural <- mat.array.Yang.8.natural <- mat.array.Slater.8.natural <- array(NA,dim=c(183,183,reptotal))
    
    for(rep in 1:reptotal){
        ## Ploidy2
        load(paste0("biallelicdata/mat_",rep,".Rdata"))
        temp <- t(biallelic.ploidy2[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.2[,,rep] <- Gmatrix(temp,
                                               ploidy=2,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)
?Gmatrix
        mat.array.Yang.2[,,rep] <- Gmatrix(temp,
                                           ploidy=2,
                                           method="Yang",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        ## Ploidy4
        temp <- t(biallelic.ploidy4[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.4[,,rep] <- Gmatrix(temp,
                                               ploidy=4,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.4[,,rep] <- Gmatrix(temp,
                                           ploidy=4,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.4[,,rep] <- Gmatrix(temp,
                                           ploidy=4,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)
        
        ## Ploidy6
        temp <- t(biallelic.ploidy6[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.6[,,rep] <- Gmatrix(temp,
                                               ploidy=6,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.6[,,rep] <- Gmatrix(temp,
                                           ploidy=6,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.6[,,rep] <- Gmatrix(temp,
                                           ploidy=6,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        ## Ploidy8
        temp <- t(biallelic.ploidy8[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.8[,,rep] <- Gmatrix(temp,
                                               ploidy=8,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.8[,,rep] <- Gmatrix(temp,
                                           ploidy=8,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.8[,,rep] <- Gmatrix(temp,
                                           ploidy=8,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        ## Ploidy4
        temp <- t(biallelic.ploidy4.natural[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.4.natural[,,rep] <- Gmatrix(temp,
                                               ploidy=4,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.4.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=4,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.4.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=4,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)
        
        ## Ploidy6
        temp <- t(biallelic.ploidy6.natural[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.6.natural[,,rep] <- Gmatrix(temp,
                                               ploidy=6,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.6.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=6,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.6.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=6,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        ## Ploidy8
        temp <- t(biallelic.ploidy8.natural[sample(1:10000,loci[loci.n]),])
        mat.array.VanRaden.8.natural[,,rep] <- Gmatrix(temp,
                                               ploidy=8,
                                               thresh.missing=1,
                                               verify.posdef=FALSE)

        mat.array.Yang.8.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=8,
                                           method="Yang",
                                           pseudo.diploid=TRUE,
                                           thresh.missing=1,
                                           verify.posdef=FALSE)

        mat.array.Slater.8.natural[,,rep] <- Gmatrix(temp,
                                           ploidy=8,
                                           method="Slater",
                                           thresh.missing=1,
                                           verify.posdef=FALSE)
    }
    
    save(mat.array.VanRaden.2,file=paste0("biallelic_results/2ploidy/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.4,file=paste0("biallelic_results/4ploidy/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.6,file=paste0("biallelic_results/6ploidy/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.8,file=paste0("biallelic_results/8ploidy/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.4.natural,file=paste0("biallelic_results/4ploidy_natural/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.6.natural,file=paste0("biallelic_results/6ploidy_natural/VanRaden_loci",loci[loci.n],".Rdata"))
    save(mat.array.VanRaden.8.natural,file=paste0("biallelic_results/8ploidy_natural/VanRaden_loci",loci[loci.n],".Rdata"))

    save(mat.array.Yang.2,file=paste0("biallelic_results/2ploidy/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.4,file=paste0("biallelic_results/4ploidy/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.6,file=paste0("biallelic_results/6ploidy/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.8,file=paste0("biallelic_results/8ploidy/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.4.natural,file=paste0("biallelic_results/4ploidy_natural/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.6.natural,file=paste0("biallelic_results/6ploidy_natural/Yang_loci",loci[loci.n],".Rdata"))
    save(mat.array.Yang.8.natural,file=paste0("biallelic_results/8ploidy_natural/Yang_loci",loci[loci.n],".Rdata"))

    save(mat.array.Slater.2,file=paste0("biallelic_results/2ploidy/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.4,file=paste0("biallelic_results/4ploidy/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.6,file=paste0("biallelic_results/6ploidy/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.8,file=paste0("biallelic_results/8ploidy/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.4.natural,file=paste0("biallelic_results/4ploidy_natural/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.6.natural,file=paste0("biallelic_results/6ploidy_natural/Slater_loci",loci[loci.n],".Rdata"))
    save(mat.array.Slater.8.natural,file=paste0("biallelic_results/8ploidy_natural/Slater_loci",loci[loci.n],".Rdata"))
    
    Time = as.matrix(proc.time() - Time)
    cat("Completed! Time =", Time[3]/60, " minutes \n Loci", loci[loci.n])
}
