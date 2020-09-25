## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## computes the metrics between observed and estimated relationships and save them in a table format.

library(DescTools) #for CCC coefficient
reps = 3

## Creating a data.frame to organize the results
stat.PolyRelatedness <- expand.grid(
    methods = c("LO","RI","WE","MM","ML"),
    ploidy = c("2ploidy","4ploidy","6ploidy","8ploidy","4ploidy_natural","6ploidy_natural","8ploidy_natural"),
    loci = c("5","10","15","20","50","100","500","1000"),
    alleles = c("2alleles","3alleles","4alleles","5alleles","10alleles","15alleles","20alleles"),
    RMSE.PO = NA,
    RMSE.FS = NA,
    RMSE.HS = NA,
    RMSE.GPGO = NA,
    RMSE.UN = NA,
    RMSE.GUGN = NA,
    RMSE.UR = NA,
    IC.PO = NA,
    IC.FS = NA,
    IC.HS = NA,
    IC.GPGO = NA,
    IC.UN = NA,
    IC.GUGN = NA,
    IC.UR = NA,
    meancor.Farthing = NA,
    meancor.IAC = NA,
    RMSE.Farthing = NA,
    RMSE.IAC = NA,
    CCC.Farthing = NA,
    CCC.IAC = NA,
    CCC.lwrci.Farthing = NA,
    CCC.uprci.Farthing = NA,
    CCC.lwrci.IAC = NA,
    CCC.uprci.IAC = NA)

stat.AGHmatrix <- expand.grid(
    methods = c("VanRaden","Yang","Slater"),
    ploidy = c("2ploidy","4ploidy","6ploidy","8ploidy","4ploidy_natural","6ploidy_natural","8ploidy_natural"),
    loci = c("5","10","15","20","50","100","500","1000"),
    alleles = c("2alleles"),
    RMSE.PO = NA,
    RMSE.FS = NA,
    RMSE.HS = NA,
    RMSE.GPGO = NA,
    RMSE.UN = NA,
    RMSE.GUGN = NA,
    RMSE.UR = NA,
    IC.PO = NA,
    IC.FS = NA,
    IC.HS = NA,
    IC.GPGO = NA,
    IC.UN = NA,
    IC.GUGN = NA,
    IC.UR = NA,
    meancor.Farthing = NA,
    meancor.IAC = NA,
    RMSE.Farthing = NA,
    RMSE.IAC = NA,
    CCC.Farthing = NA,
    CCC.IAC = NA,
    CCC.lwrci.Farthing = NA,
    CCC.uprci.Farthing = NA,
    CCC.lwrci.IAC = NA,
    CCC.uprci.IAC = NA)


## Extracting observed relatedness
obs.r = list()
load("ObservedRelatedness/r2ploidy.Rdata")
obs.r$'2ploidy' = mat.array
load("ObservedRelatedness/r4ploidy.Rdata")
obs.r$'4ploidy' = mat.array
load("ObservedRelatedness/r6ploidy.Rdata")
obs.r$'6ploidy' = mat.array
load("ObservedRelatedness/r8ploidy.Rdata")
obs.r$'8ploidy' = mat.array
load("ObservedRelatedness/rnatural4ploidy.Rdata")
obs.r$'4ploidy_natural' = mat.array
load("ObservedRelatedness/rnatural6ploidy.Rdata")
obs.r$'6ploidy_natural' = mat.array
load("ObservedRelatedness/rnatural8ploidy.Rdata")
obs.r$'8ploidy_natural' = mat.array

## Extracting results and computing metrics for PolyRelatedness metrics
## "LO","RI","WE","MM","ML"
for(i in 1:nrow(stat.PolyRelatedness)){
    print(paste0(i,"/",nrow(stat.PolyRelatedness)))
    est.r = obs.r$'2ploidy'*0 #creates an empity array

    ## loop to extract the PolyRelatedness results
    for(j in 1:reps){
        mat = read.table(paste0("PolyRelatedness_output_",stat.PolyRelatedness$methods[i],"/",
                                stat.PolyRelatedness$ploidy[i],"/",
                                stat.PolyRelatedness$loci[i],"loci/",
                                stat.PolyRelatedness$alleles[i],"/",
                                j,".txt"),skip=7)
        mat = as.matrix(mat)
        est.r[,,j] = mat
    }

    ploidy = as.character(stat.PolyRelatedness$ploidy[i])

    ## gather observed relatedness (r)
    PO.obs <- obs.r[[ploidy]][181,182,1:reps]
    FS.obs <- obs.r[[ploidy]][179,180,1:reps]
    HS.obs <- obs.r[[ploidy]][180,181,1:reps]
    GPGO.obs <- obs.r[[ploidy]][181,183,1:reps]
    UN.obs <- obs.r[[ploidy]][180,182,1:reps]
    GUGN.obs <- obs.r[[ploidy]][180,183,1:reps]
    Farthing.obs <- obs.r[[ploidy]][74:149,150,1:reps]
    IAC.obs <- obs.r[[ploidy]][151:176,177,1:reps]

    ## gather estimated relatedness (r)
    PO.est <- est.r[181,182,1:reps]
    FS.est <- est.r[179,180,1:reps]
    HS.est <- est.r[180,181,1:reps]
    GPGO.est <- est.r[181,183,1:reps]
    UN.est <- est.r[180,182,1:reps]
    UR.est <- est.r[1:10,1:10,1:reps] #the first 10 unrelated
    GUGN.est <- est.r[180,183,1:reps]
    Farthing.est <- est.r[74:149,150,1:reps]
    IAC.est <- est.r[151:176,177,1:reps]
    
    ## statistic for common relationships
    ## PO: parent-offspring, FS: fullsibs, HS: halfsibs
    ## GPGO: grandparent-grandoffspring, UN: uncle-nephew
    ## GUGN: granduncle-grandnephew, UR: unrelated (first 10 ancestrals)

    ## RMSE: root means square error
    stat.PolyRelatedness$RMSE.PO[i] = sqrt(mean(PO.est-PO.obs)^2)
    stat.PolyRelatedness$RMSE.FS[i] = sqrt(mean(FS.est-FS.obs)^2)
    stat.PolyRelatedness$RMSE.HS[i] = sqrt(mean(HS.est-HS.obs)^2)
    stat.PolyRelatedness$RMSE.GPGO[i] = sqrt(mean(GPGO.est-GPGO.obs)^2)
    stat.PolyRelatedness$RMSE.UN[i] = sqrt(mean(UN.est-UN.obs)^2)
    stat.PolyRelatedness$RMSE.GUGN[i] = sqrt(mean(GUGN.est-GUGN.obs)^2)
    stat.PolyRelatedness$RMSE.UR[i] = sqrt(mean(UR.est)^2) #observed is 0

    ## IC: Interval of confidence (how many values are within +- of the observed value)
    stat.PolyRelatedness$IC.PO[i] = mean(abs(PO.est-PO.obs)<0.05)
    stat.PolyRelatedness$IC.FS[i] = mean(abs(FS.est-FS.obs)<0.05)
    stat.PolyRelatedness$IC.HS[i] = mean(abs(HS.est-HS.obs)<0.05)
    stat.PolyRelatedness$IC.GPGO[i] = mean(abs(GPGO.est-GPGO.obs)<0.05)
    stat.PolyRelatedness$IC.UN[i] = mean(abs(UN.est-UN.obs)<0.05)
    stat.PolyRelatedness$IC.GUGN[i] = mean(abs(GUGN.est-GUGN.obs)<0.05)
    stat.PolyRelatedness$IC.UR[i] = mean(abs(UR.est)<0.05)

    ## Pearson's correlation for each pedigree (Farthing and IAC)
    stat.PolyRelatedness$meancor.Farthing[i] =  cor(as.vector(Farthing.obs),
                                                    as.vector(Farthing.est),method="pearson")
    stat.PolyRelatedness$meancor.IAC[i] =  cor(as.vector(IAC.obs),
                                               as.vector(IAC.est),method="pearson")
    
    ## RMSE for each pedigree (Farthing and IAC)
    stat.PolyRelatedness$RMSE.Farthing[i] =  sqrt(mean((as.vector(Farthing.obs) -
                                                        as.vector(Farthing.est))^2))
    stat.PolyRelatedness$RMSE.IAC[i] =  sqrt(mean((as.vector(IAC.obs) -
                                                   as.vector(IAC.est))^2))

    ## Lin's CCC for each pedigree (Farthing and IAC)
    stat.PolyRelatedness$CCC.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[1])
    stat.PolyRelatedness$CCC.lwrci.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[2])
    stat.PolyRelatedness$CCC.uprci.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[3])
    
    stat.PolyRelatedness$CCC.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                      as.vector(IAC.est))$rho.c[1])  
    stat.PolyRelatedness$CCC.lwrci.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                           as.vector(IAC.est))$rho.c[2])
    stat.PolyRelatedness$CCC.uprci.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                            as.vector(IAC.est))$rho.c[3])
}

## Extracting results and computing metrics for AGHmatrix metrics
## "VR","PD","FA" 
for(i in 1:nrow(stat.AGHmatrix)){
    print(paste0(i,"/",nrow(stat.AGHmatrix)))
     ## extract the AGHmatrix results
    mat.name = load(paste0("biallelic_results/",stat.AGHmatrix$ploidy[i],"/",
                           stat.AGHmatrix$methods[i],"_loci",stat.AGHmatrix$loci[i],".Rdata"))
    est.r = get(mat.name)
    rm(mat.name)

    ploidy = as.character(stat.AGHmatrix$ploidy[i])

    ## gather observed relatedness (r)
    PO.obs <- obs.r[[ploidy]][181,182,1:reps]
    FS.obs <- obs.r[[ploidy]][179,180,1:reps]
    HS.obs <- obs.r[[ploidy]][180,181,1:reps]
    GPGO.obs <- obs.r[[ploidy]][181,183,1:reps]
    UN.obs <- obs.r[[ploidy]][180,182,1:reps]
    GUGN.obs <- obs.r[[ploidy]][180,183,1:reps]
    Farthing.obs <- obs.r[[ploidy]][74:149,150,1:reps]
    IAC.obs <- obs.r[[ploidy]][151:176,177,1:reps]

    ## gather estimated relatedness (r)
    PO.est <- est.r[181,182,1:reps]
    FS.est <- est.r[179,180,1:reps]
    HS.est <- est.r[180,181,1:reps]
    GPGO.est <- est.r[181,183,1:reps]
    UN.est <- est.r[180,182,1:reps]
    UR.est <- est.r[1:10,1:10,1:reps] #the first 10 unrelated
    GUGN.est <- est.r[180,183,1:reps]
    Farthing.est <- est.r[74:149,150,1:reps]
    IAC.est <- est.r[151:176,177,1:reps]
    
    ## statistic for common relationships
    ## PO: parent-offspring, FS: fullsibs, HS: halfsibs
    ## GPGO: grandparent-grandoffspring, UN: uncle-nephew
    ## GUGN: granduncle-grandnephew, UR: unrelated (first 10 ancestrals)

    ## RMSE: root means square error
    stat.AGHmatrix$RMSE.PO[i] = sqrt(mean(PO.est-PO.obs)^2)
    stat.AGHmatrix$RMSE.FS[i] = sqrt(mean(FS.est-FS.obs)^2)
    stat.AGHmatrix$RMSE.HS[i] = sqrt(mean(HS.est-HS.obs)^2)
    stat.AGHmatrix$RMSE.GPGO[i] = sqrt(mean(GPGO.est-GPGO.obs)^2)
    stat.AGHmatrix$RMSE.UN[i] = sqrt(mean(UN.est-UN.obs)^2)
    stat.AGHmatrix$RMSE.GUGN[i] = sqrt(mean(GUGN.est-GUGN.obs)^2)
    stat.AGHmatrix$RMSE.UR[i] = sqrt(mean(UR.est)^2) #observed is 0

    ## IC: Interval of confidence (how many values are within +- of the observed value)
    stat.AGHmatrix$IC.PO[i] = mean(abs(PO.est-PO.obs)<0.05)
    stat.AGHmatrix$IC.FS[i] = mean(abs(FS.est-FS.obs)<0.05)
    stat.AGHmatrix$IC.HS[i] = mean(abs(HS.est-HS.obs)<0.05)
    stat.AGHmatrix$IC.GPGO[i] = mean(abs(GPGO.est-GPGO.obs)<0.05)
    stat.AGHmatrix$IC.UN[i] = mean(abs(UN.est-UN.obs)<0.05)
    stat.AGHmatrix$IC.GUGN[i] = mean(abs(GUGN.est-GUGN.obs)<0.05)
    stat.AGHmatrix$IC.UR[i] = mean(abs(UR.est)<0.05)

    ## Pearson's correlation for each pedigree (Farthing and IAC)
    stat.AGHmatrix$meancor.Farthing[i] =  cor(as.vector(Farthing.obs),
                                                    as.vector(Farthing.est),method="pearson")
    stat.AGHmatrix$meancor.IAC[i] =  cor(as.vector(IAC.obs),
                                               as.vector(IAC.est),method="pearson")
    
    ## RMSE for each pedigree (Farthing and IAC)
    stat.AGHmatrix$RMSE.Farthing[i] =  sqrt(mean((as.vector(Farthing.obs) -
                                                        as.vector(Farthing.est))^2))
    stat.AGHmatrix$RMSE.IAC[i] =  sqrt(mean((as.vector(IAC.obs) -
                                                   as.vector(IAC.est))^2))

    ## Lin's CCC for each pedigree (Farthing and IAC)
    stat.AGHmatrix$CCC.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[1])
    stat.AGHmatrix$CCC.lwrci.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[2])
    stat.AGHmatrix$CCC.uprci.Farthing[i] =  as.numeric(CCC(as.vector(Farthing.obs),
                                                           as.vector(Farthing.est))$rho.c[3])
    
    stat.AGHmatrix$CCC.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                      as.vector(IAC.est))$rho.c[1])  
    stat.AGHmatrix$CCC.lwrci.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                           as.vector(IAC.est))$rho.c[2])
    stat.AGHmatrix$CCC.uprci.IAC[i] =  as.numeric(CCC(as.vector(IAC.obs),
                                                            as.vector(IAC.est))$rho.c[3])
}

## merging results in a single data.frame and exporting in a csv
stats = rbind(stat.PolyRelatedness,stat.AGHmatrix)
write.csv(stats,"StatResults.csv",quote=FALSE,row.names=FALSE)




