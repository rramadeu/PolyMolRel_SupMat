## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## Reads PedigreeSim output and sample biallelic markers (0 or 1) to all the loci. It saves the data as a matrix in Rdata format.

system("mkdir biallelicdata")
indnames <- read.table("PedigreeSimInput/0.ped",header=TRUE,colClasses="character")[,1]

for(repindex in 1:3){
    ## Ploidy 2
    rawdata <- read.table(paste0("output/",repindex,"/2_out_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)
    
    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)
    
    ## Ploidy 2 Biallelic
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }

    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/2)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1
    while(i < 366){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1]
        j <- j+1
        i <- i+2
    }
    biallelic.ploidy2 <- biallelic.matrix.merged

    ## Ploidy 4
    rawdata <- read.table(paste0("output/",repindex,"/4_out_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 4 Biallelic
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/4)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 732){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3]
        j <- j+1
        i <- i+4
    }

    biallelic.ploidy4 <- biallelic.matrix.merged

    ## Ploidy 6
    rawdata <- read.table(paste0("output/",repindex,"/6_out_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 6 Biallelic
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/6)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 1098){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3] + biallelic.matrix[,i+4] + biallelic.matrix[,i+5]
        j <- j+1
        i <- i+6
    }

    biallelic.ploidy6 <- biallelic.matrix.merged

    ## Ploidy 8
    rawdata <- read.table(paste0("output/",repindex,"/8_out_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 8 Biallelic
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/8)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 1464){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3] + biallelic.matrix[,i+4] + biallelic.matrix[,i+5] + biallelic.matrix[,i+6] + biallelic.matrix[,i+7]
        j <- j+1
        i <- i+8
    }

    biallelic.ploidy8 <- biallelic.matrix.merged

    ## Ploidy 4 - Natural Pairing
    rawdata <- read.table(paste0("output/",repindex,"/4_out_natural_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 4 Biallelic - Natural Pairing
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/4)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 732){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3]
        j <- j+1
        i <- i+4
    }

    biallelic.ploidy4.natural <- biallelic.matrix.merged
    
    ## Ploidy 6 - Natural Pairing
    rawdata <- read.table(paste0("output/",repindex,"/6_out_natural_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 6 Biallelic - Natural Pairing
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/6)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 1098){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3] + biallelic.matrix[,i+4] + biallelic.matrix[,i+5]
        j <- j+1
        i <- i+6
    }

    biallelic.ploidy6.natural <- biallelic.matrix.merged

    ## Ploidy 8 Natural
    rawdata <- read.table(paste0("output/",repindex,"/8_out_natural_founderalleles.dat"),header=TRUE)
    rownames(rawdata) <- rawdata[,1]
    rawdata <- rawdata[,-1]
    rawdata <- t(rawdata)

    data <- data.frame(rawdata)
    data <- lapply(data,as.factor)

    ## Ploidy 8 Biallelic
    biallelic.data <- data
    biallelic.matrix <- matrix(NA,nrow=length(biallelic.data),ncol=length(biallelic.data[[1]]))
    rownames(biallelic.matrix) <- colnames(rawdata)
    colnames(biallelic.matrix) <- rownames(rawdata)
    tot.alleles <- length(levels(data[[1]]))
    for(i in 1:length(data)){
        levels(biallelic.data[[i]]) <- sample(0:1,size=tot.alleles,replace=TRUE)
        biallelic.matrix[i,] <- as.numeric(as.character(biallelic.data[[i]]))
    }


    biallelic.matrix.merged <- matrix(NA,nrow=nrow(biallelic.matrix),ncol=ncol(biallelic.matrix)/8)
    rownames(biallelic.matrix.merged) <- rownames(biallelic.matrix)
    colnames(biallelic.matrix.merged) <- indnames
    j <- i <- 1

    while(i < 1464){
        biallelic.matrix.merged[,j] <- biallelic.matrix[,i] + biallelic.matrix[,i+1] + biallelic.matrix[,i+2] + biallelic.matrix[,i+3] + biallelic.matrix[,i+4] + biallelic.matrix[,i+5] + biallelic.matrix[,i+6] + biallelic.matrix[,i+7]
        j <- j+1
        i <- i+8
    }

    biallelic.ploidy8.natural <- biallelic.matrix.merged
   
    save(biallelic.ploidy2,
         biallelic.ploidy4,
         biallelic.ploidy6,
         biallelic.ploidy8,
         biallelic.ploidy4.natural,
         biallelic.ploidy6.natural,
         biallelic.ploidy8.natural,
         file=paste0("biallelicdata/mat_",repindex,".Rdata"))

    print(repindex)
}
