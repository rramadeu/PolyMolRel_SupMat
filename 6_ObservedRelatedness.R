## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Sep-2020

## This script:
## computes the observed relationship matrix on the original data following delta computation (Equation 2 from the paper)

system("mkdir ObservedRelatedness")

r.compute <- function(x, ploidy){
    if(ploidy==2){
        all.deltas <- x[c(1,2)] %in% x[c(3,4)]
        r <- (all.deltas == 1)/2 + (all.deltas == 2)
        return(r)
    }

    if(ploidy==4){
        all.deltas <- x[c(1,2,3,4)] %in% x[c(5,6,7,8)]
        r <- (all.deltas == 1)/4 + (all.deltas == 2)/2 +
            (all.deltas == 3)*3/4 + (all.deltas == 4)
        return(r)
    }

    if(ploidy==6){
        all.deltas <- x[c(1,2,3,4,5,6)] %in% x[c(7,8,9,10,11,12)]
        r <- (all.deltas == 1)/6 + (all.deltas == 2)*2/6 +
            (all.deltas == 3)*3/6 + (all.deltas == 4)*4/6 +
            (all.deltas == 5)*5/6 + (all.deltas == 6)
        return(r)
    }

    if(ploidy==8){
        all.deltas <- x[c(1,2,3,4,5,6,7,8)] %in% x[c(9,10,11,12,13,14,15,16)]
        r <- (all.deltas == 1)/8 + (all.deltas == 2)*2/8 +
            (all.deltas == 3)*3/8 + (all.deltas == 4)*4/8 +
            (all.deltas == 5)*5/8 + (all.deltas == 6)*6/8 +
            (all.deltas == 7)*7/8 + (all.deltas == 8)
        return(r)
    }
    
}

r.tot <- function(x, ploidy){
    return(mean(colSums(apply(x,1,r.compute,ploidy))))
}



inds <- read.table("PedigreeSimInput/0.ped",colClasses="character",header=TRUE)[,1]
reps <- 3

## ploidy 2
mat.array <- array(NA,dim=c(183,183,reps))
index <- combn(1:183,2)

for(j in 1:reps){
    print(paste("Ploidy 2 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/2_out_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*2-1,index[1,i]*2)
        ind2 <- c(index[2,i]*2-1,index[2,i]*2)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=2)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/r2ploidy.Rdata")

## ploidy 4
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:reps){
    print(paste("Ploidy 4 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/4_out_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*4-3,index[1,i]*4-2,index[1,i]*4-1,index[1,i]*4)
        ind2 <- c(index[2,i]*4-3,index[2,i]*4-2,index[2,i]*4-1,index[2,i]*4)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=4)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/r4ploidy.Rdata")

## ploidy 6
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:reps){
    print(paste("Ploidy 6 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/6_out_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*6-5,index[1,i]*6-4,index[1,i]*6-3,index[1,i]*6-2,index[1,i]*6-1,index[1,i]*6)
        ind2 <- c(index[2,i]*6-5,index[2,i]*6-4,index[2,i]*6-3,index[2,i]*6-2,index[2,i]*6-1,index[2,i]*6)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=6)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/r6ploidy.Rdata")

## ploidy 8
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:100){
    print(paste("Ploidy 8 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/8_out_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*8-7,index[1,i]*8-6,index[1,i]*8-5,index[1,i]*8-4,index[1,i]*8-3,index[1,i]*8-2,index[1,i]*8-1,index[1,i]*8)
        ind2 <- c(index[2,i]*8-7,index[2,i]*8-6,index[2,i]*8-5,index[2,i]*8-4,index[2,i]*8-3,index[2,i]*8-2,index[2,i]*8-1,index[2,i]*8)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=8)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/r8ploidy.Rdata")



## ploidy 4
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:reps){
    print(paste("Ploidy 4 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/4_out_natural_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*4-3,index[1,i]*4-2,index[1,i]*4-1,index[1,i]*4)
        ind2 <- c(index[2,i]*4-3,index[2,i]*4-2,index[2,i]*4-1,index[2,i]*4)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=4)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/rnatural4ploidy.Rdata")

## ploidy 6
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:reps){
    print(paste("Ploidy 6 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/6_out_natural_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*6-5,index[1,i]*6-4,index[1,i]*6-3,index[1,i]*6-2,index[1,i]*6-1,index[1,i]*6)
        ind2 <- c(index[2,i]*6-5,index[2,i]*6-4,index[2,i]*6-3,index[2,i]*6-2,index[2,i]*6-1,index[2,i]*6)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=6)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/rnatural6ploidy.Rdata")

## ploidy 8
mat.array <- array(NA,dim=c(183,183,100))
index <- combn(1:183,2)
  
for(j in 1:reps){
    print(paste("Ploidy 8 - it",j))
    df <- as.matrix(read.table(
        paste0("output/",
               j,"/8_out_natural_founderalleles.dat"),header=TRUE))
    rownames(df) <- df[,1]
    df <- df[,-1]
    class(df) <- "numeric"
    
    r.mat <- matrix(NA,183,183)
    diag(r.mat) <- 1

    for(i in 1:ncol(index)){
        ind1 <- c(index[1,i]*8-7,index[1,i]*8-6,index[1,i]*8-5,index[1,i]*8-4,index[1,i]*8-3,index[1,i]*8-2,index[1,i]*8-1,index[1,i]*8)
        ind2 <- c(index[2,i]*8-7,index[2,i]*8-6,index[2,i]*8-5,index[2,i]*8-4,index[2,i]*8-3,index[2,i]*8-2,index[2,i]*8-1,index[2,i]*8)
        r.mat[index[1,i],index[2,i]] <- r.mat[index[2,i],index[1,i]] <- r.tot(df[,c(ind1,ind2)],ploidy=8)
    }
    mat.array[,,j] <- r.mat
}

save(mat.array,file="ObservedRelatedness/rnatural8ploidy.Rdata")
