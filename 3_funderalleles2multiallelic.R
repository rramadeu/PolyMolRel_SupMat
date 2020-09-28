## Companion script for the manuscript:
## "Estimation of Molecular Pairwise Relatedness in Autopolyploid Crops"
## by Rodrigo Amadeu, Leticia Lara, Patricio Munoz, Antonio Garcia
## Issue/Year: TBA

## This R script:
## reads PedigreeSim output and sample biallelic and multiallelic markers to all loci. It saves the data as a .txt following PolyRelatedness V1.8 format

system("mkdir multiallelicdata")
system("cd multiallelicdata; mkdir 2ploidy 4ploidy 6ploidy 8ploidy 4ploidy_natural 6ploidy_natural 8ploidy_natural")

indnames <- read.table("PedigreeSimInput/0.ped",header=TRUE,colClasses="character")[,1]
alleles <- c(2:5,10,15,20)
loci <- c(5,10,15,20,50,100,500,1000)
reptotal <- 3

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/2ploidy/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/2ploidy/", loci[i],"loci/",alleles[j],"alleles"))
}

## Diploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/2_out_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/2)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 366){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1])
            j <- j+1
            i <- i+2
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/2ploidy/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}


for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/4ploidy/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/4ploidy/", loci[i],"loci/",alleles[j],"alleles"))
}

## Tetraploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/4_out_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/4)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 732){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3])
            j <- j+1
            i <- i+4
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/4ploidy/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/6ploidy/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/6ploidy/", loci[i],"loci/",alleles[j],"alleles"))
}

## Hexaploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/6_out_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/6)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 1098){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3],
                                              multi.matrix[,i+4],multi.matrix[i+5])
            j <- j+1
            i <- i+6
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/6ploidy/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/8ploidy/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/8ploidy/", loci[i],"loci/",alleles[j],"alleles"))
}

## Octaploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/8_out_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/8)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 1464){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3],
                                              multi.matrix[,i+4],multi.matrix[i+5],
                                              multi.matrix[,i+6],multi.matrix[i+7])
            j <- j+1
            i <- i+8
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/8ploidy/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}

## Tetraploid Natural

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/4ploidy_natural/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/4ploidy_natural/", loci[i],"loci/",alleles[j],"alleles"))
}


## Tetraploid loop Natural
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/4_out_natural_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/4)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 732){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3])
            j <- j+1
            i <- i+4
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/4ploidy_natural/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/6ploidy_natural/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/6ploidy_natural/", loci[i],"loci/",alleles[j],"alleles"))
}

## Hexaploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/6_out_natural_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/6)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 1098){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3],
                                              multi.matrix[,i+4],multi.matrix[i+5])
            j <- j+1
            i <- i+6
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/6ploidy_natural/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}

## Octaploid Natural

for(i in 1:length(loci)){
    system(paste0("mkdir multiallelicdata/8ploidy_natural/",loci[i],"loci"))
    for(j in 1:length(alleles))
        system(paste0("mkdir multiallelicdata/8ploidy_natural/", loci[i],"loci/",alleles[j],"alleles"))
}

## Octaploid loop
for(repindex in 1:reptotal){
    rawdata.rep <- read.table(paste0("output/",repindex,"/8_out_natural_founderalleles.dat"),header=TRUE)
    for(allelesindex in 1:length(alleles)){
        rawdata <- rawdata.rep
        rownames(rawdata) <- rawdata[,1]
        rawdata <- rawdata[,-1]
        rawdata <- t(rawdata)
        data <- data.frame(rawdata)
        data <- lapply(data,as.factor)
        
        multi.data <- data
        multi.matrix <- matrix(NA,nrow=length(multi.data),ncol=length(multi.data[[1]]))
        rownames(multi.matrix) <- colnames(rawdata)
        colnames(multi.matrix) <- rownames(rawdata)
        tot.alleles <- length(levels(data[[1]]))
        for(i in 1:length(data)){
            levels(multi.data[[i]]) <- sample(0:(alleles[allelesindex]-1),size=tot.alleles,replace=TRUE)
            multi.matrix[i,] <- as.numeric(as.character(multi.data[[i]]))
        }
        
        multi.matrix <- multi.matrix+100
        
        multi.matrix.merged <- matrix(NA,nrow=nrow(multi.matrix),ncol=ncol(multi.matrix)/8)
        rownames(multi.matrix.merged) <- rownames(multi.matrix)
        colnames(multi.matrix.merged) <- indnames
        j <- i <- 1
        while(i < 1464){
            multi.matrix.merged[,j] <- paste0(multi.matrix[,i],multi.matrix[,i+1],
                                              multi.matrix[,i+2],multi.matrix[i+3],
                                              multi.matrix[,i+4],multi.matrix[i+5],
                                              multi.matrix[,i+6],multi.matrix[i+7])
            j <- j+1
            i <- i+8
        }
        
        multi.matrix.merged <- t(multi.matrix.merged)
        multi.matrix.merged <- cbind("AA",multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Pop"
        multi.matrix.merged <- cbind(rownames(multi.matrix.merged),multi.matrix.merged)
        colnames(multi.matrix.merged)[1] <- "Ind"

        for(lociindex in 1:length(loci)){
            out.matrix <- multi.matrix.merged[,c(1,2,sample(1:10000,loci[lociindex])+2)]       
            filename <- paste0(
                "multiallelicdata/8ploidy_natural/",
                loci[lociindex],"loci/",
                alleles[allelesindex],"alleles/",
                repindex,".txt")
            
            ## Building the file in PolyReladtness format
            write.table(out.matrix,"geno.txt",quote=FALSE,row.names=FALSE,sep="\t")
            system(paste0("echo '//configuration' > ",filename))
            system(paste0("echo '//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)' >> ",filename))
            system(paste0("echo '3    8	0	999	4' >> ",filename))
            system(paste0("echo '//genotypes' >> ",filename))
            system(paste0("cat geno.txt >> ",filename))
            system(paste0("echo '//end of file' >> ",filename))
            print(filename)
        }
    }
}
