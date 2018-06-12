#Set the working directory with source files, kinship and genotypes
setwd("C:/Users/brice6/Desktop/Thesis WorkFile/Simulations with constant seed number/Maize")
home.dir <- getwd()
#Read in the GAPIT code and all prerequisite files
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
library(rrBLUP)
library('MASS')
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/previous/gapit_functions20160408.txt")
source("http://zzlab.net/GAPIT/emma.txt")


#parameters for folder with simulated traits
hert<-0.5
a<-0.1
n.QTN<-10 #how many QTn did we simulate?
wd<-(paste("C:/Users/brice6/Desktop/Thesis WorkFile/Simulations with constant seed number/Maize/",n.QTN,"_Add_QTN2_Epi_QTN_h.2_",hert,"_add.eff_",a,"_epis.eff_0_reps_50",sep=""))
number.of.folds <- 5
r=c(0,1,2,3,5,10,25,50,100)#Pick up the top "m" (m is user inputted) SNPs with lowest SNPs

Pval_threshold=0.05#Set a P-value threshold
filter.by.Pvalues <- FALSE

seed.number <- -84134 #sample(-1000000:1000000,1) 

#Set Files names for genotype files
file.G="SNP55K_maize282_AGPv2_20100513_" 
file.Ext.G = "hmp.txt"
file.from=1
file.to=1
SNP.fraction=1
group.from = 201 
group.to = 201
SNP.fraction = 1
file.fragment = 150000 #needs to be greater than number of markers for myFRG to work
PCA.total=3

setwd(wd)
myY <- read.table(paste("Simulated.Data.50.Reps.Herit.",hert,".txt",sep="" ),head = TRUE)
myY <- myY[,1:2]

setwd("C:/Users/brice6/Desktop/Thesis WorkFile/Simulations with constant seed number/Maize")
#setwd("C:/Users/brice6/Desktop/Current Simulations")
#setwd("E:/Thesis/Simulating_Traits_Nam/282 55k snp simulation")
home.dir <- getwd()
#Read in the kinship matrix
myKI <- read.csv("Comprehensive_K.csv", head = FALSE)
#names<-colnames(myKI)
#myKI<-cbind(myKI[,1],myKI[,2:ncol(myKI)]-1)
#colnames(myKI)<-names
#rownames(myKI) <-myKI[,1]
#myKI<-myKI[,-1]
myPC<-read.csv("GAPIT.PCA.csv",head=TRUE)
myPC<-myPC[,3]
#In the near future, we want to conduct a preliminary
# run of GAPIT so that we can find an optimal comression level (using all of the data).
# It is a good idea to do this becuase i.) GAPIT will not be searching for an 
# optimal compression level at each fold, and ii.) thus this will substantially
# speed up computational time.

#####run calculate the kinship matrix that is going to be used in rrBLUP.

#This will enable us to read in multiple genotype files
#Needed if no kinship matrix
#myFRG=GAPIT.Fragment(file.path=NULL,file.from=file.G, file.to=file.to,file.total=NULL,file.G=file.G,
#                     file.Ext.G=file.Ext.G,seed=123,SNP.fraction=SNP.fraction,SNP.effect="Add",SNP.impute="Middle",
#                     genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=file.fragment,
#                     LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)

myFRG=GAPIT.Fragment(file.path=NULL,file.from=file.G, file.to=file.to,file.total=NULL,file.G=file.G,
                     file.Ext.G=file.Ext.G,seed=123,SNP.effect="Add",SNP.impute="Middle",
                     genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=file.fragment,
                     LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)
#

hm=GAPIT.HapMap(G = myFRG$G,SNP.effect="Add",SNP.impute="Major")

#####################################
#Obtain the mafs of all SNPs

#Total number of lines
ns <- nrow(hm$GD)

#Sum of the allele scores for each SNP
ss <- apply(hm$GD, 2, sum)

#Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))

#Copy the minor allele frequencies for all SNPs
maf <- apply(maf.matrix, 2, min)

#Find out which SNPs have MAF < 0.05
snps.below.0.05.maf <- which(maf < 0.05)

# Remove these SNPs from hm$GD

#need to add an if else statement for when no markers are below 0.05
hm.GD.without.snps.below.0.05.maf <- hm$GD#[,-snps.below.0.05.maf]

CV=myY[,1:2]
CV[,2]=1
###########add code to name correct trait
colnames(CV)=c("SampleID","traittrans6")


GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)

qc=GAPIT.QC(Y = myY, GT = hm$GT, CV = CV, GK = GK)

y <- as.matrix(qc$Y[,-1])

G <- as.numeric(qc$GK[,-1])

G <- matrix(G, nrow(y), ncol(qc$GK[,-1]))

G <- G - 1

cv <- (as.matrix(qc$CV[,-1]))

#Calculate the kinship matrix in rrBLUP #why is this a new step?-Brian
A1 <- A.mat(G,shrink=TRUE)


#Use the objects from hm (calcualted in Line 70) to obtain the genotype files for GAPIT
myGD <- data.frame(hm$GT, hm$GD)
colnames(myGD) <- c("taxa", as.character(hm$GI[,1]))
myGM <- hm$GI

#Save all of the above work into an object
#save.image("Workspace_20170815.Rdata")

sample.size <- length(y)
set.seed(seed.number)
sequence.sample <- rep(1:sample.size)
random.sample <- sample(1:sample.size, replace = FALSE)
increment <- ceiling(length(random.sample)/number.of.folds) 


#have a "for" loop, start it at 0, and end it at 4
#I am setting up "k" to denote the nubmer of folds - 1. This is done
# so that the for loop will work correctly.
this.genotype.file <- read.delim(paste(file.G,file.from,".",file.Ext.G,sep = ""), head = FALSE)

r.gy.m<-matrix(NA,nrow=1,ncol=(number.of.folds+4))
setwd(wd)


#Read in a phenotype
for (p in 2:51){
  pheno <- read.table(paste("Simulated.Data.50.Reps.Herit.",hert,".txt",sep="" ),head = TRUE)
  myY <- pheno[,c(1,p)]
  k <- number.of.folds - 1
  for (i in 0:k){ 
    #pick out of he training set (t.s.) (k-1) folds
    
    print(paste("----------------Starting ", (i+1), " out of ", number.of.folds, " folds-----------------", sep = ""))
    
    pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
    train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
    
    #Perform GWAS on T.S.
    myY.train <- myY[train,]
    
    #Step 2: Run GAPIT
    myGAPIT <- GAPIT(
      Y=myY.train, 
      #PCA.total=PCA.total,
      #CV=myPC,
      GD = myGD,
      GM = myGM,
      K = myKI,
      SNP.fraction=SNP.fraction,
      group.from = group.from,
      group.to = group.to,
      file.fragment = file.fragment,
      
      Geno.View.output=FALSE
    )
    
    #Read in the GWAS output file
    
    #Note to Brian and Alex: we need to get this to work for multiple traits
    Gwas.output<-read.csv(paste("GAPIT..",colnames(myY.train)[2],".GWAS.Results.csv", sep = ""))
    
    #Write table to be used for making manhattan plot of just fixed markers
    write.table(Gwas.output, paste("GWAS.output.with.causative.markers.fold=",i+1,".txt", sep = ""),quote = FALSE, row.names = FALSE,col.names = TRUE)
    #I won't save just m markers. When I make plots I can set m then rather than have a b
    #of files now (this way I only need 5 rather than 5 for each value of m)
    
    
    #This is where we need to remove casual mutation from GWAS.output #This will create more realistic simualted cirrcumstances
    causative.mutations<-read.table(paste("Genotypic.information.for.",n.QTN,".Additive.QTN.txt",sep = ""),header=T)
    names<-causative.mutations["Snp"] #names of causative snps
    Q<-c()
    for (l in 1:nrow(names)){
      Q[l]<-as.character(names[l,1])}
    #want to save to me() the row number in Gwas.output that matches causative snps names in Q
    me<-c() 
    for (h in 1:length(Q)){
      w<-Q[h]
      me[h]<-which(Gwas.output$SNP==w)
    }
    Gwas.output<-Gwas.output[-c(me),]
    #now extract the first m markers; which will be the "top m" most significant marker
    write.table(Gwas.output, paste("GWAS.output.with.markers.removed.fold=",i+1,".txt",sep =""),quote = FALSE, row.names = FALSE,col.names = TRUE)
  }  ####### end loop for GWAS######
  
  
  
  #start new loop here for rrBLUP with fixed effects
  #adding loop to conduct for each level of m in vector r
  
  k <- number.of.folds-1
  for (m in r){
    r.gy<-NULL
    for (i in 0:k){
      
      pred <- random.sample[((increment*i)+1):min(((increment*i)+increment) , sample.size)]
      train <- random.sample[-(((increment*i)+1):min(((increment*i)+increment) , sample.size))] 
      if(m != 0){
        #read back in GWAS output and use it in in new loop for K-Fold cv
        #read in file with only m rows
        Gwas.output<-read.table(paste("GWAS.output.with.markers.removed.fold=",i+1,".txt",sep="" ),header=T,nrows = m)
        #Filter out SNPs that have an FDR-adjusted P-value greater than 0.05
        if(filter.by.Pvalues) Gwas.output <- Gwas.output[which(Gwas.output[,ncol(Gwas.output)] <= Pval_threshold),]
        #Extract the SNP name; also extract the chromosome and bp information
        
        print(paste("------- Now obtaining SNP information of the peak ", m, " SNPs from fold ", (i+1), " -----------",sep = ""))
        #Read in the chromosome file that has the SNP of interest
        if(file.from != file.to){
          count.chr <- 0
          for(j in unique(Gwas.output[,2])){#{2
            this.genotype.file <- read.delim(paste(file.G,j,".",file.Ext.G,sep = ""), head = FALSE)  
            the.SNP.ids.on.this.chr <- Gwas.output[which(Gwas.output[,2] == j),1]
            
            #Pick out the SNPs that we want
            these.specific.genotypes <- this.genotype.file[which(this.genotype.file[,1] %in% the.SNP.ids.on.this.chr),]
            
            #This knocks out two birds with one stone. 
            # stone 1: if there is only one chromosome with peak SNPs of interest, it will create an output file
            #          called "the specific genotypes"
            # Stone 2: if there are >1 chromosome with peak SNPs of interest, it will initialize an output file containing
            #          associated SNPs from multiple chromosomes
            if(count.chr == 0) the.specific.genotypes <- rbind(this.genotype.file[1,],these.specific.genotypes)
            
            if(length(unique(Gwas.output[,2])) > 1){
              if(count.chr > 1){#{4
                the.specific.genotypes <- rbind(the.specific.genotypes, these.specific.genotypes)
              }
              count.chr <- count.chr + 1
            }
            
            #End product: "the.specific.genotypes" - one HapMap-formatted genotype file,
            # where there are potentially peak SNPs from multiple chromsomes.
            
          } #for(j in unique(Gwas.output[,2]))
        }else{
          #lets put the following line earlier
          #this.genotype.file <- read.delim(paste(file.G,file.from,".",file.Ext.G,sep = ""), head = FALSE)
          the.SNP.ids <- Gwas.output[,1]
          
          #Pick out the SNPs that we want
          the.specific.genotypes <- this.genotype.file[c(1,which(this.genotype.file[,1] %in% the.SNP.ids)),]
        }#End else
        
        #the.specific.genotypes <- rbind(the.specific.genotypes, the.specific.genotypes[1,])
        
        #get into hapmap format
        hm=GAPIT.HapMap(G = the.specific.genotypes,SNP.effect="Add",SNP.impute="Major")
        
        #What we need is hm$GD - these are the genotypes of the peak SNPs in numeric format
        this.cv <- data.frame(hm$GT, hm$GD)
        
        #Assign appropriate column names to this.Cv
        colnames(this.cv) <- c("SampleID", as.character(hm$GI[,1]))
        
        
        this.qc=GAPIT.QC(Y = myY, GT = hm$GT, CV = this.cv, GK = GK)
        
        cv.for.rrBLUP <- (as.matrix(this.qc$CV[,-1]))
        
      }else{
        cv.for.rrBLUP <- as.vector(rep(1, length(y)))
      }#end if(m != 0)
      #Run them as fixed effect covariates fitting in rrBLUP model in T.S
      #######Below code taht will do the CV in one fold
      
      print(paste("-------Now fitting the RR-BLUP model for fold ", (i+1), " -----------", sep = ""))
      
      yNA <- y
      yNA <- as.vector(yNA)
      yNA[pred] <- NA
      
      data1 <- data.frame(y=yNA,gid=1:length(y), cv = cv.for.rrBLUP)
      the.cv.names <- NULL
      for(j in 1:ncol(cv)) the.cv.names <- c(the.cv.names, paste("CV_",j,sep = ""))
      
      colnames(data1) <- c("y","gid", the.cv.names)
      
      rownames(A1) <- 1:nrow(A1) #A1 is created on line 114
      ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y", covariate = the.cv.names)
      #Measure correclation between OBS and Pred in validation set (V.S.)
      r.gy <- c(r.gy, cor(ans1$g[pred], y[pred]) )
      
    } 
    r.gy <- c(r.gy, mean(r.gy), sd(r.gy),m,p)
    #calcualte the average and std. deviation
    r.gy.output <- t(as.matrix(r.gy))
    
    r.gy.m<-rbind(r.gy.m,r.gy.output) 
    
    colnames(r.gy.m)<-c("fold 1","fold 2","fold 3","fold 4","fold 5","mean","sd","m","rep")
    write.table(r.gy.m,paste( number.of.folds,"folds_CV_results_",n.QTN,"QTN.txt"),quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  }
}#for (p in 1:50){



