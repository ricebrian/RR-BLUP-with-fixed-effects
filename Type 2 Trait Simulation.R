#This function will create simulated data for a trait with genes that following a geometric series.
# The user will input the number of QTL, the heritabilities to be explored, the additive effect,
# and the output directory

create.simluated.data <- function(){
  #Experimental code
  setwd(home.dir)
  count.filenum <- 0
  for(filenum in file.from:file.to){
    myFRG=GAPIT.Fragment(file.path=NULL,file.from=NULL, file.to=NULL,file.total=1,file.G=file.G,
                         file.Ext.G=file.Ext.G,seed=123,SNP.effect="Add",SNP.impute="Middle",
                         genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file = filenum, file.fragment=file.fragment,
                         LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)
    
    if(count.filenum == 0){
      all.FRGGs <- myFRG$G
    }else{
      all.FRGGs <- rbind(all.FRGGs, myFRG$G[-1,])
    }#End if(count.filenum == 0)
    count.filenum <- count.filenum+1
  }#end for(filenum in file.from:file.to) 
  
  hm=GAPIT.HapMap(G = all.FRGGs,SNP.effect="Add",SNP.impute="Major")
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
  
  #hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
  
  ###############################
  #temp code#
  hm.GD.without.snps.below.0.05.maf <- hm$GD
  
  genotypes <- data.frame(hm$GI[,1], rep(NA,nrow(hm$GI)),hm$GI[,2:3],rep(NA,nrow(hm$GI)),
                          t(hm.GD.without.snps.below.0.05.maf))
  
  colnames(genotypes) <- c("Snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
  
  #End experimental code
  
  #Create a working directory for the output results:
  dir.create(paste(home.dir,"/", output.dir, sep = ""))
  
  #Set the working directory
  setwd(paste(home.dir,"/", output.dir, sep = ""))
  
  #Randomly select (without replacement) k additive QTN, and assign an effect size
  seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- 130164
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.add.QTN <- sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
  Add.QTN.genotypic.information <- genotypes[vector.of.add.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Additive.QTN.number,"Add.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Add.QTN.genotypic.information, paste("Genotypic.information.for.", Additive.QTN.number,".Additive.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
  seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- 130164
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.epi.QTN <- sample(1:nrow(genotypes), (2*Epistatic.QTN.number), replace = FALSE)
  Epi.QTN.genotypic.information <- genotypes[vector.of.epi.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Epistatic.QTN.number,"Epi.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Epi.QTN.genotypic.information, paste("Genotypic.information.for.", Epistatic.QTN.number,".Epistatic.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Create a "base line" trait, which is basically just the additive effects; this is what we would see if 
  # the heritability of the simulated trait were 1
  additive.effect.trait.object <- t(Add.QTN.genotypic.information[,-c(1:5)]) #this was originally the base.line.trait.object
  
  epistatic.effect.trait.object <-t(Epi.QTN.genotypic.information[,-c(1:5)])
  #AEL Changed: - epistatic.effect.trait.object<- epistatic.effect.trait.object[,number.of.epistasis]
  
  colnames(additive.effect.trait.object) <- paste("Chr_",Add.QTN.genotypic.information[,3], "_",Add.QTN.genotypic.information[,4],sep = "")
  colnames(epistatic.effect.trait.object) <- paste("Chr_",Epi.QTN.genotypic.information[,3], "_",Epi.QTN.genotypic.information[,4],sep = "")
  
  #make base.line.trait additive.component and epistatic.component
  additive.component<- as.data.frame(matrix(0, nrow = nrow(additive.effect.trait.object), ncol = 1))
  epistatic.component<- as.data.frame(matrix(0, nrow = nrow(epistatic.effect.trait.object), ncol = 1))
  #base.line.trait <- as.data.frame(matrix(0, nrow = nrow(base.line.trait.object), ncol = 1)) 
  
  additive.component <- additive.component + (additive.effect.trait.object[,1]*(big.additive.QTN.effect))
  if(Additive.QTN.number >= 2){
    for(i in 2:Additive.QTN.number) additive.component <- additive.component + (additive.effect.trait.object[,i]*(additive.effect^(i-1)))
  }#end if(Additive.QTN.number >= 2)
  rownames(additive.component) <- rownames(additive.effect.trait.object)
  colnames(additive.component) <- "Additive.effect"
  additive.genetic.variance <- var(additive.component)
  
  last.number.of.this.loop <- Epistatic.QTN.number - 1
  for(i in 0:last.number.of.this.loop) epistatic.component <- epistatic.component + ((epistatic.effect.trait.object[,((2*i)+1)]*epistatic.effect.trait.object[,((2*i)+2)])*(epistatic.effect^(i+1)))
  rownames(epistatic.component) <- rownames(epistatic.effect.trait.object)
  colnames(epistatic.component) <- "Epistatic.effect"
  epistatic.genetic.variance<- var(epistatic.component)
  
  #Set the row names of the base.line.trait object to the new names
  base.line.trait <- additive.component+epistatic.component
  #base.line.trait.with.new.taxa <- merge(base.line.trait, taxa.name.converter, by.x = "row.names", 
  #                                       by.y = "Old_Taxa_ID")
  
  #the.new.taxa.ids <- as.character(base.line.trait.with.new.taxa[,2])
  #base.line.trait <- as.matrix(base.line.trait.with.new.taxa[,2], nrow = nrow(base.line.trait.with.new.taxa))
  # rownames(base.line.trait) <- as.character(base.line.trait[,3])
  
  
  #For loop through the vector of heritabilities
  for(i in heritabilities.vector){
    #If heritability is zero
    if(i == 0){ 
      #Simulate m replicates of n N(0,b) random variables, where b = additive.genetic.variance
      the.seed.number.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = 1)
        if(j == 1){
          simulated.data <- the.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
      }
      
      #Format the output file for the simulated phenotypes
      simulated.data <- cbind(rownames(base.line.trait),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      colnames(simulated.data)[2] <- "the.normal.random.variables.1"
      
      
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
      
    }else{
      #Calcualte V_e, the residual variance
      #residual.variance <- (additive.genetic.variance*(1-i))/i
      residual.variance <-( ((additive.genetic.variance+epistatic.genetic.variance)/i) - additive.genetic.variance - epistatic.genetic.variance)
      #Use V_e to generate n replicates of N(0, Ve) random variables
      the.seed.number.vector <- NULL
      col.name.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = sqrt(residual.variance))
        the.base.line.trait.plus.normal.random.variables <- base.line.trait+the.normal.random.variables
        if(j == 1){
          simulated.data <- the.base.line.trait.plus.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.base.line.trait.plus.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
        col.name.vector <- c(col.name.vector, paste("Heritability_",i, "_Rep_", j, sep = ""))
      }
      
      colnames(simulated.data)  <- col.name.vector
      
      #Format the output file for the simulated phenotypes
      simulated.data <- cbind(rownames(base.line.trait),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
    }   
  }#End for(i in heritabilities.vector)
  
}#end "create.simluated.data()"





###########################################################################################
###########################################################################################
###########################################################################################
#setwd("/Users/adminuser/Box Sync/Lipka_Mainzer_Chen_Epistasis_Shared_Folder/Simulation_Study")
setwd("/Users/adminuser/Box Sync/Rice_Lipka_Shared_Folder/Thesis/Code_for_Simulating_Traits/From 282")
home.dir <- getwd()
#dir.of.GBS.SNPs <- "/Users/adminuser/Desktop/Work/Tocos_NAM_2009_2010/Joint_Linkage_Analysis/GBS_SNPs/"


#Read in the 1,106 markers that are genotyped on the NAM familes scored for kernel color in Chandler et al. 2013

#setwd(dir.of.GBS.SNPs)
#genotypes <- read.table("4K_SNPsmdp_genotype_test1.hmp.txt", head = TRUE)
#setwd(home.dir)
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
#source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.R")
source("http://zzlab.net/GAPIT/emma.txt")
#library(rrBLUP)

j<-c(25)
w<-c(0.9)
s<-c(0.9)
for (i in j){
  for (t in s){
    for (k in w){
      ###############
      #User input
      setwd("C:/Users/brice6/Desktop/Thesis WorkFile/Sorghum/Sorghum/Type 2 Simulations")
      home.dir <- getwd()
      #Number of additive QTN (k)
      Additive.QTN.number <- i
      
      
      #Number of epistatic QTN (m)
      Epistatic.QTN.number <- 2
      
      #Vector of heritabilities to investigate
      heritabilities.vector <- k
      
      #Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
      big.additive.QTN.effect <- t
      additive.effect <- 0.1
      
      
      #Size of the epistatic effect of the largest QTL (must be (-1,1) but preferably (0,1)|
      epistatic.effect <- 0
      
      #Number of replicates of simulated phenotypes for each heritability (m)
      replicates <- 50
      
      
      file.G="GBS_Markers_Sorghum_116128_n=320_Inds.txt" 
      file.Ext.G = "hmp.txt"
      
      
      file.from = 1
      file.to = 1
      
      file.fragment = 1000000
      
      
      #Output directory
      output.dir <- paste(Additive.QTN.number,"_Add_QTN",Epistatic.QTN.number,"_Epi_QTN_h.2_",
                          heritabilities.vector,"_add.eff_", big.additive.QTN.effect,"_epis.eff_", epistatic.effect,"_reps_", replicates, sep = "")
      
      ################
      #Create the simulated data
      create.simluated.data()
      
      
      print(i)
      print(t)
      print(k)
    }
  }
}






