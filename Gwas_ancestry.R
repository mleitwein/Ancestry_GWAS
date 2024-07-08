### GWAS avec Gemma ###

## AquaExcel ##

# mars 2022

#libraries
library(data.table)
library("lme4")
library(MuMIn)
    
        
########## GEMMA Haploytpes ################
options(scipen = 999) 
# recode genotypes in ancestry info: 012 : 0=00 (atl), 1=01,  2=11 (est)
01_utility_script/loter_012_info.R

pheno<- fread("GWAS/phenotype_AquaExcel.txt", header=T, sep="\t") 

#### GWAS with Fat musculat content as phenotype

mod7<-lm( fat ~ poids_vif + SEXE, na.action = na.exclude ,data=pheno)  
  
  summary(mod7)
  head(resid(mod7))
  res<-pheno[,2]
  
  res<-cbind(res, as.data.frame(resid(mod7)))
  write.table(res, "gras_seul.txt", quote = FALSE, col.names=T,row.names=FALSE, sep="\t")
  
 		##File for Gemma
		pheno<- fread("GWAS/phenotype_AquaExcel.txt", header=T, sep="\t") 
		  
		#haplotype on garde que les individus qui sont dans le phenotype      
		  dt<-fread("GWAS/genotyp_imp_haplotype.txt", header=F)%>%
			dplyr::filter(V1%in%pheno$CODE)
		  
		  #phenotype on garde que les individus present dans le genotype et on les ordonnes dans le bon sens
		  pheno <- pheno %>%
			dplyr::filter(CODE%in%dt$V1) %>%
			dplyr::arrange(match(CODE,dt$V1))
		  
		  FirstCol <- as.data.frame(dt[,1])
		  #on reorganise le format genotype      
		  
		  dt <- t(dt) #transposition
		  dt <- dt[-nrow(dt),]
		  dt <- dt[-1,]
		  #SNP file info        
		  map <- read.table("../snp_info.txt",h=T, sep="\t",stringsAsFactors = F)
		  map <- map[-nrow(map),]
		  #correspondance avec les alleles         
		  dt2 <- read.table("../SNP_allele.txt",h=F, sep=",",stringsAsFactors = F)
		  dt2$V2<-"T"
		  dt2$V3<-"C"
		  dt <- cbind(dt2[,1:3],dt)
		  #verification 
		  dt[1:10,1:10]
		  dim(dt)
		  dim(pheno)
		  dim(map)
		  # save files         
		  write.table(dt,"gemma_genotypes_Haplo_7.txt", col.names = F, row.names = F, sep=", ", quote=F)
		  write.table(map[,c("SNPID","BPPos","Chr")],"gemma_map_Haplo_7.txt", col.names = F, row.names = F, sep=", ", quote=F)
		  write.table(pheno[,!colnames(pheno)%in%"V1"],"gemma_gras_seul_7.txt", col.names = F, row.names = F, sep="\t", quote=F)
		  
		  
		### Run Gemma sur cluster
				

		## Estimate Relatedness Matrix 
		$GEM_HOME/gemma \
		  -g  gemma_genotypes_Haplo_7.txt \
		  -p  gemma_gras_seul_7.txt \
		  -gk \
		  -o kindship_matrix
		  
		  
		### Association Tests with Multivariate Linear Mixed Models
		$GEM_HOME/gemma \
		  -g  gemma_genotypes_Haplo_7.txt \
		  -p  gemma_gras_seul_7.txt \
		  -n 1 \
		  -a gemma_map_Haplo_7.txt \
		  -maf 0.1 \
		  -lmm  4 \
		  -k output/kindship_matrix.cXX.txt \
		  -o gwas_univariate_gras_Haplo_mod7 
 
 
 
 
#### GWAS with Atlantic ancestry as trait 
  
  ancestry<-fread("GWAS/ancestry_by_ind.txt") 
  names(ancestry)[1]<-"CODE"
  pheno2<-merge(ancestry, pheno, by="CODE")
  
  mod9<-lm(ATL ~ poids_vif + SEXE , na.action = na.exclude ,data=pheno2)  
  
  res<-pheno[,2]
  
 
  res<-cbind(res, as.data.frame(resid(mod9)))
  write.table(res, "Ancestry_mod9.txt", quote = FALSE, col.names=T,row.names=FALSE, sep="\t")  
  write.table(pheno2[,1:2], "Ancestry_9.txt", quote = FALSE, col.names=T,row.names=FALSE, sep="\t")  
  
  
		##File for Gemma
		#phenotype    
		  pheno <- read.table("Ancestry_9.txt",h=T,sep="\t", stringsAsFactors = F)
		  
		#haplotype on garde que les individus qui sont dans le phenotype      
		  dt<-fread("../genotyp_imp_haplotype.txt", header=F)%>%
			dplyr::filter(V1%in%pheno$CODE)
		  
		  #phenotype on garde que les individus present dans le genotype et on les ordonnes dans le bon sens
		  pheno <- pheno %>%
			dplyr::filter(CODE%in%dt$V1) %>%
			dplyr::arrange(match(CODE,dt$V1))
		  
		  FirstCol <- as.data.frame(dt[,1])
		  #on reorganise le format genotype      
		  
		  dt <- t(dt) #transposition
		  dt <- dt[-nrow(dt),]
		  dt <- dt[-1,]
		  #SNP file info        
		  map <- read.table("../snp_info.txt",h=T, sep="\t",stringsAsFactors = F)
		  map <- map[-nrow(map),]
		  #correspondance avec les alleles         
		  dt2 <- read.table("../SNP_allele.txt",h=F, sep=",",stringsAsFactors = F)
		  dt2$V2<-"T"
		  dt2$V3<-"C"
		  dt <- cbind(dt2[,1:3],dt)
		  #verification 
		  dt[1:10,1:10]
		  dim(dt)
		  dim(pheno)
		  dim(map)
		  # save files         
		  write.table(dt,"gemma_genotypes_Haplo_9.txt", col.names = F, row.names = F, sep=", ", quote=F)
		  write.table(map[,c("SNPID","BPPos","Chr")],"gemma_map_Haplo_9.txt", col.names = F, row.names = F, sep=", ", quote=F)
		  write.table(pheno[,!colnames(pheno)%in%"V1"],"gemma_ATL_9.txt", col.names = F, row.names = F, sep="\t", quote=F)
		  
		  
		### Run Gemma sur cluster
				

		## Estimate Relatedness Matrix from Genotypes
		$GEM_HOME/gemma \
		  -g  gemma_genotypes_Haplo_9.txt \
		  -p  gemma_ATL_9.txt \
		  -gk \
		  -o kindship_matrix
		  
		  
		### Association Tests with Multivariate Linear Mixed Models
		$GEM_HOME/gemma \
		  -g  gemma_genotypes_Haplo_9.txt.txt \
		  -p  gemma_ATL_9.txt \
		  -n 1 \
		  -a gemma_map_Haplo_9.txt \
		  -maf 0.1 \
		  -lmm  4 \
		  -k output/kindship_matrix.cXX.txt \
		  -o gwas_univariate_gras_Haplo_mod9  