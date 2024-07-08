################################################################
#### Creation animaux sans genotypes pour FImpute
#### R version 3.6.2
################################################################

# install.packages("tidyverse")

library(tidyverse)

### Nombre de colonnes = nombre de SNP sur ta puce
nb_snp <- 47681

### Nombre d'individus avec un genotype manquant
nb_ind <- 6

### Creation d'une amtrice de "5"
dd <- matrix(5, ncol=nb_snp, nrow=nb_ind)
dim(dd)
dd[1:6,1:10]

### Identifiant de tes individus a imputer
dd_anim <- data.frame(Ind=paste0("*",1:5), Chip=1,stringsAsFactors = F)
dd_anim2 <- data.frame(Ind=paste0("#",1), Chip=1,stringsAsFactors = F)
dd_anim3 <- rbind(dd_anim,dd_anim2)
dd_anim3

### Accoler les infos ID avec les genotypes
dd <- cbind(dd_anim3,dd)

### Accoler les genotypes les uns aux autres
DD_final <- tidyr::unite(as.data.frame(dd,stringsAsFactors=F),genotype_calls,3:ncol(dd),sep="",remove=T)
write.table(DD_final, paste0('05_Fimpute/missing_geno.txt'), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
