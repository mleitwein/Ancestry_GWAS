# Loter encodage 0 1 2
require(dplyr)
require(data.table)


setwd("~/PostDoc_Ifremer/02_AquaExcel/06_Loter/03_result/01_Loter_res")

# load data
for (i in 1:25){
 
ancestry_chr <- read.csv(paste0("chr", i,"/ancestry_chr", i, ".txt"), sep = " ", header = FALSE)

# sum haplotype ancestry
merged_ancestry_id <- ancestry_chr[c(FALSE,TRUE),] + ancestry_chr[c(TRUE, FALSE),]

#write table

filename <- paste0("tmp/ancestry_012_chr", i, ".txt")
write.table(merged_ancestry_id, filename, quote = FALSE, col.names=FALSE,row.names=FALSE, sep=" ")
}


# cbind tables 
setwd("~/PostDoc_Ifremer/02_AquaExcel/06_Loter/03_result/01_Loter_res")
files <- list.files(path="tmp/", pattern=paste0(".txt")) 

df    <- do.call(cbind,lapply(paste0("tmp/",files),function(fn) fread(fn,header=FALSE, sep=" ")))

  #Id ind
id<-read.table("chr1/admixed_id.txt", header=F)
df1 <- cbind(id,df)
write.table(df1, "genotyp_imp_haplotype.txt", quote = FALSE, col.names=F,row.names=FALSE, sep="\t")




