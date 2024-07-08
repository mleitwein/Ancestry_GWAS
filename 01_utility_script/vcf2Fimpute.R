## Format data: VCF to Fimpute

# deps
library(stringr, quietly = TRUE)
library(vcfR, quietly = TRUE)
library(parallel, quietly = TRUE)

## input
data_dir1 <- file.path("02_data_puces/")
data_dir<-file.path("05_Fimpute/")
vcf_file <- file.path(data_dir, "reseq_polyhigh_filtered_0.05.recode.vcf")
pedigree_file <- file.path(data_dir, "pedigree_Fimpute_Aquaexcel.txt")

## output
output_dir <- file.path(data_dir, "result")
output_file <- file.path(output_dir, "AquaExcel_Fimpute.geno")
indiv_file <- file.path(output_dir, "AquaExcel_Fimpute.indv")
pos_file <- file.path(output_dir, "AquaExcel_Fimpute.pos")

## pedigree data
pedigree_data <- read.csv(pedigree_file, header = TRUE, sep = "\t")

## read VCF file
vcf <- read.vcfR(vcf_file)

## indiv id
indiv_id <- colnames(vcf@gt)[-1]

# correct indiv id
indiv_id <- str_remove(indiv_id, "_.*\\.CEL")

# extract missing parents
missing_parent_id <- sort(c(
    unique(pedigree_data$sire[!pedigree_data$sire %in% indiv_id]),
    unique(pedigree_data$dam[!pedigree_data$dam %in% indiv_id])
))

## function to convert genotype (vectorized version)
convert_gt <- Vectorize(
    vectorize.args = "gt",
    FUN = function(gt) {
        return(switch(
            gt,
            "0/0" = "0",
            "0/1" = "1",
            "1/1" = "2",
            "5"
        ))
    }
)

# convert data from VCF
vcf_data <- Reduce(
    "rbind",
    mclapply(
        colnames(vcf@gt)[-1],
        function(item) {
            return(data.frame(
                IID = item,
                Chip = "1",
                Call = str_c(unname(convert_gt(vcf@gt[,item])), collapse = ""),
                stringsAsFactors = FALSE
            ))
        },
        mc.cores = 1 #4 on linux
    )
)

if(length(unique(sapply(vcf_data$Call, str_length))) != 1)
    stop("issue with data extracted from the VCF")

# number of loci
n_loci <- unique(sapply(vcf_data$Call, str_length))

# fix indiv id
vcf_data$IID <- str_remove(vcf_data$IID, "_.*\\.CEL")

# add missing parent
vcf_data <- rbind(
    vcf_data,
    data.frame(
        IID = missing_parent_id,
        Chip = "1",
        Call = strrep("5", n_loci),
        stringsAsFactors = FALSE
    )
)

# file.geno
write.table(
    vcf_data, 
    file = output_file,
    col.names = TRUE, row.names = FALSE, quote = FALSE
)

# file.indv
write.table(
    data.frame(indiv_id), 
    file = indiv_file,
    col.names = FALSE, row.names = FALSE, quote = FALSE
)

# file.pos
write.table(
    vcf@fix[,c("CHROM", "POS")], 
    file = pos_file,
    col.names = FALSE, row.names = FALSE, quote = FALSE
)
