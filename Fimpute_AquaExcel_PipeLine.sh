###############################################################
###
###	                 	Fimpute Pipeline 
###
##############################################################
#16/07/2021
#Projet: AquaExcel
#Maeva
##########


### Inpute File


	#filtre vcf with maf<0.05
	#qsub -I -l walltime=03:00:00 -l ncpus=1 -l mem=8g

	 vcftools --vcf 02_data_puces/reseq_polyhigh.vcf --maf 0.05  --recode --out 05_Fimpute/01_data_file/reseq_polyhigh_filtered_0.05
	
	# conversion du vcf en 012 avec filtre
	
	 vcftools --vcf 05_Fimpute/01_data_file/reseq_polyhigh_filtered_0.05.recode.vcf --012 --out 05_Fimpute/01_data_file/AquaExcel_filtered_0.05

## Genotype file
	
	# recode in 0,1, 2 and 5 for missing 
	
	01_utility_script/vcf2Fimpute.R 
	
	### Ajouter pere 1 a 5 manquant + 1 mere avec que des 5 pour les genotypes manquant 
	   	#R script 01_utility_script//create_missing_parents_geno.R
	
	
	
## SNPs info file 	

	gawk 'BEGIN {print "SNP_ID" "\t" "Chrom" "\t" "pos" "\t" "Chip"}  {print $1"_"$2 "\t" $1 "\t" $2 "\t" NR}' AquaExcel_filtered_0.05.012.pos > ../02_input_file/AquaExcel_SNP_Fimpute.txt
	

## Pedigree File from APIS + add sex info (convert 1 into 'M' and 2 into 'F')
	# pedigree_Fimpute_AquaExcel.txt

## parameters file
	#save_genotype pour uniquement inference des genotype,si phasing enlever l'option

input=02_input
out=03_output

	echo "title=\"population based imputation\";
	genotype_file=\"${input}/AquaExcel_Fimpute.geno\";
	snp_info_file=\"${input}/AquaExcel_SNP_Fimpute.txt\";
	ped_file=\"${input}/pedigree_Fimpute_AquaExcel.txt\";
	output_folder=\"${out}/AquaExcel_FIMPUTE_OUTPUT\";
	parentage_test /find_match_cnflt;
	save_genotype;
	njob=10;" > AquaExcel_genotype_FImpute_param.ctr


## Run Fimpute
	
	Fimpute3 AquaExcel_genotype_FImpute_param.ctr


	