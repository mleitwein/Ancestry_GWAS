# Ancestry_GWAS
The fate of a polygenic phenotype within the genomic landscapes of introgression in the European seabass hybrid zone. 
Maeva Leitwein, G. Durif, E. Delpuech, P.A Gagnaire, B. Ernande, M. Vandeputte, A. Vergnet, M. Duranton, F. Clota and F. Allal


# Warning
The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.



## Imputation & phasing

01_Fimpute_AquaExcel_PipeLine.sh


## Ancestry analysis with Lotter
 
	# Utility script file format Fimpute to Loter
	01_utility_script/Fimpute2Loter.R

	#Run Loter for each chromosome 
	#run_loter.sh 
	for i in {1..25}
	do
	loter_cli -r chr$i/ATL_ref_hap.txt chr$i/EST_ref_hap.txt -a chr$i/admixed_hap.txt -f txt -o chr$i/ancestry_chr$i.txt -n 8 -pc -v
	done

## Introgression rate and tracts analysis see script on https://github.com/mleitwein/local_ancestry_inference_with_ELAI

## GLM models 

02_GLM_ancestry_rec.R
			
## Gwas with ancestry 
	#GEMMA file and run
	03_Gwas_ancestry.sh	
	#plot
	01_utility_script/script_gemma_qqplot_manhattan_plot.R

	
