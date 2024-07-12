########################################################################
# SCRIPT script_gemma_qqplot_manhattan_plot.R
# Version R : 4.1.0
#
# Utilisation : Analyses des resultats de GWAS
#			          Graphiques
# Emilie Delpuech
#########################################################################

library(data.table)
library(ggplot2)
library(tidyverse)
library(viridis)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(dplyr)

set.seed(42)

setwd("~/PostDoc_Ifremer/02_AquaExcel/12_gwas_phenotype/02_resultats/gwas_haplotypes/test_residus")

RES_GWAS <- as.data.frame(fread("gwas_univariate_gras_Haplo_mod7.assoc.txt" , header=T,sep="\t",stringsAsFactors = F)) %>%
  dplyr::mutate(log10_pval = -log10(p_wald),
                chr_LG=sapply(strsplit(rs, split='_', fixed=TRUE), function(x) (x[1])),
                chisq = qchisq(.$p_wald,df=1,lower.tail = F)) 

##########################################################################
### QQplot
##########################################################################


GGP_qqplot <- ggplot(data.frame(x=-log10(sort(RES_GWAS$p_wald)),y=-log10(ppoints(length(RES_GWAS$p_wald)))),aes(y,x)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, alpha = 0.85,col="darkred") +
  scale_x_continuous(limits=c(0,ceiling(max(-log10(RES_GWAS$p_wald),-log10(ppoints(length(RES_GWAS$p_wald))))))) +
  scale_y_continuous(limits=c(0,ceiling(max(-log10(RES_GWAS$p_wald),-log10(ppoints(length(RES_GWAS$p_wald))))))) +
  xlab("p-value attendue") +
  ylab("p-value observee") +
  # ggtitle("Ifremer") +
  theme_bw(10) +
  theme(text = element_text(face = "bold"))

 ggsave(
   filename = paste0("qqplot_residu_ATL_mod7.png"),
   GGP_qqplot,
   width = 300/300, height = 300/300, dpi = 400,
   scale = 7
 )
print("QQplot OK")

##########################################################################
### Manhattan plot
##########################################################################




value_max <- ceiling(max(RES_GWAS$log10_pval))

## Manhattan plot avec ggplot
don <- RES_GWAS %>%
  # Compute chromosome size
  dplyr::group_by(chr) %>% 
  dplyr::filter(!chr%in%25) %>% 
  dplyr::summarise(chr_len=max(ps)) %>%
  # Calculate cumulative position of each chromosome
  dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  # Add this info to the initial dataset
  dplyr::left_join(RES_GWAS, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  dplyr::arrange(chr, ps) %>%
  dplyr::mutate(BPcum=ps+tot) %>% 
  # Add highlight and annotation information
  # dplyr::mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  dplyr::mutate( is_annotate=ifelse(log10_pval>-log10(0.05/nrow(RES_GWAS)), "yes", "no"))
  #dplyr::mutate( is_annotate=ifelse(log10_pval>-log10(0.05/mean(c(table(RES_GWAS$chr)))), "yes", "no"))
# %>% 
# Remove snp with pvalue lower than 1
# dplyr::filter(log10_pval>1)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) %>% left_join(RES_GWAS %>% dplyr::select(chr,chr_LG), by="chr")
log10P <- expression(paste("-log"[10], plain(Pvalue)))

RES_GWAS %>%
  # Compute chromosome size
  dplyr::group_by(chr) %>% nrow


### QTL ###
# QTL ancestry ATL 
RES_GWAS_ATL <- as.data.frame(fread("gwas_univariate_gras_Haplo_mod9.assoc.txt" , header=T,sep="\t",stringsAsFactors = F)) %>%
  dplyr::mutate(log10_pval = -log10(p_wald),
                chr_LG=sapply(strsplit(rs, split='_', fixed=TRUE), function(x) (x[1])),
                chisq = qchisq(.$p_wald,df=1,lower.tail = F)) 
qtl <- RES_GWAS_ATL[ which(RES_GWAS_ATL$log10_pval > 5.9793 & RES_GWAS_ATL$chr_LG!="UN"), ]
qtl<- qtl[,c("chr_LG", "ps" )]

m<-merge(qtl, don[,c("chr_LG", "ps", "BPcum")])

min<-aggregate( BPcum~chr_LG, m, function(x) min(x))
max<-aggregate( BPcum~chr_LG, m, function(x) max(x))

suqtl<-merge(min, max, by="chr_LG")
names(suqtl)[2]<-"xmin"
names(suqtl)[3]<-"xmax"
suqtl$ymin<-"0"
suqtl$ymax<-"15"
suqtl$x<-(suqtl$xmax+suqtl$xmin)/2


## Ajout des regions RI Reproductive Isolation = region genome avec % ancestry ATL pop WEM > 0.5 % 

WEM<-fread("~/GLM/Intro_ATL_in_WEM.txt", header=T)
WEM$chr_LG=sapply(strsplit(WEM$SNP_ID, split='_', fixed=TRUE), function(x) (x[1]))
  
 #Keep value ancestry <0.05
WEM1<-WEM[which(ancestry_freq < 0.05),]
WEM2<-WEM1[,c("SNP_ID")]
names(WEM2)[1]<-"rs"
ri<-merge(WEM2, don[,c("rs", "BPcum")], all.x=T)



ggP_global <- ggplot(don, aes(x=BPcum, y=log10_pval)) +
  ## Show all points:
  # geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.5) +
  geom_vline(xintercept = ri$BPcum,colour = "mistyrose2", size=0.9, alpha=0.9) +
  geom_point( aes(color=as.factor(chr))) +
  # geom_point(aes(color=as.factor(chr)), size=0.8) +
  scale_color_manual(values = rep(viridis(2,begin=0.1,end=0.9),20)) +
  geom_hline(yintercept = -log10(0.05/nrow(don)),colour = "red", linetype="dashed") +
  #geom_hline(yintercept = -log10(0.05/mean(c(table(don$chr)))),colour = "black", linetype="dashed") +
  ## custom X & Y axis:
  scale_x_continuous( label = unique(axisdf$chr_LG), breaks= unique(axisdf$center)) +
  ## Value max choisen :
  # scale_y_continuous(limits=c(0,MAX_PVAL_POP),expand = c(0, 0.1),breaks = c(seq(0,20,5)), labels = c(seq(0,20,5))) + # remove space between plot area and x axis, and last point (y axis) and the title
  # expand_limits(y=c(0, MAX_PVAL_POP)) +
  scale_y_continuous(limits=c(0,value_max+1),expand = c(0, 0.1) ) + # remove space between plot area and x axis, and last point (y axis) and the title
  expand_limits(y=c(1, value_max)) +
  
  ## custom names X & Y axis:
  labs(x="Chromosomes", y=log10P , title="") +
  # Add highlighted points
   geom_point(data=subset(don, is_annotate=="yes"), color="darkorange2", size=1.5) +
  # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=rs), size=2)+
  geom_vline(xintercept = suqtl$x,colour = "dodgerblue", size=0.9, alpha=0.6) +
    
  #geom_rect(data=suqtl, inherit.aes=FALSE,
   #         aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
    #        ), alpha=0.3) +
  
  ## Custom the theme:
  theme_minimal(10) +
  theme(
    # text = element_text(face = "bold", size=12),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) 
  )

# plot_list[[pop]] <- ggdraw() +
plot.with.inset <- ggdraw() +
  draw_plot(ggP_global) +
  # draw_plot(GGP_qqplot, x = 0.05, y = .6, width = .3, height = .3)
  draw_plot(GGP_qqplot, x = 0.3, y = .9, width = .3, height = .3)

ggsave("Manhattan_mod7_qtl_ancestry_RI.pdf",
       # plot.with.inset,
       ggP_global,
       width = 8,
       height = 4)
	   
	   
print("Manhattan OK")


## grep outliers and save file
#y=(-log10(0.05/mean(c(table(don$chr)))))
y=(-log10(0.05/nrow(don)))
qtl <- don[ which(don$log10_pval >= y & don$chr_LG!="UN"),]

write.table(qtl, "qtl_gras_seul.txt", quote = FALSE, col.names=TRUE,row.names=FALSE, sep="\t")




###### PLOT correlation QTL value vs RI ###

#region RI
RI<-WEM[which(ancestry_freq < 0.05),]

#region sans RI
sansRI<-WEM[which(ancestry_freq > 0.05),]
m=mean(sansRI$ancestry_freq)

# calcul du me

WEM$RIint<-1-(WEM$ancestry_freq/m)
wemMe<-as.data.frame(WEM[,c(1,8)])
names(wemMe)[1]<-"rs"

# merge data avec value QTL
df<-merge(wemMe, don, by="rs")

# plot 

ggplot(df, aes(RIint, log10_pval)) + geom_point(show.legend = T, size=0.7) +
    theme_bw()
