
## AquaExcel 
## Fevrier 2022

# working directory

rm(list=ls())

# Librairies
library(plyr)
library(dplyr)
library(ggplot2)
library("lme4")
library(MuMIn)
library(reshape2)
library(data.table)
library(tidyverse)
library(here)
library(knitr)
library(MASS)
library(zoo)

#### input data

  # Ancestry rate


EST<-fread("/Intro_ATL_in_EST.txt", header=T)
WEM<-fread("/Intro_ATL_in_WEM.txt", header=T)


  #recombination rate

est_rec<-fread("~/est_windows_2e6.txt", header=T)
names(est_rec)[1]<-"SNP_ID"
names(est_rec)[3]<-"SNP_pos"

  #merge files

est<-merge(est_rec, EST, by=c("SNP_ID", "SNP_pos"))
est<- est[ order(est$chr, est$SNP_pos), ]
write.table(est, "01_data/est_ancestry_rec_windows.txt", quote = FALSE, col.names=TRUE,row.names=FALSE, sep="\t")

wem<-merge(est_rec, WEM, by=c("SNP_ID", "SNP_pos"))
wem<- wem[ order(wem$chr, wem$SNP_pos), ]
write.table(wem, "01_data/wem_ancestry_rec_windows.txt", quote = FALSE, col.names=TRUE,row.names=FALSE, sep="\t")



#################### GLM by LGs ###

https://delladata.fr/introduction-aux-glmm-avec-donnees-de-proportion/
  

#est
    est=na.omit(est)
    est$scaleRec<-scale(est$Mean_EST_Rec)
    m<-mean(est$Mean_EST_Rec)
    sd<-sd(est$Mean_EST_Rec)
  
   # glmer with binomiale     
    y = cbind(est$nATL, est$nEST)
    mod00 <- glmmPQL(y~scaleRec, random=~1 + scaleRec| chr, family=quasibinomial, data=est)
    mod11 <- glmmPQL(y~1, random=~1|chr, family=quasibinomial, data=est)
    
    #extraire effet aleatoire + fixef : effet fixe du model : pour avoir leffet de chaque chromosome (pente de correlation /chr)
    ranef(mod00)$scaleRec+fixef(mod00) #deviation / slpe
    ranef(mod00)$scaleRec+fixef(mod00)[2] 
    
    summary(mod00)
    
    #slope prediction for plot

      predict(mod00, newdata = scaleRec, level=0, type="response")
      scaleRec<-unique(est$scaleRec)
      uniq_rec<-unique(est$scaleRec)
   
    
     # conf interval  + prediction de la variable reponse
     newdat1<-data.frame(scaleRec=unique(est$scaleRec))
     mm<-model.matrix(~scaleRec,newdat1)
     y<-mm%*%mod00$coefficients$fixed
     pvar1 <- diag(mm %*% tcrossprod(mod00$varFix,mm))
     tvar1 <- pvar1+sd(mod00$coefficients$random$chr) ## must be adapted for more complex models
     newdat1 <- data.frame(
       scaleRec=newdat1$scaleRec,
       y= 1/(1+exp(-(y))),
       plo =  1/(1+exp(-(y-1.96*sqrt(pvar1))))
       , phi = 1/(1+exp(-(y+1.96*sqrt(pvar1))))
       , tlo =  1/(1+exp(-(y-1.96*sqrt(tvar1))))
       , thi =  1/(1+exp(-(y+1.96*sqrt(tvar1))))
     )
     
     
  #plot Est    
    
pest<- ggplot(est, aes(scaleRec, ancestry_freq)) +
          geom_point( show.legend = T, size=0.8, color= "darkorange") +
          #geom_line(data=modele, aes(x=uniq_rec,y=perdict), color="black")+
          geom_line(data=newdat1, aes(x=uniq_rec,y=y), color="black")+
          geom_line(data=newdat1, aes(x=uniq_rec,y=tlo), color="grey")+
          geom_line(data=newdat1, aes(x=uniq_rec,y=thi), color="grey")+
          ggtitle("EST population") + theme_bw() + ylab("Atlantic ancestry \n") +
          ylim(0, 0.2)+
          xlab("\n Mean EST recombination rate in cM/Mb")+
          scale_x_discrete(limits=c(-1.17,-0.301, 0.569, 1.44,2.31, 3.18,4.05 , 4.927 ), labels=c("0", "2", "4", "6","8",  "10", "12", "14"))+
          annotate("text", x= 4, y = 0.12, label = "p < 0.001", size=4)

#### modeles LG/LG 
datalist = list()
res <- data.frame(matrix(ncol = 8, nrow = 24))
names(res)<-c("LG", "chr", "estimates", "IC_2.5","IC_97.5", "p_value", "R2c", "R2m")
res$LG<- seq(1:24)

for (i in 1:24) {
      df <- subset(est, chr == i)
      y1 = cbind(df$nATL, df$nEST)
      mod0<-glm(y1~scaleRec , family=binomial, data=df) 
      summary(mod0)
        #model value
       res[i, 2]<-as.character(unique(df$map))
       res[i, 3]<-as.character(summary(mod0)$coefficients[2, 1])
       res[i, 4]<- as.character(confint(mod0, method="Wald")[2])
       res[i, 5]<- as.character(confint(mod0, method="Wald")[4])
       res[i, 6]<-as.character(summary(mod0)$coefficients[2, 4])
       res[i, 7]<-as.character(r.squaredGLMM(mod0)[3])
       res[i, 8]<-as.character(r.squaredGLMM(mod0)[2])
       res[i, 9]<-"EST"
       #predit model pour courbe
      newdata=data.frame(scaleRec=unlist(unique(df$scaleRec)))
    
      preddat <- predict(mod0,
                         type = "link",
                         newdata=newdata,
                         se.fit=TRUE) %>% 
        as.data.frame()
               # model object mod0 has a component called linkinv that 
               # is a function that inverts the link function of the GLM:
      preddat$lower = mod0$family$linkinv(preddat$fit - 1.96*preddat$se.fit) 
      preddat$point.estimate = mod0$family$linkinv(preddat$fit) 
      preddat$upper = mod0$family$linkinv(preddat$fit + 1.96*preddat$se.fit)
      preddat$chr=i
      preddat$map=as.character(unique(df$map))
      preddat$scale_rec=unique(df$scaleRec)
      datalist[[i]] <- preddat
}

    big_data = do.call(rbind, datalist)
    write.table(res, "GLM_EST_estimates_byLGs.txt", quote = FALSE, col.names=TRUE,row.names=FALSE, sep="\t")
    
    #plot LG by LG  
    est$map <- factor(est$map, levels = c(unique(est$map)))
    big_data$map <- factor(big_data$map, levels = c(unique(big_data$map)))
    
    lg_pest<-ggplot(est, aes(scaleRec, ancestry_freq)) +
                geom_point( show.legend = T, size=0.8, color= "darkorange") +
                
                geom_line(data=big_data, aes(x=scale_rec,y=lower), color="grey")+
                geom_line(data=big_data, aes(x=scale_rec,y=upper), color="grey")+
                geom_line(data=big_data, aes(x=scale_rec,y=point.estimate), color="black")+
                ggtitle("EST population") + theme_bw() + ylab("Atlantic ancestry \n") +
                xlab("\n Mean EST recombination rate in cM/Mb")+
                scale_x_discrete(limits=c(-1.17,-0.301, 0.569, 1.44,2.31, 3.18,4.05 , 4.927 ), labels=c("0", "2", "4", "6","8",  "10", "12", "14"))+
      facet_wrap(~ map, nrow=3)

 
### wem
    wem=na.omit(wem)
    wem$scaleRec<-scale(wem$Mean_EST_Rec)
    m<-mean(wem$Mean_EST_Rec)
    sd<-sd(wem$Mean_EST_Rec)
    
    y = cbind(wem$nATL, wem$nEST)
       
    #model 
    mod00 <- glmmPQL(y~scaleRec, random=~1 + scaleRec| chr, family=quasibinomial, data=wem)
    mod11 <- glmmPQL(y~1, random=~1|chr, family=quasibinomial, data=wem)
    
    #extraire effet aleatoire + fixef : effet fixe du model : pour avoir leffet de chaque chromosome (pente de correlation /chr)
    ranef(mod00)$scaleRec+fixef(mod00) #deviation / pente
    ranef(mod00)$scaleRec+fixef(mod00)[2] #valeur recalculer de la pente / chr
    
    summary(mod00)
    
    # interval conf + prediction de la variable reponse
    newdat<-data.frame(scaleRec=unique(wem$scaleRec))
    mm<-model.matrix(~scaleRec,newdat)
    y<-mm%*%mod00$coefficients$fixed
    pvar1 <- diag(mm %*% tcrossprod(mod00$varFix,mm))
    tvar1 <- pvar1+sd(mod00$coefficients$random$chr) ## must be adapted for more complex models
    newdat <- data.frame(
      scaleRec=newdat$scaleRec,
      y= 1/(1+exp(-(y))),
      plo =  1/(1+exp(-(y-1.96*sqrt(pvar1))))
      , phi = 1/(1+exp(-(y+1.96*sqrt(pvar1))))
      , tlo =  1/(1+exp(-(y-1.96*sqrt(tvar1))))
      , thi =  1/(1+exp(-(y+1.96*sqrt(tvar1))))
    )
    
    #plot wem
pwem<- ggplot(wem, aes(scaleRec, ancestry_freq)) + geom_point( show.legend = T, size=0.8, color= "limegreen") +
      #geom_line(data=modele, aes(x=uniq_rec,y=perdict), color="black")+
      geom_line(data=newdat, aes(x=uniq_rec,y=y), color="black")+
      geom_line(data=newdat, aes(x=uniq_rec,y=tlo), color="grey")+
      geom_line(data=newdat, aes(x=uniq_rec,y=thi), color="grey")+
      ggtitle("WEM population") + theme_bw() + ylab("Atlantic ancestry \n") +
      ylim(0, 0.8)+
      xlab("\n Mean EST recombination rate in cM/Mb")+
      scale_x_discrete(limits=c(-1.17,-0.301, 0.569, 1.44,2.31, 3.18,4.05 , 4.927 ), labels=c("0", "2", "4", "6","8",  "10", "12", "14"))+
      annotate("text", x= 3.8, y = 0.75, label = "p < 0.001", size=4)
    
#### modeles LG/LG 
datalist = list()
res <- data.frame(matrix(ncol = 8, nrow = 24))
names(res)<-c("LG", "chr", "estimates", "IC_2.5","IC_97.5", "p_value", "R2c", "R2m")
res$LG<- seq(1:24)

    for (i in 1:24) {
       
      df <- subset(wem, chr == i)
      y1 = cbind(df$nATL, df$nEST)
      mod0<-glm(y1~scaleRec , family=binomial, data=df) 
      summary(mod0)
      res[i, 2]<-as.character(unique(df$map))
      res[i, 3]<-as.character(summary(mod0)$coefficients[2, 1])
      res[i, 4]<- as.character(confint(mod0, method="Wald")[2])
      res[i, 5]<- as.character(confint(mod0, method="Wald")[4])
      res[i, 6]<-as.character(summary(mod0)$coefficients[2, 4])
      res[i, 7]<-as.character(r.squaredGLMM(mod0)[3])
      res[i, 8]<-as.character(r.squaredGLMM(mod0)[2])
      res[i, 9]<-"WEM"
      
      
      newdata=data.frame(scaleRec=unlist(unique(df$scaleRec)))
      
      preddat <- predict(mod0,
                         type = "link",
                         newdata=newdata,
                         se.fit=TRUE) %>% 
        as.data.frame()
      # model object mod0 has a component called linkinv that 
      # is a function that inverts the link function of the GLM:
      preddat$lower = mod0$family$linkinv(preddat$fit - 1.96*preddat$se.fit) 
      preddat$point.estimate = mod0$family$linkinv(preddat$fit) 
      preddat$upper = mod0$family$linkinv(preddat$fit + 1.96*preddat$se.fit)
      preddat$chr=i
      preddat$map=as.character(unique(df$map))
      preddat$scale_rec=unique(df$scaleRec)
      datalist[[i]] <- preddat
    }
    
    big_data = do.call(rbind, datalist)
    write.table(res, "GLM_WEM_estimates_byLGs.txt", quote = FALSE, col.names=TRUE,row.names=FALSE, sep="\t")
    
#plot LG by LG  
    
    wem$map <- factor(wem$map, levels = c(unique(wem$map)))
    big_data$map <- factor(big_data$map, levels = c(unique(big_data$map)))

lg_pwem<-ggplot(wem, aes(scaleRec, ancestry_freq)) +
              geom_point( show.legend = T, size=0.8, color= "limegreen") +
              facet_wrap(~map, nrow=3)+
              geom_line(data=big_data, aes(x=scale_rec,y=lower), color="grey")+
              geom_line(data=big_data, aes(x=scale_rec,y=upper), color="grey")+
              geom_line(data=big_data, aes(x=scale_rec,y=point.estimate), color="black")+
              ggtitle("WEM population") + theme_bw() + ylab("Atlantic ancestry \n") +
              xlab("\n Mean EST recombination rate in cM/Mb")+
              scale_x_discrete(limits=c(-1.17,-0.301, 0.569, 1.44,2.31, 3.18,4.05 , 4.927 ), labels=c("0", "2", "4", "6","8",  "10", "12", "14"))
            ggsave("GLM_LG_WEM_Rec_Ancestry.pdf", lg_pwem, width = 15, height = 8, limitsize = F)              
            
