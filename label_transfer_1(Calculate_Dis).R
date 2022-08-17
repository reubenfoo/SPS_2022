CCLE <-read.csv("./ccle_80SPS_4Q.csv",header = T,row.names = 1,check.names=F)
clinical <-read.csv('C:/Users/SIQI/Desktop/619MOREGUANFANG.csv',header = T,row.names = 1,check.names=F)
CCLE_Basal<-subset(CCLE,(CCLE$subtype == "Luminal"))#| (CCLE$subtype == "Basal A"))
clinical_Basal<-subset(clinical,clinical$Subtype == 'LumA')

CCLE_Basal_TR<-subset(CCLE,CCLE$Q == "TR")                      
CCLE_Basal_TL<-subset(CCLE,CCLE$Q == "TL")
CCLE_Basal_BR<-subset(CCLE,CCLE$Q == "BR")
CCLE_Basal_BL<-subset(CCLE,CCLE$Q == "BL")

clinical_Basal_ge <- clinical[,1:80]
CCLE_ge<-CCLE[,1:80]
#only gene expression
CCLE_Basal_TR_ge <- CCLE[,1:80]

CCLE_Basal_TL_ge <- CCLE_Basal_TL[,1:80]

CCLE_Basal_BR_ge <- CCLE_Basal_BR[,1:80]
CCLE_Basal_BL_ge <- CCLE_Basal_BL[,1:80]

library(rdist)
dis<-cdist(clinical_Basal_ge,CCLE_ge)
dis_Basal_TR <-cdist(clinical_Basal_ge,CCLE_Basal_TR_ge)
dis_Basal_TL <-cdist(clinical_Basal_ge,CCLE_Basal_TL_ge)

dis_Basal_BR <-cdist(clinical_Basal_ge,CCLE_Basal_BR_ge)
dis_Basal_BL <-cdist(clinical_Basal_ge,CCLE_Basal_BL_ge)


dis_mean_BR <- rowMeans(dis_Basal_BR)
dis_mean_BL <- rowMeans(dis_Basal_BL)
dis_mean_TR <- rowMeans(dis_Basal_TR)
dis_mean_TL <- rowMeans(dis_Basal_TL)
df <- data.frame(dis)
df <- data.frame(dis_mean_BR, dis_mean_BL,dis_mean_TR,dis_mean_TL)
write.csv(df,file = "C:/Users/SIQI/Desktop/all.csv")
