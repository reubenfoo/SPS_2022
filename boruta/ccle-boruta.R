library(Boruta)
library(ggplot2)
library(reshape2)
library(dplyr)
library(caret)
library(gridExtra)
library('testit')
library(MLeval)

set.seed(111)

ccle.processed <- read.table('./CCLE_processed_dedup.txt', header=TRUE, check.names=FALSE)

##### Importing SPS genes #####

sps.genes <- unlist(read.table('SPS_genesymbol.txt', header=TRUE))
sps.genes[c(25:26)] <- c('H2AX','H2AZ1') ## alternative gene naming
assert(length(which(sps.genes %in% ccle.processed$gene_name)) == 81)

##### Importing PAM genes #####

pam.genes <- unlist(read.table('PAM50_genesymbol.txt'))
pam.genes[c(11,26,40)] <- c('NUF2','NDC80','ORC6') ## alternative gene naming
assert(length(which(pam.genes %in% ccle.processed$gene_name)) == 50)

metadata <- read.table('./CCLE_metadata_reordered.txt', header=TRUE)
data1 <- read.csv('./Expression_22Q1_Public_subsetted.csv')[1:6]
metadata <- left_join(metadata, data1, by=c('BROAD_ID' = 'depmap_id'))

##### Filter out breast cancer cell lines, transpose, remove zero-variance features and add labels  #####

breast.cancer.ids <- grep('BREAST',colnames(ccle.processed))
ccle.processed <- ccle.processed[,c(1:2,breast.cancer.ids,ncol(ccle.processed)-1,ncol(ccle.processed))]

dataset <- data.frame(t(ccle.processed[,3:length(ccle.processed),drop=F]),check.names=F)
dataset <- dataset[,apply(dataset[1:51,], 2, var) != 0,drop=F] # drop zero-variance genes
dataset.meta <- data.frame(t(dataset[52:53,]))
dataset <- dataset[-c(52:53),]
dataset[,1:ncol(dataset)] <- lapply(dataset[,1:ncol(dataset)], FUN = function(y){as.numeric(y)})

##### Min-max normalisation of dataset  #####

normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}
dataset <- data.frame(t(apply(dataset,1,normalize)))


dataset.meta$entrez_id <- as.integer(dataset.meta$entrez_id)
na.ids <- which(is.na(dataset.meta$gene_name) | dataset.meta$gene_name == '')
dataset.meta$gene_name <- replace(dataset.meta$gene_name, na.ids, ccle.processed$gene_id)


dataset$labels <- metadata$tissue_status[(breast.cancer.ids-2)]
dataset$labels[which(dataset$labels %in% c('Normal','Tumour'))] <- 'Non.Metastasis'
dataset$labels <- as.factor(dataset$labels)

dataset$subtype <- metadata$lineage_4[(breast.cancer.ids-2)]
dataset$subtype[which(dataset$subtype %in% c('Basal A','Basal B'))] <- 'Basal'
dataset$subtype[which(!dataset$subtype %in% c('Basal'))] <- 'Not.Basal'
dataset$subtype <- as.factor(dataset$subtype)

sps.ids <- which(dataset.meta$gene_name %in% sps.genes)
pam.ids <- which(dataset.meta$gene_name %in% pam.genes)


##################### Boruta on ALL samples #####################


# boruta.bank_train <- Boruta(dataset[,1:(ncol(dataset)-2)], dataset$labels, doTrace = 2, maxRuns=200)
# boruta.bank_train
# 
# boruta_summary <- attStats(boruta.bank_train)
# # write.csv(boruta_summary,'boruta-ccle-breast-200runs-second.csv')
# # write.csv(boruta.bank_train$ImpHistory,'boruta-ccle-breast-200runs-history-second.csv')
# boruta_summary <- read.csv('boruta-ccle-breast-200runs-second.csv')
# imp.history <-read.csv('boruta-ccle-breast-200runs-history-second.csv')
# imp.history <- imp.history[,-1,drop=F]
# 
# selected.features.ids <- which(boruta_summary$decision %in% c('Confirmed'))
# selected.features <- boruta_summary[selected.features.ids,,drop=F]
# 
# selected.features$gene <- dataset.meta$gene_name[selected.features.ids]
# 
# boxplot <- data.frame(imp.history[,c(selected.features.ids,(ncol(imp.history)-2):ncol(imp.history)),drop=F])
# colnames(boxplot)[1:(ncol(boxplot)-3)] <- dataset.meta$gene_name[selected.features.ids]
# boxplot <- melt(boxplot[,c(ncol(boxplot),ncol(boxplot)-1,ncol(boxplot)-2,order(selected.features$medianImp))])
# 
# ggplot(boxplot,aes(x=variable,y=value)) + geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##################### ML Model #####################


# trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 100)
# tunegrid <- expand.grid(
#   .mtry = c(2, 105, 208),
#   .splitrule = "gini",
#   .min.node.size = c(1))
# 
# nb_fit1 <- train(x=dataset[,1:43651,drop=F], y=dataset$labels, method='ranger',trControl=trctrl, tuneGrid=tunegrid, trace=T)
# nb_fit2 <- train(x=dataset[,selected.features.ids,drop=F], y=dataset$labels, method='ranger',trControl=trctrl)
# nb_fit3 <- train(x=dataset[,sps.ids,drop=F], y=dataset$labels, method='ranger',trControl=trctrl)
# 
# print(nb_fit1)
# print(nb_fit2)
# print(nb_fit3)
# 
# model.res1 <- nb_fit1$results
# model.res2 <- nb_fit2$results
# model.res3 <- nb_fit3$results
# 
# write.table(model.res1, 'prelim-all.txt', row.names = F)
# write.table(model.res2, 'prelim-boruta.txt', row.names = F)
# write.table(model.res3, 'prelim-sps.txt', row.names = F)
# 
# p1 <- ggplot(model.res1, aes(x=as.factor(mtry), y=Accuracy)) + geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=Accuracy-AccuracySD, ymax=Accuracy+AccuracySD), width=.2, position=position_dodge(.9)) +
#   coord_cartesian(ylim=c(0,1)) + xlab('All Genes') + ylab('CV Accuracy')
# 
# p2 <- ggplot(model.res2[which(model.res2$splitrule == 'gini'),], aes(x=as.factor(mtry), y=Accuracy)) + geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=Accuracy-AccuracySD, ymax=Accuracy+AccuracySD), width=.2, position=position_dodge(.9)) +
#   coord_cartesian(ylim=c(0,1)) + xlab('Boruta') + ylab('')
# 
# p3 <- ggplot(model.res3[which(model.res3$splitrule == 'gini'),], aes(x=as.factor(mtry), y=Accuracy)) + geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=Accuracy-AccuracySD, ymax=Accuracy+AccuracySD), width=.2, position=position_dodge(.9)) +
#   coord_cartesian(ylim=c(0,1)) + xlab('SPS') + ylab('')
# 
# grid.arrange(p1,p2,p3,ncol=3)


##################### Boruta on 80% train, then 20% test #####################


set.seed(111)

rm(data.partition.ids, boruta.normfeatures.cfms, boruta.normfeatures.tent)

##### Import precalculated partitions (1000 iterations)

data.partition.ids <- read.table('data-partition-ids-survival.txt')

##### Import Boruta features (with Random Forest & 100 runs)

no_col3 <- max(count.fields('boruta-impforest-80train-100runs-normfeatures-confirmed.txt'))
boruta.normfeatures.cfm <- read.table('boruta-impforest-80train-100runs-normfeatures-confirmed.txt',fill=TRUE,header=F,col.names=paste("col", 1:no_col3, sep="_"))
no_col4 <- max(count.fields('boruta-impforest-80train-100runs-normfeatures-tentative.txt'))
boruta.normfeatures.tent <- read.table('boruta-impforest-80train-100runs-normfeatures-tentative.txt',fill=TRUE,header=F,col.names=paste("col", 1:no_col4, sep="_"))

##### Explore Boruta features

# feature.maths <- function(df.features, func){
#   lengths <- apply(df.features[,-1],1,function(y){length(as.numeric(na.omit(unlist(y))))})
#   return(func(lengths))
# }
# feature.maths(boruta.normfeatures.tent, mean)
# 
# lengths <- apply(boruta.normfeatures.cfm[,-1],1,function(y){length(as.numeric(na.omit(unlist(y))))})
# lengths2 <- apply(boruta.normfeatures.tent[,-1],1,function(y){length(as.numeric(na.omit(unlist(y))))})
# mean(lengths + lengths2)

##### Generate 1000 random partitions of 80% train - 20% test

normal.sample.id <- grep('HMEL', rownames(dataset))
# dataset.part <- dataset[-normal.sample.id,]
# train.ids <- createDataPartition(dataset.part$subtype, times=1000, p=0.8)
# for (i in 1:1000){
#   resample <- train.ids[[i]]
#   resample[which(resample >= 32)] <- resample[which(resample >= 32)] + 1
#   if (i > 0){
#     write.table(t(resample), file='data-partition-ids-subtype-nonormal-2.txt', append=T, col.names=F, row.names=F)
#   }
# }
# rm(normal.sample.id, dataset.part, train.ids, resample)

##### Obtain 81 random genes

# random.ids <- sample(1:43651, 81)
# write.table(random.ids, '81-random-genes.txt',row.names=F,col.names=F)

for (i in 1:1){
  print(i)
  
  # train-test 80-20 split with precalculated partitions
  current.partition.ids <- as.integer(data.partition.ids[i,])
  print(current.partition.ids)
  dataset.train <- dataset[current.partition.ids,]
  dataset.test  <- dataset[-current.partition.ids,] ## if splitting by survival label
  # dataset.test  <- dataset[-c(32,current.partition.ids),] ## if splitting by subtype label
  
  ##### Perform Boruta feature selection on training data only (100 runs)
  # boruta.eighty <- Boruta(x=dataset.train[,1:(ncol(dataset.train)-2)], y=dataset.train$subtype, doTrace = 0, maxRuns=100, holdHistory=T)
  # print(boruta.eighty)
  # boruta_summary <- attStats(boruta.eighty)
  # confirmed.features.ids <- which(boruta_summary$decision %in% c('Confirmed'))
  # tentative.features.ids <- which(boruta_summary$decision %in% c('Tentative'))
  # selected.features.ids <- c(confirmed.features.ids,tentative.features.ids)

  ##### Get selected Boruta features
  confirmed.features.ids <- as.numeric(na.omit(unlist(boruta.normfeatures.cfm[i,-1])))
  tentative.features.ids <- as.numeric(na.omit(unlist(boruta.normfeatures.tent[i,-1])))
  selected.features.ids <- c(confirmed.features.ids,tentative.features.ids)
  print(selected.features.ids)
  
  ##### Random Forest training
  tunegrid <- expand.grid(
    .mtry = 208,
    .splitrule = "gini",
    .min.node.size = c(1)
    )

  fitControl <- trainControl(
                              # verboseIter = TRUE,
                              classProbs = TRUE
                              )

  all.fit <- train(x=dataset.train[,1:(ncol(dataset.train)-2),drop=F],y=dataset.train$labels,method='ranger',tuneGrid=tunegrid, trControl=fitControl)
  boruta.fit <- train(x=dataset.train[,selected.features.ids,drop=F], y=dataset.train$labels,method='ranger',tuneLength=5, trControl=fitControl)
  sps.fit <- train(x=dataset.train[,sps.ids,drop=F],y=dataset.train$labels,method='ranger',tuneLength=5, trControl=fitControl)
  
  ##### Get the predictions of your model in the test set
  # predictions1 <- predict(all.fit, dataset.test[,1:(ncol(dataset.test)-2),drop=F])
  # predictions2 <- predict(boruta.fit, dataset.test[,selected.features.ids,drop=F])
  # predictions3 <- predict(sps.fit, dataset.test[,sps.ids,drop=F])
  
  ##### Predictions with probabilities
  predict.probs.1 <- predict(all.fit, dataset.test[,1:(ncol(dataset.test)-2),drop=F],type="prob")
  predict.probs.2 <- predict(boruta.fit, dataset.test[,selected.features.ids,drop=F],type="prob")
  predict.probs.3 <- predict(sps.fit, dataset.test[,sps.ids,drop=F],type="prob")
  
  ##### MLeval package to obtain AUC-ROC
  test1 <- evalm(data.frame(predict.probs.1, dataset.test$labels))
  test2 <- evalm(data.frame(predict.probs.2, dataset.test$labels))
  test3 <- evalm(data.frame(predict.probs.3, dataset.test$labels))
  auc1 <- test1$stdres$Group1['AUC-ROC',1]
  auc2 <- test2$stdres$Group1['AUC-ROC',1]
  auc3 <- test3$stdres$Group1['AUC-ROC',1]

  ##### Get prediction accuracy
  accuracy1 <- sum(max.col(predict.probs.1) == as.numeric(dataset.test$labels))/length(predictions1)
  accuracy2 <- sum(max.col(predict.probs.2) == as.numeric(dataset.test$labels))/length(predictions2)
  accuracy3 <- sum(max.col(predict.probs.3) == as.numeric(dataset.test$labels))/length(predictions3)

  ##### Write scores
  # write.table(t(c(i,auc1, auc2, auc3)), file='survival-auc-everything.txt', append=T, col.names=F, row.names=F)
}

library(ggbeeswarm)

rm(acc.list)

acc.list <- read.table('survival-auc-everything.txt', header=T)[,-1,drop=F]
# acc.list2 <- read.table('subtype-boruta-normrun-pam.txt', header=F)[1:100,-1,drop=F]
# acc.list3 <- read.table('subtype-boruta-normrun-random.txt', header=F)[1:100,-1,drop=F]
# acc.list$PAM <- acc.list2$V2
# acc.list$Random <- acc.list3$V2

wilcox.test(x = acc.list$All, y = acc.list$SPS)

# p1 <- ggplot(melt(acc.list), aes(x=value, fill=variable,colour=variable)) + geom_density(alpha=0.2) +
#   coord_cartesian(xlim=c(0.4,1)) + xlab('Test Accuracy') + ylab('Density') +
#   geom_vline(xintercept=mean(acc.list$SPS), color='green', alpha=0.7, size=0.9) +
#   geom_vline(xintercept=mean(acc.list$Boruta), color='yellow', alpha=0.6, size=0.9) +
#   geom_vline(xintercept=mean(acc.list$All), color='red', alpha=0.7, size=0.9) +
#   geom_vline(xintercept=mean(acc.list$PAM), color='purple', alpha=0.7, size=0.9) +
#   geom_vline(xintercept=mean(acc.list$Random), color='blue', alpha=0.7, size=0.9) +
#   theme(legend.position='none')
  
p2 <- ggplot(melt(acc.list), aes(x=variable, y=value, color=variable)) + 
  # geom_jitter() +
  coord_cartesian(ylim=c(0.15,1)) +
  geom_violin(width=1, alpha=0.4, aes(fill=variable), color="black") +
  geom_boxplot(width=0.1, fill='white', color="black") +
  xlab('Model') + ylab('Test AUC') +
  labs(color='Model', fill='Model')

p2
# grid.arrange(p1,p2,ncol=2,widths=c(6,7))

##### Plot frequency of occurrence of genes in Boruta features

nocol <- max(count.fields('boruta-impforest-80train-100runs-normfeatures-confirmed.txt'))
features.list.cfm <- read.table('boruta-impforest-80train-100runs-normfeatures-confirmed.txt',fill=TRUE,header=F,col.names=paste("col", 1:nocol, sep="_"))
nocol2 <- max(count.fields('boruta-impforest-80train-100runs-normfeatures-tentative.txt'))
features.list.tent <- read.table('boruta-impforest-80train-100runs-normfeatures-tentative.txt',fill=TRUE,header=F,col.names=paste("col", 1:nocol2, sep="_"))
features.list <- dplyr::bind_rows(features.list.cfm, features.list.tent)

features.freq <- data.frame(table(unlist(features.list[,2:ncol(features.list)])))
features.freq <- features.freq[order(features.freq$Freq, decreasing=T),,drop=F]
features.freq$id <- dataset.meta$gene_name[as.integer(as.character(features.freq$Var1))]
features.freq$id[which(features.freq$id == '')] <- paste('ID',apply(features.freq[which(features.freq$id == ''),],1,function(y){y[1]}),sep='.')

ggplot(features.freq, aes(x=id, y=Freq)) + geom_bar(stat='identity', width=1) + xlab('Gene') + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_x_discrete(limits=features.freq$id)

##### Check if data partitions are unbiased (i.e. each sample of the 51 occurs approx. same # of times)

# stuff <- unlist(lapply(c(1:51),function(x) length(which(data.partition.ids == x))))
# ggplot(data.frame(stuff), aes(x=stuff)) + geom_boxplot()
