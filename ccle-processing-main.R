library(ggplot2)
library(gridExtra)
library(factoextra)
library(reshape2)
library(pheatmap)
library('dplyr')
library(testit)

ccle.processed <- read.table('./CCLE_processed_dedup.txt', header=TRUE, check.names=FALSE)
metadata <- read.table('./CCLE_metadata_reordered.txt', header=TRUE)

## 
valid.ids <- which(metadata$tissue_status %in% c('Tumour','Metastasis','Normal'))
metadata <- metadata[valid.ids,,drop=F]
ccle.processed <- ccle.processed[,-(setdiff(1:1019,valid.ids)+2),drop=F]

## importing breast cancer subtype metadata
data1 <- read.csv('./Expression_22Q1_Public_subsetted.csv')
metadata <- left_join(metadata, data1[1:6], by=c('BROAD_ID' = 'depmap_id'))

## importing SPS and PAM genes
sps.genes <- unlist(read.table('SPS_genesymbol.txt', header=TRUE))
sps.genes[c(25:26)] <- c('H2AX','H2AZ1') ## alternative gene naming
assert(length(which(sps.genes %in% ccle.processed$gene_name)) == 81)

pam.genes <- unlist(read.table('PAM50_genesymbol.txt'))
pam.genes[c(11,26,40)] <- c('NUF2','NDC80','ORC6') ## alternative gene naming
assert(length(which(pam.genes %in% ccle.processed$gene_name)) == 50)

# boruta_summary <- read.csv('boruta-ccle-breast-200runs-second.csv')
# boruta.features <- boruta_summary$X[which(boruta_summary$decision %in% c('Confirmed'))]
# ccle.boruta.ids <- substring(boruta.features,2)

# no_col2 <- max(count.fields('boruta-impforest-80train-100runs-normfeatures.txt'))
# boruta.rf.normfeatures <- read.table('boruta-impforest-80train-100runs-normfeatures.txt',fill=TRUE,header=F,col.names=paste("col", 1:no_col2, sep="_"))
# selected.features.ids <- as.numeric(na.omit(unlist(boruta.rf.normfeatures[700,2:no_col2])))

## 26 different types of tissue - sorted in decreasing order
cell.lines <- gsub('^[^_]*_','',metadata$CCLE_ID)
cell.line.labels <- as.data.frame(sort(table(cell.lines),decreasing=TRUE))


######################## PCA on selected cell lines ############################


set.seed(111)

for (k in 1:1){
  print(k)
  vector <- c(5) ## which cell types to use (5 = breast cancer)
  
  for (j in 1:length(vector)) {
    if (j == 1) {
      string <- cell.line.labels$cell.lines[[vector[[j]]]]
    } else {
      string <- paste(string, cell.line.labels$cell.lines[[vector[[j]]]] , sep = '|')
    }
  }
  
  ## filter out cell lines of specified type
  filtered.ids <- grep(string,colnames(ccle.processed))
  ccle.subset <- ccle.processed[,c(1:2,filtered.ids,ncol(ccle.processed)-1,ncol(ccle.processed))]
  ccle.subset[,3:(ncol(ccle.subset)-2)] <- lapply(ccle.subset[,3:(ncol(ccle.subset)-2)], function(x) as.numeric(x))
  ccle.subset <- ccle.subset[apply(ccle.subset[,3:(ncol(ccle.subset)-2)], 1, var) != 0,,drop=F]

  ## filter out genes (SPS, PAM, Boruta, random-81 etc.)
  ccle.sub <- ccle.subset[ccle.subset$gene_name %in% sps.genes,]
  # ccle.sub <- ccle.subset[ccle.subset$gene_name %in% pam.genes,]
  # ccle.sub <- ccle.subset[ccle.boruta.ids,]
  # ccle.sub <- ccle.subset[selected.features.ids,]
  # ccle.sub <- ccle.subset[c(sample(which(ccle.subset$gene_name != '' & !(ccle.subset$entrez %in% sps.genes$Entrez)),81)),]
  # ccle.sub <- ccle.subset # pass-through line

  ## PCA with scaling and centering (and removal of zero-variance rows if necessary)
  ccle.sub <- ccle.sub[apply(ccle.sub[3:(ncol(ccle.sub)-2)], 1, var) != 0,]
  subset.pca <- prcomp(t(ccle.sub[,3:(ncol(ccle.sub)-2)]), scale = TRUE)
    
  # first.two.var <- fviz_eig(subset.pca)$data[c(1:2),2]
  # write.table(first.two.var[[1]], './first-variances-breast-100000-secondrun.txt', append=T, col.names=F, row.names=F)
  # write.table(first.two.var[[2]], './second-variances-breast-100000-secondrun.txt', append=T, col.names=F, row.names=F)
  # curr.p <- plotPCA(subset.pca, "SPS Genes - Breast Only")[[5]]
  # write.table(curr.p, './fisher-p-values-breast-100000-secondrun.txt', append=T, col.names=F, row.names=F)
}


######################## 100000 random runs ############################


# p.val.all <- read.table('./fisher-p-values-breast-100000-secondrun.txt')
# first.var.all <- read.table('./first-variances-breast-100000-secondrun.txt')
# second.var.all <- read.table('./second-variances-breast-100000-secondrun.txt')
# 
# sum(p.val.all <= 0.021821) / 100000      # 4.10% of p-values fall below SPS p-value
# sum(first.var.all > 42.355) / 100000      # 0% exceed PC1 variance
# sum(second.var.all > 8.626) / 100000      # 7.251% exceed PC2 variance
# 
# p1 <- ggplot(p.val.all, aes(x=V1)) + geom_histogram(binwidth=0.05) +
#   xlab('p-value') + ggtitle('Fisher Test p-values') +
#   geom_vline(xintercept = 0.021821, colour='red') +
#   geom_vline(aes(xintercept=median(V1)), color="blue", linetype="dashed")
# 
# p1
# 
# p2 <- ggplot(first.var.all, aes(x=V1)) + geom_histogram(binwidth=0.2) +
#   xlab('PC1 Proportion (%)') + ggtitle('Proportion of Variance Explained by PC1') +
#   geom_vline(xintercept = 42.355, colour='red') +
#   geom_vline(aes(xintercept=median(V1)), color="blue", linetype="dashed")
# 
# p3 <- ggplot(second.var.all, aes(x=V1)) + geom_histogram(binwidth=0.1) +
#   xlab('PC2 Proportion (%)') + ggtitle('Proportion of Variance Explained by PC2') +
#   geom_vline(xintercept = 8.626, colour='red') +
#   geom_vline(aes(xintercept=median(V1)), color="blue", linetype="dashed")
# 
# median(first.var.all$V1)
# 
# grid.arrange(p2,p3,ncol=2)


######################## Plot PCA ############################


# subtypes <- c('Acute Myeloid Leukemia','B-Cell Non-Hodgkin\'s Lymphoma', 'B-Lymphoblastic Leukemia',
#               'Burkitt\'s Lymphoma', 'Chronic Myelogenous Leukemia','Hodgkin\'s Lymphoma','Non-Cancerous',
#               'Other Blood Cancers', 'Plasma Cell Myeloma', 'T-Cell Non-Hodgkin\'s Lymphoma', 'T-Lymphoblastic Leukemia')
            
plotPCA <- function(pca, title){
  pca.plot <- as.data.frame(pca$x[,c('PC1','PC2')])
  pca.plot$labels <- gsub('\\..+$', '', rownames(pca.plot))
  pca.plot$labels2[pca.plot$labels != 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'] <- 'Solid'
  pca.plot$labels2[pca.plot$labels == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'] <- 'Non-Solid'
  
  pca.plot$status <- make.names(metadata$cancer_type_detail[(filtered.ids-2)])
  pca.plot$status2 <- metadata$model_name[(filtered.ids-2)]
  pca.plot$status3 <- metadata$tissue_status[(filtered.ids-2)]
  pca.plot$status4 <- metadata$lineage_3[(filtered.ids-2)]
  pca.plot$status5 <- metadata$lineage_4[(filtered.ids-2)]
  
  median1 <- median(pca.plot$PC1)
  median2 <- median(pca.plot$PC2)
  ids1 <- which(pca.plot$PC1 <= median1)
  ids2 <- which(pca.plot$PC2 <= median2)
  
  # pca.plot <- pca.plot[which(pca.plot$status %in% c('B-Cell Non-Hodgkin\'s Lymphoma','T-Cell Non-Hodgkin\'s Lymphoma',
  #                                                   'Hodgkin\'s Lymphoma','Burkitt\'s Lymphoma')),]
  
  p.eig <- get_eig(pca)
  pc1.var <- round(p.eig[1,2],2)
  pc2.var <- round(p.eig[2,2],2)

  plot1 <- ggplot(pca.plot,aes(x=PC1,y=PC2,label=status2)) + 
    geom_point(size=3, aes(colour=status3)) +
    # geom_text(hjust=0, vjust=0, size=2) +
    geom_hline(aes(yintercept = median2), colour='red') +
    geom_vline(aes(xintercept = median1), colour='red') +
    ggtitle(title) +
    xlab(paste('PC1 -',pc1.var,'%')) +
    ylab(paste('PC2 -',pc2.var,'%')) +
    theme(legend.text = element_text(size=10),
          legend.position = 'top',
          legend.key.height=unit(1.2,"cm")) +
    labs(colour='Survival') 
  
  plot2 <- ggplot(pca.plot,aes(x=PC1,y=PC2,label=status2)) + 
    geom_point(size=3, aes(colour=status5)) +
    # geom_text(hjust=0, vjust=0, size=2) +
    geom_hline(aes(yintercept = median2), colour='red') +
    geom_vline(aes(xintercept = median1), colour='red') +
    ggtitle(title) +
    xlab(paste('PC1 -',pc1.var,'%')) +
    ylab(paste('PC2 -',pc2.var,'%')) +
    theme(legend.text = element_text(size=10),
          legend.position = 'top') +
    labs(colour='Subtype') 
    
  # plot2 <- ggplot(pca.plot,aes(x=PC1,y=PC2)) + 
  #   geom_point(size=2, aes(colour=status5)) +
  #   geom_hline(aes(yintercept = median2), colour='red') +
  #   geom_vline(aes(xintercept = median1), colour='red') +
  #   facet_wrap(~pca.plot$status5) +
  #   ggtitle(title) +
  #   theme(legend.text = element_text(size=10),
  #         legend.position = 'top')
  
  tl <- which(pca.plot$PC1 < median(pca.plot$PC1) & pca.plot$PC2 >= median(pca.plot$PC2))
  tr <- which(pca.plot$PC1 >= median(pca.plot$PC1) & pca.plot$PC2 >= median(pca.plot$PC2))
  bl <- which(pca.plot$PC1 < median(pca.plot$PC1) & pca.plot$PC2 < median(pca.plot$PC2))
  br <- which(pca.plot$PC1 >= median(pca.plot$PC1) & pca.plot$PC2 < median(pca.plot$PC2))
  
  quad.list <- data.frame(matrix(NA, nrow = nrow(pca.plot), ncol = 2))
  quad.list[tl,1] <- 'TL'
  quad.list[tr,1] <- 'TR'
  quad.list[bl,1] <- 'BL'
  quad.list[br,1] <- 'BR'
  quad.list[,2] <- metadata$model_name[(filtered.ids-2)]
  
  quadrants_cw <- c(length(tl),length(tr), length(br), length(bl))
  
  percentages <- c(sum(pca.plot$PC1 < median(pca.plot$PC1) & pca.plot$PC2 >= median(pca.plot$PC2) & pca.plot$status3 %in% c('Metastasis')),
                   sum(pca.plot$PC1 >= median(pca.plot$PC1) & pca.plot$PC2 >= median(pca.plot$PC2) & pca.plot$status3 %in% c('Metastasis')),
                   sum(pca.plot$PC1 >= median(pca.plot$PC1) & pca.plot$PC2 < median(pca.plot$PC2) & pca.plot$status3 %in% c('Metastasis')),
                   sum(pca.plot$PC1 < median(pca.plot$PC1) & pca.plot$PC2 < median(pca.plot$PC2) & pca.plot$status3 %in% c('Metastasis')))
  
  p.val <- fisher.test(rbind(percentages,quadrants_cw-percentages))
  
  return(list(plot1, plot2,
              quadrants_cw, percentages, p.val$p.value, quad.list,
              ids1, ids2))
}

## generate PCA plots and specify title of plots
list <- plotPCA(subset.pca, "All Genes")

## plot PCA with metastatic label and subtype label respectively
grid.arrange(list[[1]],list[[2]],ncol=2)

## other metrics (quadrant count, percentage, and p-value of Fisher Exact test)
quads <- list[[3]]
percent <- list[[4]]
p.value <- list[[5]]


######################## Quadrant Significant Genes ############################


quad.labels <- list[[6]]
quad.labels[,1] <- factor(quad.labels[,1], levels=c('TL','BL','TR','BR'))
order.ids <- order(quad.labels[,1])
quad.labels <- quad.labels[order.ids,]
colnames(quad.labels) <- 'Quadrant'
sps.gene.labels <- as.data.frame(ccle.sub$gene_name) ## SPS genes in order of increasing index
all.gene.labels <- as.data.frame(ccle.subset$gene_name) ## all 46674 genes in order of increasing index
ccle.sorted <- t(scale(t(ccle.sub[,c(order.ids+2)])))
rownames(quad.labels) <- colnames(ccle.sorted)

quads.pairs <- list(c('TL','TR'),c('BL','BR'),c('TL','BL'),c('TR','BR'))

for (i in 2){ ## change the quad-pair HERE
  ## select which quads to compare (TL, TR, BL, BR)
  quads <- which(quad.labels$Quadrant %in% quads.pairs[[i]])
  
  ## Wilcox test with Bonferroni correction
  result0 <- as.data.frame(apply(ccle.sorted[,quads],1,function(x) wilcox.test(x ~ quad.labels$Quadrant[quads], exact=FALSE)$p.value))
  quadrant.sig.genes <- sum(result0[,1] <= 0.05)
  gene.order0 <- order(result0[,1])[1:quadrant.sig.genes] ## select lowest p-values in increasing order
  
  # write.table(t(c(paste(quads.pairs[[i]], collapse = ' & '),sps.gene.labels[gene.order0,1])), './quadrants-significant-genes-list.txt', append=T, col.names=F, row.names=F)
}
  
## heatmap for all SPS genes
pheatmap(ccle.sorted[gene.order0,quads], labels_row = sps.gene.labels[gene.order0,1],
         labels_col = quad.labels[quads,2],
         annotation_col = quad.labels[quads,1,drop=F], 
         cluster_cols = F, cluster_rows = T,
         clustering_distance_rows = 'euclidean')


######################## Create File with Quadrant Labels ############################


# new.row <- rep(NA,55)
# ccle.sub2 <- ccle.sub
# 
# for (i in 1:nrow(quad.labels)){
#   new.row[[i+2]] <- as.character(quad.labels$Quadrant[which(rownames(quad.labels) == colnames(ccle.sub)[[i+2]])])
# }
# 
# ccle.sub2 <- rbind(ccle.sub2, new.row)
# 
# write.csv(ccle.sub2, 'ccle-breast-cancer-with-quadrant-labels.csv')


####################################################


# # K-W test with Bonferroni correction
# res <- as.data.frame(apply(ccle.test,1,function(x) kruskal.test(x ~ quad.labels[,1])$p.value)*nrow(ccle.test))
# res.order <- order(res[,1])
# 
# ccle.test2 <- ccle.test[which(res < 0.05),]
# 
# # heatmap for only significant K-W genes
# pheatmap(ccle.test[which(res < 0.05),], labels_row = gene.labels[which(res < 0.05),1],
#          labels_col = metadata$CCLE_ID[filtered.ids-2][order.ids],
#          annotation_col = quad.labels, cluster_cols = F, cluster_rows = T,
#          clustering_distance_rows = 'euclidean',
#          display_numbers = F)
# 
# # Wilcox test between quadrant 3 and 4
# res3 <- data.frame(apply(ccle.test,1,function(x) wilcox.test(x = x[27:42], y = x[43:51], alternative='less',exact=FALSE)$p.value))
# res3.order <- order(res3[,1])
# 
# # heatmap for only top 10 significant genes between quad 3 and 4
# pheatmap(ccle.test[res3.order,27:51][1:10,], labels_row = gene.labels[res3.order,1][1:10], 
#          labels_col = metadata$CCLE_ID[filtered.ids-2][order.ids][27:51],
#          annotation_col = quad.labels[27:51,,drop=F], cluster_cols = F, cluster_rows = F,
#          clustering_distance_rows = 'euclidean',
#          display_numbers = F)


######################## Plot PCA with PC loadings ############################


sub.pca.var <- get_pca_var(subset.pca)
sub.contrib <- as.data.frame(sub.pca.var$contrib)
sub.contrib$entrez_id <- ccle.sub$entrez_id
sub.contrib$gene_name <- ccle.sub$gene_name
sub.correl <- as.data.frame(sub.pca.var$cor)

topContributors <- function(df, cor, pc, n) {
  col <- paste('Dim.',pc,sep='')
  corr <- cor[order(df[,pc], decreasing=TRUE),col]
  df <- df[order(df[,pc], decreasing=TRUE),]
  res <- df[1:n,c(col,'entrez_id')]
  colnames(res) <- c(paste('PC',pc,sep=''),'Entrez')
  res$gene_name <- df$gene_name[1:n]
  res$correlation <- corr[1:n]
  return(list(res, rownames(res)))
}

plotStrip <- function(res, df0, label.ids, y.val, ylimits) {
  df1 <- as.data.frame(t(df0[rownames(res),3:(ncol(df0)-2)]))
  
  df2 <- as.data.frame(df0[rownames(res),ncol(df0)])
  colnames(df2) <- 'variable'
  df2$ycoord <- apply(df1, 2, max) + 0.8
  df2$label <- as.factor(round(res$correlation, digits=2))
  
  df1$tissue <- rownames(df1)
  df1$tissue <- gsub('\\..+$', '', df1$tissue)
  # ids <- which(df1$tissue == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
  df1$label <- 'Above'
  df1$label[label.ids] <- 'Below'
  colnames(df1)[1:(ncol(df1)-2)] <- df0[rownames(res),ncol(df0)]
  
  test <- melt(df1)
  
  # ggplot(test, aes_string(x='variable',y='value',color='label')) + 
  #   geom_boxplot(outlier.size=1) + 
  #   # geom_text(data = df2, aes_string(x = 'variable', y = y.val, label = 'label'), inherit.aes = FALSE, angle=45, size=3.5) +
  #   xlab('Gene Name') + ylab('log (Expression)') + labs(title = paste('Top',nrow(res),'genes in',names(res)[1]),
  #                                                       caption = "Numbers indicate correlation with PC") +
  #   theme(legend.position='none',
  #         axis.text.x=element_text(angle=45,hjust=1)) +
  #   coord_cartesian(ylim=ylimits)
  
  ggplot(test, aes_string(x='variable',y='value',color='label')) +
    geom_boxplot(outlier.size=1) +
    # geom_text(data = df2, aes_string(x = 'variable', y = y.val, label = 'label'), inherit.aes = FALSE, angle=45, size=3.5) +
    xlab('Gene Name') + ylab('log (Expression)') + labs(title = paste('Top',nrow(res),'genes in',names(res)[1]),
                                                        caption = "Numbers indicate correlation with PC") +
    theme(legend.position='none',
          axis.text.x=element_text(angle=45,hjust=1)) +
    coord_cartesian(ylim=ylimits)
}

result1 <- topContributors(sub.contrib, cor=sub.correl, pc=1, n=20)
result2 <- topContributors(sub.contrib, cor=sub.correl, pc=2, n=20)

grid.arrange(
  plotStrip(result1[[1]], ccle.processed[,c(1:2,filtered.ids,ncol(ccle.processed)-1,ncol(ccle.processed)),drop=F], list[[7]], 7, c(0.5,6)), 
  plotStrip(result2[[1]], ccle.processed[,c(1:2,filtered.ids,ncol(ccle.processed)-1,ncol(ccle.processed)),drop=F], list[[8]], 9, c(-5,9)), ncol=2)


pc1.top <- result1[[1]]
pc1.top$label <- 'FALSE'
levels(pc1.top$label) <- c('TRUE','FALSE')
pc2.top <- result2[[1]]
pc2.top$label <- as.factor(pc2.top$correlation <= 0)
pc2.top$label <- relevel(pc2.top$label,'TRUE')

top1 <- ggplot(pc1.top, aes(x=reorder(gene_name, -correlation),y=correlation)) +
  geom_bar(stat='identity', fill='#619CFF',width = 0.8) +
  theme(legend.position='none',
        axis.text.x=element_text(angle=45,hjust=1)) +
  coord_cartesian(ylim=c(0.3,0.9)) +
  xlab('Gene Name') + ylab('abs (Correlation)') + labs(title = 'Top 20 genes in PC1')

top2 <- ggplot(pc2.top, aes(x=reorder(gene_name, -abs(correlation)),y=abs(correlation))) +
  geom_bar(stat='identity', aes(fill=label,width = 0.8)) +
  theme(legend.position='none',
        axis.text.x=element_text(angle=45,hjust=1)) +
  coord_cartesian(ylim=c(0.3,0.9)) +
  scale_fill_manual(values=c("#F8766D", "#619CFF")) +
  xlab('Gene Name') + ylab('') + labs(title = 'Top 20 genes in PC2') +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

grid.arrange(top1,top2,ncol=2)
