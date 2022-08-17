library('biomaRt')
library('dplyr')
library(ggplot2)
library(factoextra)
library(gridExtra)
library(reshape2)
library(pheatmap)

set.seed(500)

ccle <- read.table('CCLE_RNAseq_rsem_genes_tpm_20180929.txt', header=T, check.names=F)

ccle.id <- data.frame(colnames(ccle[3:1021]))
colnames(ccle.id) <- 'CCLE_ID'
ccle.id$model_name <- gsub('_.+$', '', ccle.id$CCLE_ID)

# metadata <- read.csv('./model_list_20220124.csv')
# metadata$model_name <- toupper(gsub("[^[:alnum:][:space:]]","",metadata$model_name)) ## retain alphanumeric characters only
# metadata <- left_join(ccle.id, metadata, by='CCLE_ID')
# write.table(metadata,'./CCLE_metadata.txt', row.names=FALSE)

ccle[,3:1021] <- lapply(ccle[,3:1021], function(x) as.numeric(x))

ccle <- ccle[which(rowSums(ccle[,3:1021]) != 0),]
ccle[,3:ncol(ccle)] <- log(ccle[,3:ncol(ccle)]+0.000001)# log-normalising the data

# remove version number from gene_id
ccle$gene_id <- gsub('\\..+$', '', ccle$gene_id)

# remove cell line # from column headings
# colnames(ccle)[3:ncol(ccle)] <- gsub('^[^_]*_','',colnames(ccle)[3:ncol(ccle)])

# using biomaRt to transform Ensembl IDs to Entrez IDs and Gene Name
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_id_list <- getBM(filters="ensembl_gene_id",
                attributes=c("ensembl_gene_id","entrezgene_id",'external_gene_name'),
                values=ccle$gene_id,
                mart=mart)

# appending Gene Name and Entrez ID to CCLE dataset
colnames(gene_id_list) <- c('gene_id','entrez_id','gene_name')
gene_id_list_dedup <- gene_id_list[!duplicated(gene_id_list$gene_id),]
ccle <- left_join(ccle,gene_id_list_dedup,by='gene_id')

# write.table(ccle,'./CCLE_processed_dedup.txt', row.names = FALSE)

# 26 different types of tissue - sorted in decreasing order
colnames(ccle)[3:1021] <- gsub('^[^_]*_', '', colnames(ccle)[3:1021])
cell.line.labels <- as.data.frame(sort(table(colnames(ccle)[3:1021]),decreasing=TRUE))

################################

# random.ids <- vector(mode = "list", length = 1000000)
# 
# for (j in 1:1000000) {
#   random <- sample(rownames(ccle),81)
#   random.ids[[j]] <- random
# }
# 
# write.table(random.ids, './random-ids-1000000.txt')
# 
# set.seed(500)
# 
# row.names <- rownames(ccle)
# 
# for (j in 1:1000) {
#   print(j)
#   
#   random.ids <- sample(row.names,81)
#   # random.ids <- rownames(ccle.sps)
#   
#   ccle.subset <- ccle[random.ids,3:1021]
#   subset.pca <- prcomp(t(ccle.subset), scale = TRUE)
#   
#   first.var <- fviz_eig(subset.pca)$data[1,2]
#   
#   write(first.var, './first-variances-test.txt', append = TRUE)
# }
# 
# test <- read.table('./first-variances-50000.txt')
# 
# ggplot(as.data.frame(first.var), aes(x=first.var)) + geom_histogram()

################################

# picking out SPS rows from CCLE dataset
# sps.genes <- read.table('SPS_entrez.txt', header=TRUE)
# ccle.sps <- ccle[ccle$entrez_id %in% sps.genes$Entrez,]

# PCA with scaling and centering
# sps.pca <- prcomp(t(ccle.sps[,3:1021]), scale = TRUE)
# ccle.pca <- prcomp(t(ccle[,3:1021]), scale = TRUE)

# plotPCA <- function(pca, title, topn = 1:10){
#   pca.plot <- as.data.frame(pca$x[,c('PC1','PC2')])
#   pca.plot$labels <- gsub('\\..+$', '', rownames(pca.plot))
#   pca.plot$labels2[pca.plot$labels != 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'] <- 'Solid'
#   pca.plot$labels2[pca.plot$labels == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'] <- 'Non-Solid'
#   
#   pca.plot$status <- metadata$tissue_status
#   pca.plot <- pca.plot[which(pca.plot$status %in% c('Tumour','Normal','Metastasis','Unknown','NA') &
#                                pca.plot$labels %in% cell.line.labels$Var1[topn]),]
#   
#   plot1 <- ggplot(pca.plot,aes(x=PC1,y=PC2)) + 
#     geom_point(aes(colour = labels),size=1) + 
#     ggtitle(title) +
#     theme(legend.text = element_text(size=5),
#           legend.position = 'right') +
#     geom_hline(yintercept = median(pca.plot$PC2), col = "red") +
#     geom_vline(xintercept = median(pca.plot$PC1), col = 'red')
#   
#   plot2 <- ggplot(pca.plot,aes(x=PC1,y=PC2)) + 
#     geom_point(aes(colour = status),size=1) +
#     ggtitle(title) +
#     theme(legend.text = element_text(size=15),
#           legend.position = 'right') +
#     geom_hline(yintercept = median(pca.plot$PC2), col = "red") +
#     geom_vline(xintercept = median(pca.plot$PC1), col = 'red')
#   
#   quadrants_cw <- c(sum(pca.plot$PC1 <= median(pca.plot$PC1) & pca.plot$PC2 > median(pca.plot$PC2)),
#                  sum(pca.plot$PC1 > median(pca.plot$PC1) & pca.plot$PC2 > median(pca.plot$PC2)),
#                  sum(pca.plot$PC1 > median(pca.plot$PC1) & pca.plot$PC2 <= median(pca.plot$PC2)),
#                  sum(pca.plot$PC1 <= median(pca.plot$PC1) & pca.plot$PC2 <= median(pca.plot$PC2)))
#   
#   percentages <- c(sum(pca.plot$PC1 <= median(pca.plot$PC1) & pca.plot$PC2 > median(pca.plot$PC2) & pca.plot$status %in% c('Metastasis')),
#                    sum(pca.plot$PC1 > median(pca.plot$PC1) & pca.plot$PC2 > median(pca.plot$PC2) & pca.plot$status %in% c('Metastasis')),
#                    sum(pca.plot$PC1 > median(pca.plot$PC1) & pca.plot$PC2 <= median(pca.plot$PC2) & pca.plot$status %in% c('Metastasis')),
#                    sum(pca.plot$PC1 <= median(pca.plot$PC1) & pca.plot$PC2 <= median(pca.plot$PC2) & pca.plot$status %in% c('Metastasis')))
#   
#   return(list(plot1,plot2, quadrants_cw))
# }
# 
# list <- plotPCA(sps.pca, "SPS PCA - 81 Genes", 1:26)
# p1 <- list[[1]]
# p2 <- list[[2]]
# percent <- list[[3]]
# 
# list <- plotPCA(ccle.pca, "CCLE PCA - 46674 Genes", 1:26)
# p3 <- list[[1]]
# p4 <- list[[2]]
# 
# f1 <- fviz_eig(sps.pca)
# f2 <- fviz_eig(ccle.pca)
# 
# grid.arrange(p1,p3,p2,p4,ncol=2)

#######################################################

# sps.pca.var <- get_pca_var(sps.pca)
# sps.contrib <- as.data.frame(sps.pca.var$contrib)
# sps.contrib$entrez_id <- ccle.sps$entrez_id
# sps.contrib$gene_name <- ccle.sps$gene_name
# sps.correl <- as.data.frame(sps.pca.var$cor)
# 
# ccle.pca.var <- get_pca_var(ccle.pca)
# ccle.contrib <- as.data.frame(ccle.pca.var$contrib)
# ccle.contrib$entrez_id <- ccle$entrez_id
# ccle.contrib$gene_name <- ccle$gene_name
# ccle.correl <- as.data.frame(ccle.pca.var$cor)
# 
# topContributors <- function(df, cor, pc, n) {
#   col <- paste('Dim.',pc,sep='')
#   corr <- cor[order(df[,pc], decreasing=TRUE),col]
#   df <- df[order(df[,pc], decreasing=TRUE),]
#   res <- df[1:n,c(col,'entrez_id')]
#   colnames(res) <- c(paste('PC',pc,sep=''),'Entrez')
#   res$gene_name <- df$gene_name[1:n]
#   res$correlation <- corr[1:n]
#   return(list(res, rownames(res)))
# }
# 
# plotStrip <- function(res, df0) {
#   df1 <- as.data.frame(t(df0[rownames(res),3:1021]))
#   
#   df2 <- as.data.frame(df0[rownames(res),1023])
#   colnames(df2) <- 'variable'
#   df2$ycoord <- apply(df1, 2, max) + 0.5
#   df2$label <- as.factor(round(res$correlation, digits=2))
#   
#   df1$tissue <- rownames(df1)
#   df1$tissue <- gsub('\\..+$', '', df1$tissue)
#   ids <- which(df1$tissue == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
#   df1$label <- 'Solid'
#   df1$label[ids] <- 'Non-Solid'
#   colnames(df1)[1:(ncol(df1)-2)] <- df0[rownames(res),1023]
#   
#   test <- melt(df1)
#   
#   ggplot(test, aes_string(x='variable',y='value',color='label')) + 
#     geom_boxplot(outlier.size=0) + scale_x_discrete(guide = guide_axis(n.dodge=2)) +
#     geom_text(data = df2, aes_string(x = 'variable', y = 'ycoord', label = 'label'), inherit.aes = FALSE) +
#     xlab('Gene Name') + ylab('log (Expression)') + labs(title = paste('Top',nrow(res),'genes in',names(res)[1]),
#                                                         caption = "Numbers indicate correlation with PC")
# }
# 
# result <- topContributors(sps.contrib, cor=sps.correl, pc=1, n=20)
# plotStrip(result[[1]], ccle)

#######################################################

# ccle.sps2 <- ccle.sps[result[[2]],3:1021]
# rownames(ccle.sps2) <- ccle.sps[result[[2]],1023]
# ccle.sps2 <- t(scale(t(ccle.sps2)))
# ccle.sps2 <- ccle.sps2[,order(colnames(ccle.sps2))]
# ids <- which(startsWith(colnames(ccle.sps2),'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'))
# 
# mat_col <- data.frame(matrix(NA, nrow=ncol(ccle.sps2), ncol=2))
# mat_col[ids,1] <- 'NON-SOLID'
# mat_col[setdiff(1:ncol(ccle.sps2), ids),1] <- 'SOLID'
# mat_col[,2] <- gsub('\\..+$', '', colnames(ccle.sps2))
# rownames(mat_col) <- colnames(ccle.sps2)
# colnames(mat_col) <- c('Label1','Label2')
# 
# ccle.sps2 <- ccle.sps2[,order(mat_col[,1])]
# 
# ccle.sps2[which(ccle.sps2 > 4)] <- 4
# ccle.sps2[which(ccle.sps2 < -4)] <- -4
# 
# pheatmap(ccle.sps2, annotation_col = mat_col,
#                 show_colnames = FALSE, show_rownames = TRUE,
#                 cluster_cols = FALSE)

#######################################################

# test <- as.data.frame(sps.pca$x[,c(1:3)])
# test$tissue.type <- gsub('\\..+$', '', rownames(test))
# test$ccle.id <- metadata$CCLE_ID
# test$status <- metadata$tissue_status
# test$status[which(is.na(test$status))] <- 'Unknown'
# 
# test <- test[order(rownames(test)),]
# tissue.ids <- which(test$tissue.type %in% cell.line.labels$Var1[1:20])
# ids <- which(startsWith(rownames(test),'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'))
# 
# 
# mat_col2 <- data.frame(matrix(NA, nrow=nrow(test), ncol=3))
# mat_col2[ids,1] <- 'NON-SOLID'
# mat_col2[setdiff(1:nrow(mat_col2), ids),1] <- 'SOLID'
# mat_col2[,2] <- gsub('\\..+$', '', rownames(test))
# mat_col2[,3] <- test$status
# rownames(mat_col2) <- rownames(test)
# 
# pheatmap(test[tissue.ids,1:3], annotation_row = mat_col2[tissue.ids,],
#          show_rownames=FALSE)

# cut <- cutree(res$tree_row,3)
# 
# group1 <- test[tissue.ids[which(cut==1)],]
# group2 <- test[tissue.ids[which(cut==2)],]
# group3 <- test[tissue.ids[which(cut==3)],]


#######################################################

# library(genefilter)
# 
# stat.ids <- which(metadata$tissue_status %in% c('Tumour','Metastasis'))
# status.labels <- metadata[stat.ids,][,'tissue_status']
# 
# sub.ccle <- ccle[,3:1021]
# stat.factor <- as.factor(status.labels)
# 
# t.test <- rowttests(as.matrix(sub.ccle[,stat.ids]), stat.factor)
# 
# sig.ids <- which(t.test$p.value < 0.00000000001)
# sps.sig.ids <- sig.ids[which(ccle$entrez_id[sig.ids] %in% sps.genes$Entrez)]
# 
# stripchart(as.numeric(ccle[36,stat.ids])~stat.factor, method="jitter")
