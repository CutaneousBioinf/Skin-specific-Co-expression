library(data.table)
library(tidyverse)
library(preprocessCore)
library(NbClust)
library(umap)
library(cowplot)
library(reshape2)


# expression data for all tissues
data <- fread(cmd="zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f3-")
colnames(data) <- scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR==3' | cut -f3-"),what="")
rownames(data) <- scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f1"),what="")
final_data <- as.matrix(data) 


## Sample Attributes table
table1 <- read.table("GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt",
                     header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 


## Subject Phenotypes table
table2 <- read.table("GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt",
                     header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 


## How to get FPKM expression data for each tissue and get the largest subset's FPKM expression data for some tissue
Exposed_Skinid <- rownames(table1)[table1$SMTSD == "Skin - Sun Exposed (Lower leg)"]
Exposed_SkinData <- final_data[,colnames(final_data)%in%Exposed_Skinid]
New_Exposed_SkinData <- Exposed_SkinData[rowSums(Exposed_SkinData!=0)>= (473*0.1),]

Adipose_id <- rownames(table1)[table1$SMTS == "Adipose Tissue"]
Adipose_Data <- final_data[,colnames(final_data)%in%Adipose_id]
New_AdiposeData <- Adipose_Data[rownames(New_Exposed_SkinData),]

Subcutaneous_id <- rownames(table1)[table1$SMTSD == "Adipose - Subcutaneous"]
Subcutaneous_Data <- final_data[,colnames(final_data)%in%Subcutaneous_id]
New_SubcutaneousData <- Subcutaneous_Data[rownames(New_Exposed_SkinData),]

Adrenal_id <- rownames(table1)[table1$SMTS == "Adrenal Gland"]
Adrenal_Data <- final_data[,colnames(final_data)%in%Adrenal_id]
New_AdrenalData <- Adrenal_Data[rownames(New_Exposed_SkinData),]

Muscle_id <- rownames(table1)[table1$SMTS == "Muscle"]
Muscle_Data <- final_data[,colnames(final_data)%in%Muscle_id]
New_MuscleData <- Muscle_Data[rownames(New_Exposed_SkinData),]

Tibial_id <- rownames(table1)[table1$SMTSD == "Artery - Tibial"]
Tibial_Data <- final_data[,colnames(final_data)%in%Tibial_id]
New_TibialData <- Tibial_Data[rownames(New_Exposed_SkinData),]


## How to get normalized expression data for each tissue after quantiles normalization and inverse normalization. 
## Using sun-exposed skin tissue as an example
tempdata_quantile <- normalize.quantiles(New_Exposed_SkinData)
rownames(tempdata_quantile) <- rownames(New_Exposed_SkinData)
colnames(tempdata_quantile) <- colnames(New_Exposed_SkinData)

inverse_Exposed_SkinData <- t(apply(tempdata_quantile,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))}))


## How to get spearman correlation matrix of each tissue's normalized expression data
inverse_Exposed_SkinCor <- cor(t(inverse_Exposed_SkinData), method = "spearman")
inverse_AdiposeCor <- cor(t(inverse_AdiposeData), method = "spearman")
inverse_AdrenalCor <- cor(t(inverse_AdrenalData), method = "spearman")
inverse_BladderCor <- cor(t(inverse_BladderData), method = "spearman")


## How to do correlation comparison to find positive correlated gene pairs that're only significant in sun-exposed skin
inverse_Exposed_SkinCor[inverse_Exposed_SkinCor < 0 ] <- 0
inverse_Exposed_SkinCor[is.na(inverse_Exposed_SkinCor)] <- 0
diag(inverse_Exposed_SkinCor) <- 0

inverse_AdiposeCor[is.na(inverse_AdiposeCor)] <- 0
compare1 <- (inverse_Exposed_SkinCor > inverse_AdiposeCor) & (inverse_Exposed_SkinCor > 0)
inverse_Exposed_SkinCor <- ifelse(compare1, inverse_Exposed_SkinCor, 0)

inverse_AdrenalCor[is.na(inverse_AdrenalCor)] <- 0
compare2 <- (inverse_Exposed_SkinCor > inverse_AdrenalCor) & (inverse_Exposed_SkinCor > 0)
inverse_Exposed_SkinCor <- ifelse(compare2, inverse_Exposed_SkinCor, 0) 

inverse_BladderCor[is.na(inverse_BladderCor)] <- 0
compare3 <- (inverse_Exposed_SkinCor > inverse_BladderCor) & (inverse_Exposed_SkinCor > 0)
inverse_Exposed_SkinCor <- ifelse(compare3, inverse_Exposed_SkinCor, 0) # keep doing this for all other tissues 

inverse_IntestineCor[is.na(inverse_IntestineCor)] <- 0
compare29 <- (inverse_Exposed_SkinCor > inverse_IntestineCor) & (inverse_Exposed_SkinCor > 0)
inverse_final_compared <- ifelse(compare29, inverse_Exposed_SkinCor, 0)


## get the linear regression formula
index_Compared <- which(inverse_final_compared>0,arr.ind=T)
index_Compared <- index_Compared[which(index_Compared[,1]>index_Compared[,2]),]
index_Compared <- cbind(colnames(inverse_final_compared)[index_Compared[,1]],
                        colnames(inverse_final_compared)[index_Compared[,2]])

form <- list()
for(i in 1:2042499){ 
  form[[i]] = as.formula(paste(index_Compared[i,1], 
                               paste(c(index_Compared[i,2],"SEX","AGE","BMI"),collapse = "+"), sep = "~"))
  }


## get the linear regression source data
table2_simplified <- table2[,c("SEX","AGE","BMI")]
table2_simplified <- table2_simplified[which(!apply(table2_simplified,1,function(x)any(is.na(x)))),]
table2_simplified$SubjectID <- rownames(table2_simplified)

inverse_Exposed_SkinData <- as.data.frame(t(readRDS("inverse_data/inverse_Exposed_SkinData.rds")))
inverse_Exposed_SkinData$SubjectID <- sapply(rownames(inverse_Exposed_SkinData),
                                             function(x)paste(unlist(strsplit(x,"-"))[1:2],collapse="-"))
regression_source1 <- merge(table2_simplified,inverse_Exposed_SkinData,by="SubjectID")


inverse_subcutaneous_data <- as.data.frame(t(readRDS("inverse_data/inverse_subcutaneous_data.rds")))
inverse_subcutaneous_data$SubjectID <- sapply(rownames(inverse_subcutaneous_data),
                                             function(x)paste(unlist(strsplit(x,"-"))[1:2],collapse="-"))
regression_source2 <- merge(table2_simplified,inverse_subcutaneous_data,by="SubjectID")


inverse_MuscleData <- as.data.frame(t(readRDS("inverse_data/inverse_MuscleData.rds")))
inverse_MuscleData$SubjectID <- sapply(rownames(inverse_MuscleData),
                                             function(x)paste(unlist(strsplit(x,"-"))[1:2],collapse="-"))
regression_source3 <- merge(table2_simplified,inverse_MuscleData,by="SubjectID") 


inverse_tibial_data <- as.data.frame(t(readRDS("inverse_data/inverse_tibial_data.rds")))
inverse_tibial_data$SubjectID <- sapply(rownames(inverse_tibial_data),
                                       function(x)paste(unlist(strsplit(x,"-"))[1:2],collapse="-"))
regression_source4 <- merge(table2_simplified,inverse_tibial_data,by="SubjectID") # keep doing this for all tissue's largest subset


## run the following R scripts to perform linear regression, get summary statistics

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
j <- args[1]
form <- readRDS("form.rds")
regression_sourcej <- readRDS(paste0("regression_source", j, ".rds"))
summary_j <- list()
for(i in 1:length(form)){
  summary_j[[i]] <- summary(lm(form[[i]],data=regression_sourcej))
}
saveRDS(summary_j, paste0("summary_", j ,".rds"))



#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
j <- args[1]
summary_j <- readRDS(paste0("summary_", j, ".rds"))
beta_j <- c()
for(i in 1:length(summary_j)){
  beta_j[i] <- summary_j[[i]]$coefficients[2,1]
}
saveRDS(beta_j, paste0("beta_", j ,".rds"))



#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
j <- args[1]
summary_j <- readRDS(paste0("summary_", j, ".rds"))
se_j <- c()
for(i in 1:length(summary_j)){
  se_j[i] <- summary_j[[i]]$coefficients[2,2]
}
saveRDS(se_j, paste0("se_", j ,".rds"))


## comparing summary statistics

one_2_beta <- beta_1 - beta_2
one_2_se <- sqrt((se_1)^2 + (se_2)^2)
one_2_z <- one_2_beta / one_2_se
one_2_pnorm <- pnorm(abs(one_2_z), lower.tail=F)*2
one_2_fdr <- p.adjust(one_2_pnorm, method="fdr")

one_3_beta <- beta_1 - beta_3
one_3_se <- sqrt((se_1)^2 + (se_3)^2)
one_3_z <- one_3_beta / one_3_se
one_3_pnorm <- pnorm(abs(one_3_z), lower.tail=F)*2
one_3_fdr <- p.adjust(one_3_pnorm, method="fdr")

one_4_beta <- beta_1 - beta_4
one_4_se <- sqrt((se_1)^2 + (se_4)^2)
one_4_z <- one_4_beta / one_4_se
one_4_pnorm <- pnorm(abs(one_4_z), lower.tail=F)*2
one_4_fdr <- p.adjust(one_4_pnorm, method="fdr")

c1 <- ((one_2_beta > 0) & (one_2_fdr < 0.1))
c2 <- c1 & ((one_3_beta > 0) & (one_3_fdr < 0.1))
c3 <- c2 & ((one_4_beta > 0) & (one_4_fdr < 0.1)) ## keep doing this for every beta, fdr

c29 <- c28 & ((one_30_beta >0) & (one_30_fdr < 0.1))


## get the final significant protein-encoding gene pairs
gene_symbol <- read.table(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3{print $1\"\t\"$2}'"),
                        stringsAsFactors=F)

index_Compared_1 <- subset(index_Compared, c29)
index_Compared_1 <- as_tibble(index_Compared_1)
colnames(index_Compared_1) <- c("gene1", "gene2")
index_Compared_1 <- index_Compared_1 %>% left_join(gene_symbol, by = c("gene1" = "V1"))
colnames(index_Compared_1)[3] <- "gene1_symbol"
index_Compared_1 <- index_Compared_1 %>% left_join(gene_symbol, by = c("gene2" = "V1"))
colnames(index_Compared_1)[4] <- "gene2_symbol"
protein_symbol <- protein_symbol %>% select(gene, Category)
index_Compared_1 <- index_Compared_1 %>% left_join(protein_symbol, by = c("gene1_symbol" = "gene"))
colnames(index_Compared_1)[5] <- "gene1_category"
index_Compared_1 <- index_Compared_1 %>% left_join(protein_symbol, by = c("gene2_symbol" = "gene"))
colnames(index_Compared_1)[6] <- "gene2_category"
index_Compared_1 <- filter(index_Compared_1, gene1_category == "Protein_coding" & gene2_category == "Protein_coding")
saveRDS(index_Compared_1, "~/final_left_gene_pairs.rds")


## get the final significant genes and their number of connections
final_left_genes <- final_left_gene_pairs[,1:2] %>% 
                   pivot_longer(gene1:gene2, values_to = "gene") %>% 
                   select(-1) %>% 
                   count(gene) %>% 
                   arrange(desc(n))


## perform PCA on significant genes' normalized expression data
pca_data_source <- subset(inverse_Exposed_SkinData,
                          rownames(inverse_Exposed_SkinData) %in% final_left_genes$gene)
skin_pca <- prcomp(pca_data_source, center = T, scale = T)


## silhouette plot to decide the optimal number of clusters
fviz_nbclust(skin_pca$x[,1:50], kmeans, method = "silhouette", k.max = 15) + theme_minimal() + ggtitle("The Silhouette Plot")


## perform UMPA on first 50 PCs and colored by the output of k-means clustering 
set.seed(200)
skin_umap = umap(skin_pca$x[,1:50])
km.out = kmeans(skin_pca$x[,1:50], centers = 6, nstart = 20)
plot(skin_umap$layout, col = (km.out+2), pch = 20)


## graph significant genes' connectivity
cor_data_source <- subset(inverse_Exposed_SkinData,
                          rownames(inverse_Exposed_SkinData) %in% final_left_genes$gene)
final_genes_cor <- cor(t(cor_data_source), method = "spearman")
diag(final_genes_cor) <- 0
final_genes_cor <- as.data.frame(final_genes_cor)
final_genes_connectivity <- final_genes_cor %>% 
                            transmute(connectivity = rowSums(., na.rm = T))
final_genes_connectivity <- final_genes_connectivity %>% 
                            arrange(connectivity)
final_genes_connectivity <- final_genes_connectivity %>% 
                            mutate(interval = cut_width(connectivity, width = 100, center = 50))
final_genes_connectivity <- final_genes_connectivity %>% 
                            mutate(interval_2 = cut_width(connectivity, width = 100, center = 50, labels = FALSE))
final_genes_connectivity <- final_genes_connectivity %>% 
                            mutate(interval_2 = (interval_2-1)*100 + 50)

final_genes_connectivity <- final_genes_connectivity %>% 
                           group_by(interval_2) %>% 
                           summarise(n = n()) %>% 
                           mutate(freq = n/sum(n))

final_genes_connectivity[,c(1,3)] <- log(final_genes_connectivity[,c(1,3)],10)

ggplot(data = final_genes_connectivity) + 
      geom_point(aes(x = interval_2, y = freq)) + 
      geom_smooth() + 
      theme_bw()


## graphing skin relative expression using FPKM expression values
New_SalivaryData <- readRDS("New_SalivaryData.rds")
salivary_mean <- rowMeans(New_SalivaryData, na.rm = T)
salivary_mean_2804 <- subset(salivary_mean, names(salivary_mean) %in% final_left_genes$gene)
salivary_mean_35385 <- subset(salivary_mean, !(names(salivary_mean) %in% final_left_genes$gene))

New_Exposed_SkinData <- readRDS("New_Exposed_SkinData.rds")
skin_mean <- rowMeans(New_Exposed_SkinData, na.rm = T)
skin_mean_2804 <- subset(skin_mean, names(skin_mean) %in% final_left_genes$gene)
skin_mean_35385 <- subset(skin_mean, !(names(skin_mean) %in% final_left_genes$gene)) ## do this for all other tissues

skin_relative_expression_2804 <- skin_mean_2804 * 29 / (adipose_mean_2804+adrenal_mean_2804+bladder_mean_2804+blood_mean_2804+
                                                       vessel_mean_2804+brain_mean_2804+breast_mean_2804+cervix_mean_2804+colon_mean_2804+
                                                       esophagus_mean_2804+fallopian_mean_2804+heart_mean_2804+kidney_mean_2804+liver_mean_2804+
                                                       lung_mean_2804+muscle_mean_2804+nerve_mean_2804+ovary_mean_2804+pancreas_mean_2804+
                                                       pituitary_mean_2804+prostate_mean_2804+spleen_mean_2804+stomach_mean_2804+testis_mean_2804+
                                                       thyroid_mean_2804+uterus_mean_2804+vagina_mean_2804+salivary_mean_2804+intestine_mean_2804)

skin_relative_expression_35385 <- skin_mean_35385 * 29 / (adipose_mean_35385+adrenal_mean_35385+bladder_mean_35385+blood_mean_35385+
                                                          vessel_mean_35385+brain_mean_35385+breast_mean_35385+cervix_mean_35385+colon_mean_35385+
                                                          esophagus_mean_35385+fallopian_mean_35385+heart_mean_35385+kidney_mean_35385+liver_mean_35385+
                                                          lung_mean_35385+muscle_mean_35385+nerve_mean_35385+ovary_mean_35385+pancreas_mean_35385+
                                                          pituitary_mean_35385+prostate_mean_35385+spleen_mean_35385+stomach_mean_35385+testis_mean_35385+
                                                          thyroid_mean_35385+uterus_mean_35385+vagina_mean_35385+salivary_mean_35385+intestine_mean_35385)

vec1 <- as.data.frame(skin_relative_expression_2804)
vec2 <- as.data.frame(skin_relative_expression_35385)

vec1 <- vec1 %>% mutate(gene_type = "2804 significant genes")
names(vec1)[1] <- "relative_expression"
vec2 <- vec2 %>% mutate(gene_type = "35385 non-significant genes")
names(vec2)[1] <- "relative_expression"
vec1_2 <- rbind(vec1, vec2)

main.plot <- ggplot(vec1_2,aes(x=relative_expression, fill=gene_type)) +  
             geom_density(alpha=0.25) + theme_bw() + scale_x_continuous(trans = 'log2') +           
             scale_y_continuous(limits = c(0,1))

inset.plot <- ggplot(vec1_2,aes(x=gene_type, y=relative_expression, fill=gene_type)) + 
              geom_boxplot() + theme_bw() + scale_y_continuous(trans = 'log2') + 
              theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), 
              legend.position = "none", axis.text.x=element_blank(),axis.ticks.x=element_blank()) + coord_flip()

ggdraw() + draw_plot(main.plot) + draw_plot(inset.plot, x = 0.27, y = 0.57, width = 0.4, height = 0.4)


## graphing skin specific expression using FPKM expression values

















