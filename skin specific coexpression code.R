library(data.table)
library(tidyverse)
library(preprocessCore)

# expression data for all tissues
data <- fread(cmd="zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f3-")
colnames(data) <- scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR==3' | cut -f3-"),what="")
rownames(data) <- scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f1"),what="")
final_data <- as.matrix(data) 

## Sample Attributes table
table1 <- read.table("GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt",header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 

## Subject Phenotypes table
table2 <- read.table("GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt",header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 

## How to get FPKM expression data for each tissue and get the largest subset's FPKM expression data for some tissue
Exposed_Skinid <- rownames(table1)[table1$SMTSD == "Skin - Sun Exposed (Lower leg)"]
Exposed_SkinData <- final_data[,colnames(final_data)%in%Exposed_Skinid]
New_Exposed_SkinData <- Exposed_SkinData[rowSums(Exposed_SkinData!=0)>= (473*0.1),]

Adipose_id <- rownames(table1)[table1$SMTS == "Adipose Tissue"]
Adipose_Data <- final_data[,colnames(final_data)%in%Adipose_id]
New_AdiposeData <- Adipose_Data[rownames(New_Exposed_SkinData),]

Subcutaneous_id <- rownames(table1)[table1$SMTSD == "Adipose - Subcutaneous"]
Subcutaneous_Data <- final_data[,colnames(final_data)%in%Subcutaneous_id]
New_Subcutaneous_Data <- Subcutaneous_Data[rownames(New_Exposed_SkinData),]

Adrenal_id = rownames(table1)[table1$SMTS == "Adrenal Gland"]
Adrenal_data = final_data[,colnames(final_data)%in%Adrenal_id]
New_AdrenalData = Adrenal_data[rownames(New_Exposed_SkinData),]


## How to get normalized expression data for each tissue after quantiles normalization and inverse normalization
tempdata_quantile <- normalize.quantiles(New_Exposed_SkinData)
rownames(tempdata_quantile) <- rownames(New_Exposed_SkinData)
colnames(tempdata_quantile) <- colnames(New_Exposed_SkinData)

inverse_Exposed_SkinData <- t(apply(tempdata_quantile,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))}))

## How to get spearman correlation matrix of each tissue's normalized expression data
inverse_Exposed_SkinCor <- cor(t(inverse_Exposed_SkinData), method = "spearman")
inverse_AdiposeCor <- cor(t(inverse_AdiposeData), method = "spearman")
inverse_AdrenalCor <- cor(t(inverse_AdrenalData), method = "spearman")

## How to do correlation comparison to find positive correlated gene pairs that're only significant in sun-exposed skin
inverse_Exposed_SkinCor[inverse_Exposed_SkinCor < 0 ] <- 0
inverse_Exposed_SkinCor[is.na(inverse_Exposed_SkinCor)] <- 0
diag(inverse_Exposed_SkinCor) <- 0
inverse_AdiposeCor[is.na(inverse_AdiposeCor)] <- 0

compare1 <- (inverse_Exposed_SkinCor > inverse_AdiposeCor) & (inverse_Exposed_SkinCor > 0)
inverse_Exposed_SkinCor <- ifelse(compare1, inverse_Exposed_SkinCor, 0)
compare2 <- (inverse_Exposed_SkinCor > inverse_AdrenalCor) & (inverse_Exposed_SkinCor > 0)
inverse_Exposed_SkinCor <- ifelse(compare2, inverse_Exposed_SkinCor, 0) # keep doing this for all other tissues 

compare29 <- (inverse_Exposed_SkinCor > inverse_IntestineCor) & (inverse_Exposed_SkinCor > 0)
inverse_final_compared <- ifelse(compare29, inverse_Exposed_SkinCor, 0)


## linear regression formula
index_Compared <- which(inverse_final_compared>0,arr.ind=T)
index_Compared <- index_Compared[which(index_Compared[,1]>index_Compared[,2]),]
index_Compared <- cbind(colnames(inverse_final_compared)[index_Compared[,1]],
                        colnames(inverse_final_compared)[index_Compared[,2]])

form <- list()
for(i in 1:2042499){ 
  form[[i]] = as.formula(paste(index_Compared[i,1], 
                               paste(c(index_Compared[i,2],"SEX","AGE","BMI"),collapse = "+"), sep = "~"))
  }


## linear regression source data
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
regression_source3 <- merge(table2_simplified,inverse_MuscleData,by="SubjectID") # keep doing this for all tissue's largest subset


## perform linear regression, get beta and se values

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
j = args[1]
form = readRDS("form.rds")
regression_sourcej = readRDS(paste0("regression_source", j, ".rds"))
summary_j = list()
for(i in 1:length(form)){
  summary_j[[i]] = summary(lm(form[[i]],data=regression_sourcej))
}
saveRDS(summary_j, paste0("summary_", j ,".rds"))



#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
j = args[1]
summary_j = readRDS(paste0("summary_", j, ".rds"))
beta_j = c()
for(i in 1:length(summary_j)){
  beta_j[i] = summary_j[[i]]$coefficients[2,1]
}
saveRDS(beta_j, paste0("beta_", j ,".rds"))



#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
j = args[1]
summary_j = readRDS(paste0("summary_", j, ".rds"))
se_j = c()
for(i in 1:length(summary_j)){
  se_j[i] = summary_j[[i]]$coefficients[2,2]
}
saveRDS(se_j, paste0("se_", j ,".rds"))


## comparing summary statistics

one_2_beta = beta_1 - beta_2
one_2_se = sqrt((se_1)^2 + (se_2)^2)
one_2_z = one_2_beta / one_2_se
one_2_pnorm = pnorm(abs(one_2_z), lower.tail=F)*2
one_2_fdr = p.adjust(one_2_pnorm, method="fdr")

one_3_beta = beta_1 - beta_3
one_3_se = sqrt((se_1)^2 + (se_3)^2)
one_3_z = one_3_beta / one_3_se
one_3_pnorm = pnorm(abs(one_3_z), lower.tail=F)*2
one_3_fdr = p.adjust(one_3_pnorm, method="fdr")

c1 = ((one_2_beta > 0) & (one_2_fdr < 0.1))
c2 = c1 & ((one_3_beta > 0) & (one_3_fdr < 0.1))
c3 = c2 & ((one_4_beta > 0) & (one_4_fdr < 0.1)) ## keep doing this for every beta, se
c29 = c28 & ((one_30_beta >0) & (one_30_fdr < 0.1))


## get the final significant protein-encoding gene pairs
gene_symbol =read.table(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3{print $1\"\t\"$2}'"),
                        stringsAsFactors=F)

index_Compared_1 = subset(index_Compared, c29)
index_Compared_1 = as_tibble(index_Compared_1)
colnames(index_Compared_1) = c("gene1", "gene2")
index_Compared_1 = index_Compared_1 %>% left_join(gene_symbol, by = c("gene1" = "V1"))
colnames(index_Compared_1)[3] = "gene1_symbol"
index_Compared_1 = index_Compared_1 %>% left_join(gene_symbol, by = c("gene2" = "V1"))
colnames(index_Compared_1)[4] = "gene2_symbol"
protein_symbol = protein_symbol %>% select(gene, Category)
index_Compared_1 = index_Compared_1 %>% left_join(protein_symbol, by = c("gene1_symbol" = "gene"))
colnames(index_Compared_1)[5] = "gene1_category"
index_Compared_1 = index_Compared_1 %>% left_join(protein_symbol, by = c("gene2_symbol" = "gene"))
colnames(index_Compared_1)[6] = "gene2_category"
index_Compared_1 = filter(index_Compared_1, gene1_category == "Protein_coding" & gene2_category == "Protein_coding")
saveRDS(index_Compared_1, "~/final_left_gene_pairs.rds")


## get the final significant genes and their number of connections
final_left_genes = final_left_gene_pairs[,1:2] %>% 
                   pivot_longer(gene1:gene2, values_to = "gene") %>% 
                   select(-1) %>% 
                   count(gene) %>% 
                   arrange(desc(n))












