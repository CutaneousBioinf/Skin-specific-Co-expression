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















