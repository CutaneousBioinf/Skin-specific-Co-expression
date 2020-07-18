library(data.table)
library(tidyverse)

data=fread(cmd="zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f3-")
colnames(data)=scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR==3' | cut -f3-"),what="")
rownames(data)=scan(pipe("zcat /net/dumbo/home/xwen/ncbi/dbGaP-9060/gtex_v7_data/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz | awk 'NR>3' | cut -f1"),what="")
final_data=as.matrix(data) # expression data for all tissues

## Sample Attributes table
table1 = read.table("GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt",header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 

## Subject Phenotypes
table2 = read.table("GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt",header=T,row.names=1,sep="\t",quote=NULL,comment.char="",stringsAsFactors=F) 

Exposed_Skinid = rownames(table1)[table1$SMTSD == "Skin - Sun Exposed (Lower leg)"]
Exposed_SkinData = final_data[,colnames(final_data)%in%Exposed_Skinid]
New_Exposed_SkinData = Exposed_SkinData[rowSums(Exposed_SkinData!=0)>= (473*0.1),]