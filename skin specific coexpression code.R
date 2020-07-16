final_data = readRDS("final_data.rds") ## expression data for all tissues
table1 = readRDS("table1.rds")  ## Sample Attributes table
Exposed_Skinid = rownames(table1)[table1$SMTSD == "Skin - Sun Exposed (Lower leg)"]
Exposed_SkinData = final_data[,colnames(final_data)%in%Exposed_Skinid]
New_Exposed_SkinData = Exposed_SkinData[rowSums(Exposed_SkinData!=0)>= (473*0.1),]