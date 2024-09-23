library(GEOquery)

# download data GEO GSE38713
gse <- getGEO("GSE38713", GSEMatrix = TRUE)

# if dataset has several data you can access like
set <- gse[[1]]

# extract expression
data<-exprs(set)

# extract annotation
features<-fData(set)
library(dplyr)
features%>%select(1,11)->small
colnames(small)<-c("ID_REF","gene")

library(stringr)
infos<-as.data.frame(str_split_fixed(small$gene," /// ",2))
small$gene<-infos$V1

small[small==""]<-NA

# extract phenotype
annot <- pData(set)
head(data)
small$REF_ID<-NULL

all<-merge(small,data,by="row.names")
dim(all)

df<-all[complete.cases(all$gene),]
df$Row.names<-NULL

save(data,file="data.rda")

write.table(annot,file="phenotype.tsv",row.names=T,sep="\t")
