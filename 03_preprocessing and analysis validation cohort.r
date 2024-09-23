platform<-read.csv("plateform.csv",h=T)

data<-read.table("GSE.txt",h=T,sep="\t")

library(dplyr)

platform%>%inner_join(data,by="ID_REF")->all

all<-all[complete.cases(all$gene),]

head(all)

all$ID_REF<-NULL

library(transpipe14)

ok<-filtermatrix(all)

save(ok,file="matrix.rda")

pheno<-read.csv("pheno.csv",h=T,row.names=1)

pcatrans(ok,pheno,group="group",pal="Set1",alpha=0.7,names=F)

ok[ok=="null"]<-NA

ok %>% mutate_if(is.character, as.numeric)->ok

all(colnames(ok)==row.names(pheno))
trans<-as.data.frame(t(ok))
trans$group<-pheno$group
library(randomForest)
set.seed(1357)

sig<-read.table("limmasig.tsv",h=T,sep="\t")
sig$gene<-row.names(sig)
sig%>%filter(logFC>=1)->pos
vector<-pos$gene
inter<-intersect(colnames(trans),vector)
trans%>%select(group,all_of(inter))->small
summary(small)

## select columns without NA
df<-small[ , colSums(is.na(small)) == 0]


## save data for machine learning
save(df,file="df_for_ml.rda")






