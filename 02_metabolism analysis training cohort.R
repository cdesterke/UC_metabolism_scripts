meta<-read.csv("metabolism.csv",h=T)

load("bm110.rda")

library(dplyr)
meta%>%dplyr::rename(gene.mm="mm.gene")->meta
bm110%>%inner_join(meta,by="gene.mm")->metahs
metahs%>%dplyr::rename(gene="hs.gene")->metahs
library(transpipe15)
##df
load("data.rda")
df$ID_REF<-NULL

metahs%>%select(hs.id,gene)%>%
	inner_join(df,by="gene",relationship = "many-to-many")->meta_uc

meta_uc$hs.id<-NULL

ok<-filtermatrix(meta_uc)

pheno<-read.csv("pheno.csv",h=T,row.names=1)
all(row.names(pheno)==colnames(ok))
pcatrans(ok,pheno,group="group",pal="Dark2",alpha=1,names=F)


pcatransellipses(ok,pheno,group="group2",pal="Dark2",alpha=1,names=F,x=1,y=2,level=0.65)


res<-deg(ok,pheno$group2,control="control")

sig<-filtresig(res)
dim(sig)
vollimma(res,nb=500,fc=0.5,p=0.05,size=3,alpha=1)


pheno%>%select(group,group2,gender,disease_extension)->annot
bestheat(ok, annot, scale = "none", font = 10, rownames = F) 

process<-reducedf(sig,ok,n=145)

pcatransellipses(process,pheno,group="group2",pal="Dark2",alpha=1,names=F,x=1,y=2,level=0.65)

write.table(sig,file="limmasig.tsv",row.names=T,sep="\t")

trans<-as.data.frame(t(process))

save(trans,file="trans145.rda")


