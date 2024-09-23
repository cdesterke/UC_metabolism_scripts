
id<-pos$gene
beta<-pos$coef

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation




df%>%mutate(enet.score=(LIPG*3.98379269518892)+(PSAT1*3.69116126060057)+(PGM3*2.73987439434978)+
(CD38*2.28484414622611)+(BLVRA*1.98777489251261)+(CBR3*1.94276785064452)+
(NT5DC2*1.76136247817983)+(PHGDH*1.71099708709267)+(GPX7*1.57596477746821)+
(CASP1*1.56316238573594)+(ASRGL1*1.39975311631545)+(SOD3*1.2502038759225)+
(CHST2*0.964512346378638)+(CHST11*0.954009692731847)+(KYNU*0.941701068395825)+
(PLA2G7*0.915682328069718)+(SRM*0.871402818133839)+(PTGS2*0.797024504827082)+
(LPIN1*0.471364923366202)+(ME1*0.313674994225149)+(PTGDS*0.141580424565731)+
(ADA*0.132526286013327))->df

all<-paste(id,collapse="+")
res<-df


library(pROC)
library(Epi)
data<-as.data.frame(cbind(X,y))
data$y<-as.factor(data$y)
ROC(form = y~LIPG+PSAT1+PGM3+CD38+BLVRA+CBR3+NT5DC2+
PHGDH+GPX7+CASP1+ASRGL1+SOD3+CHST2+CHST11+
KYNU+PLA2G7+SRM+PTGS2+LPIN1+ME1+PTGDS+ADA , plot="ROC", data=data)



library(cutpointr)


cp <- cutpointr(df, enet.score,  group,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
df%>%mutate(enet.cat=ifelse(enet.score>=-9.37009,"HIGH","low"))->df
df$enet.cat<-as.factor(df$enet.cat)
df$enet.cat<-relevel(df$enet.cat,ref="low")

chisq.test(table(df$enet.cat,df$group))




library(vcd)

struct <- structable(~ enet.cat+group, data = df)
mosaic(struct, , direction = "v", pop = FALSE,colorize = T, shade = TRUE)
     #  gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)


save(df,file="enet_score.rda")

df$meta.cat<-as.factor(df$meta.cat)
df$meta.cat<-relevel(df$meta.cat,ref="low")


library(transpipe15)


data<-as.data.frame(X)

data%>%select(all_of(id))->data

trans<-as.data.frame(t(data))

df%>%select(group,enet.cat)->pheno

all(colnames(trans)==row.names(pheno))

bestheat(trans,scale="row",pheno,rownames=T,font=12)



library(ggplot2)
library(ggbeeswarm)

ggplot(df,aes(group,enet.score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Dark2")+
  geom_point(aes(fill=factor(group),size=1),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "right")+xlab("group")+ggtitle("")

groupe1 <- df$enet.score[df$group == "Normal"]
groupe2 <- df$enet.score[df$group == "UC"]

# Effectuer le test t pour échantillons indépendants
resultat <- t.test(groupe1, groupe2)
print(resultat)


ROC(form = group~enet.score , plot="ROC", data=df)


