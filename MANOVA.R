#install.packages("car")
setwd("~/Documents/ASF/")
list.files()
raw_L6<-read.csv("table_mc3566_sorted_L6.txt",sep="\t",header=T,dec=".",comment.char="#")

#identify the ASF taxa.
#transpose the data fram
#attach sample group


#10/08/2015
#Ashley Zhou


#always keep abundance numeric!
rownames(raw_L6)<-raw_L6[,1]
raw_L6<-raw_L6[,-1]
is.numeric(raw_L6)
t_L6 = t(raw_L6) #transpose
clean_l6=t_L6[,c(1,2,3,7,8,9,11)]
is.numeric(t_L6)  #numeric
colnames(clean_l6)<-c("ASF519","ASF457","ASF500","ASF360","ASF356","ASF492","E-Coli")

#import mapping file
mapping<-read.csv("F0F1_mapping.txt",sep="\t",header=T)
clean_mapping<-mapping[,c(1,6)]
rownames(clean_mapping)=clean_mapping[,1]

#sample = row.names(clean_l6)
#clean_l6<-cbind(sample,clean_l6)
full<-merge(clean_l6,clean_mapping, by=0)
full<-full[,c(-1,-9)]
is.numeric(full[,7])    #numeric
#parse group
full$Groups<-as.character(full$Groups)
full$gen<-rep(NA,nrow(full))
full$DSS<-rep(NA,nrow(full))
full$LF82<-rep(NA,nrow(full))
full$NEO<-rep(NA,nrow(full))
for (i in 1:nrow(full)){
  if (strsplit(full$Groups[i]," +")[[1]][1]=="F0" || strsplit(full$Groups[i]," +")[[1]][1]=="F1") full$gen[i] = strsplit(full$Groups[i]," +")[[1]]
  if (!(is.na(strsplit(full$Groups[i]," +")[[1]][2])) && strsplit(full$Groups[i]," +")[[1]][2]=="LF82+DSS") {
    full$DSS[i] = 1
    full$LF82[i]=1
  }
  else if(!(is.na(strsplit(full$Groups[i]," +")[[1]][3])) && strsplit(full$Groups[i]," +")[[1]][3]=="LF82+DSS"){
    full$DSS[i] = 1
    full$LF82[i]=1
  }
  else if(!(is.na(strsplit(full$Groups[i]," +")[[1]][2])) && strsplit(full$Groups[i]," +")[[1]][2]=="LF82"){
    full$LF82[i]=1
  }
  else if(!(is.na(strsplit(full$Groups[i]," +")[[1]][3])) && strsplit(full$Groups[i]," +")[[1]][3]=="LF82"){
    full$LF82[i]=1
  }
  if(!(is.na(strsplit(full$Groups[i]," +")[[1]][2])) && strsplit(full$Groups[i]," +")[[1]][2]=="NEO"){
    full$NEO[i]=1
  }
  else if (full$Groups[i]=="DSS"){
    full$DSS[i]=1;
  }
}
for (j in 1: nrow(full)){
  if (is.na(full$DSS[j])) full$DSS[j]=0
  if (is.na(full$LF82[j])) full$LF82[j]=0
  if (is.na(full$NEO[j])) full$NEO[j]=0
  if (is.na(full$gen[j])) full$gen[j]="control"
}
attach(full)
is.numeric(ASF519)
test1<-aov(c(ASF519,ASF457)~as.factor(DSS),data=full)
summary(test1)
detach(full)
is.numeric(as.matrix(full[,c(1,2,3,4,5,6)]))

#correlation between E-coli and LF82
cor(full$LF82, full$'E-Coli')
#0.4792


#before doing MANOVA, we assess the multivariate normality of our response variableS;
#1. chi's squared plot (a multivariate equivalent of a q-q plot)
#install.packages("MVN")
library(MVN)
uniPlot(full[,4],type="qqplot")

mardiaTest(as.matrix(full[,c(1,2,3,4,5,6)]),qqplot=T)
hzTest(as.matrix(full[,c(1,2,3,4,5,6)]),qqplot=T)
roystonTest(as.matrix(full[,c(1,2,3,4,5,6)]),qqplot=T)

mardiaTest(lmmodel1$residuals,qqplot=T)  #the residual is also not MVN

#the default "manova() function in R-stats" uses Type 1 SS
mosdel1<-manova(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)+as.factor(LF82)+as.factor(NEO),data=full)
summary(model1,test="Hotelling-Lawley")
#test options are
#"Pillai"
#"Wilks"
#"Hotelling-Lawley"
#"Roy"
#correlation matrix between the dependent variables
cor(as.matrix(full[,c(1,2,3,4,5,6)]))

model2<-manova(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82),data=full)
summary(model2)

model3<-manova(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full)
summary(model3)

model4<-manova(as.matrix(full[1:10,c(1,2,3,4,5,6)])~as.factor(DSS),data=full[1:10,])
summary(model4,test="Wilks")

model5<-lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82),data=full)
anova(model5)
#install.packages("car")
#install.packages("heplots")
library(car)
library(heplots)
etasq(model5)

#the Manova() function in the "car" package uses Type 2/ Type 3 SS
model6<-Manova(lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full))
summary(model6)

#Type 3 is not very reliable
model7<-Manova(lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full),type=3)
summary(model7)

etasq(model6)
etasq(model7)

lmmodel1<-lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full)
predict(lmmodel1)

#kruskal.test() doesn't work on multivariate case
kruskal.test(,data=full)

#install.packages('MNM')
library(MNM)
mnm1<-mv.l1lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full)
summary(mnm1)
mnm2<-mv.l1lm(as.matrix(full[,c(1,2,3,4,5,6)])~gen+as.factor(DSS)*as.factor(LF82)+as.factor(NEO),data=full,scores="r")
summary(mnm2)
