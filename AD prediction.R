########################################################################
load("data_pheno_AD prediction_NBB.RData")
dim(a)
cpg=read.csv('list of AD-related CpGs.csv')
cpg.sig=unique(c(cpg[,'AD.NBB'],cpg[,'Abeta.NBB'],cpg[,'Tau.NBB']))

###### regress out age, sex, cell proportion, batch variable
for (i in 1:length(cpg.sig)){
  a[,'cpg']=a[,colnames(a)==cpg.sig[i]]
  f=lm(cpg~msex+age_death+race+pmi+V1+V2+V3+prop.NeuN.pos+factor(Sample_Plate),a)
  a[,colnames(a)==cpg.sig[i]]=residuals(f)
}
###### regress out age, sex, cell proportion, batch variable

library(glmnet)
alpha.1=0 
glmnet.Training.CV=cv.glmnet(as.matrix(a[,colnames(a)%in%cpg.sig]),a[,'AD'],nfolds=10,alpha=alpha.1,family="binomial")
pred.rosmap=predict(glmnet.Training.CV,as.matrix(a[,colnames(a)%in%cpg.sig]),type="response",s=lambda.glmnet.Training)

###### prediction in NBB 
pheno=read.delim("~/NBB.txt")
dim(pheno) 
load("~/NBB/M_QN_imputed.RData")
dnam=dnam[rownames(dnam)%in%cpg.sig,]
nbb=cbind(pheno,dnam)

###### regress out age, sex, cell proportion, batch variable
for (i in 1:length(cpg.sig)){
  nbb[,'cpg']=nbb[,colnames(nbb)==cpg.sig[i]]
  f=lm(cpg~msex+age_death+pmi+V1+V2+V3+prop.NeuN.pos+factor(Sample_Plate),nbb)
  nbb[,colnames(nbb)==cpg.sig[i]]=residuals(f)
}
###### regress out age, sex, cell proportion, batch variable

dt.predict=as.matrix(nbb[,colnames(nbb)%in%cpg.sig])
pred.cpg=predict(glmnet.Training.CV,dt.predict,type="response",s=lambda.glmnet.Training)

library(pROC)
roc_1<-roc(a[,'AD'],as.vector(pred.rosmap),plot=TRUE,col="blue",lwd=1,print.auc=TRUE,print.auc.pattern = "AUROC:%.3f",legacy.axes=TRUE,main="",print.auc.y=0.8,print.auc.x=0.88)
roc_2<-roc(nbb[,'AD'],as.vector(pred.cpg),plot=TRUE,col="red",add=T,lwd=1,print.auc=TRUE,print.auc.pattern = "AUROC:%.3f",legacy.axes=TRUE,print.auc.y=0.66,print.auc.x=0.62)
legend(0.7,0.4,legend=c('RADC','NBB'),col=c("blue",'red',"darkgreen",'black'),text.col=c("blue",'red',"darkgreen",'black'),lwd=2,text.font=1,cex=1,pt.cex=1,bty='n',y.intersp=0.3,seg.len = 0.8)

