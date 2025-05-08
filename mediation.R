load("data_pheno_AD prediction_NBB.RData")  
dim(dt)
mediation_cpg=read.csv('38 CpGs.csv')
mediation_cpg=mediation_cpg[,'cpg']

re=read.delim("result/bacon_AD.txt")
dim(re)
mediation_cpg.pos=re[as.numeric(re[,'beta'])>0,'cpg'] # 10 CpGs with beta>0
mediation_cpg.neg=re[as.numeric(re[,'beta'])<0,'cpg'] # 28 CpGs with beta<0

### DNAm score
dt[,'cpg']=0
for (i in 1:length(mediation_cpg.pos)){ 
  dt[,'cpg']=dt[,'cpg']+as.numeric(dt[,colnames(dt)==mediation_cpg.pos[i]])
}
i
for (i in 1:length(mediation_cpg.neg)){ 
  dt[,'cpg']=dt[,'cpg']-as.numeric(dt[,colnames(dt)==mediation_cpg.neg[i]])
}

################ mediation analysis ################
### regress out age, sex, cell proportion, batch variable
f=lm(cpg~msex+age_death+race+pmi+V1+V2+V3+prop.NeuN.pos+factor(Sample_Plate),dt)
dt[,'cpg']=residuals(f)
### regress out age, sex, cell proportion, batch variable
library(mediation)
# residual
f1<-lm(cpg~x,dt) 
f2<-glm(y~x+cpg+age_death+msex+race,dt,family=binomial())
set.seed(1234)  
med<- mediate(f1, f2, treat = "x", mediator = "cpg",sims=10000) # boot=T
a=summary(med)
te=c(a$tau.coef,a$tau.ci,a$tau.p)
acme=c(a$d.avg,a$d.avg.ci,a$d.avg.p)
ade=c(a$z.avg,a$z.avg.ci,a$z.avg.p)
prop_med=c(a$n.avg,a$n.avg.ci,a$n.avg.p)
res=c(te,acme,ade,prop_med)
res
colnames(res)=c('beta (total effect)','lower CI (TE)','upper CI (TE)',
                'P (TE)','beta (ACME)','lower CI (ACME)','upper CI (ACME)','P (ACME)',
                'beta (ADE)','lower CI (ADE)','upper CI (ADE)','P (ADE)',
                'beta (prop med)','lower CI (prop med)','upper CI (prop med)','P (prop med)')
write.table(res,'mediation result_DNAm score.csv',sep=',',row.names = F, col.names=T)

