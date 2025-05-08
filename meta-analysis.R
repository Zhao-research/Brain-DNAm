############################# RADC 
re1=read.delim("result/bacon_AD.txt") 
dim(re1) 

############################# NBB 
output=read.delim("output_braak_nft.txt")
###### Bacon correction
library(bacon)
bc <- bacon(effectsizes=output$beta,standarderrors=output$se)
output[,'p']=pval(bc)
output[,'beta']=bacon::es(bc)
output[,'se']=bacon::se(bc)
output$fdr <- p.adjust(output$p, method = "fdr")
###### Bacon correction
dim(output) 
output=output[rownames(output)%in%re1[,'cpg'],]
dim(output)
re2=output
dim(re2) 


############################# ROSMAP 450k 
output=read.delim("output_ROSMAP_450k.txt")
dim(output) 
###### Bacon correction
bc <- bacon(effectsizes=output$beta,standarderrors=output$se)
output[,'p']=pval(bc)
output[,'beta']=bacon::es(bc)
output[,'se']=bacon::se(bc)
output=output[output[,'cpg']%in%re1[,'cpg'],]
re3=output
dim(re3) 


############################# MS 450k
output=read.delim("output_PFC_Mount.Sinai.txt")
dim(output) 
###### Bacon correction
bc <- bacon(effectsizes=output$beta,standarderrors=output$se)
output[,'p']=pval(bc)
output[,'beta']=bacon::es(bc)
output[,'se']=bacon::se(bc)
output[,'cpg']=rownames(output)
output=output[output[,'cpg']%in%re1[,'cpg'],]
re4=output
dim(re4) 


############################# London.1 450k
output=read.delim("output_PFC_London.1.txt")
dim(output) 
###### Bacon correction
bc <- bacon(effectsizes=output$beta,standarderrors=output$se)
output[,'p']=pval(bc)
output[,'beta']=bacon::es(bc)
output[,'se']=bacon::se(bc)
output[,'cpg']=rownames(output)
output=output[output[,'cpg']%in%re1[,'cpg'],]
re5=output
dim(re5) 


#################
array450=c(re3[,'cpg'],re4[,'cpg'],re5[,'cpg'])
dim(re2) 
re2=re2[rownames(re2)%in%array450,]
dim(re2) 
re2[,'cpg']=rownames(re2)
all=c(re1[,'cpg'],re2[,'cpg'],re3[,'cpg'],re4[,'cpg'],re5[,'cpg'])
summary_cpg=data.frame(table(all))
n_rep=3
cpg.list=summary_cpg[,1]
out=cbind(as.matrix(cpg.list),cohort)
n0=nrow(out)
### meta-analysis
library(meta)
n0=length(cpg.list)
res.meta=matrix(data = NA, nrow = n0, ncol = 10)
for (i in 1:n0){
  yi=c(re1[re1[,'cpg']==cpg.list[i],'beta'],re2[re2[,'cpg']==cpg.list[i],'beta'],re3[re3[,'cpg']==cpg.list[i],'beta'],re4[re4[,'cpg']==cpg.list[i],'beta'],re5[re5[,'cpg']==cpg.list[i],'beta']) 
  vi=c(re1[re1[,'cpg']==cpg.list[i],'se'],re2[re2[,'cpg']==cpg.list[i],'se'],re3[re3[,'cpg']==cpg.list[i],'se'],re4[re4[,'cpg']==cpg.list[i],'se'],re5[re5[,'cpg']==cpg.list[i],'se']) 
  out=metagen(yi,vi) 
  res.meta[i,]<-c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, out$Q, out$pval.Q, 1-pchisq(out$Q, out$df.Q))
  print(i)
}
new=data.frame(cbind(as.matrix(cpg.list),beta,se,CI_l,CI_u,p_sig,Q,Q_p,I2))
colnames(new)=colnames(new_all)
new=rbind(new_all,new)
dim(new) 
new[,'fdr']=p.adjust(as.numeric(new$p_sig), method = "fdr")
sum(as.numeric(new[,'p_sig'])<0.05) 
sum(as.numeric(new[,'fdr'])<0.05) 
write.csv(new,file='meta.csv', row.names=FALSE)

