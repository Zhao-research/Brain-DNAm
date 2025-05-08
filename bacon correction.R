library(bacon)
output=read.delim(paste0("output_AD.txt"))
dim(output) 
sum(output$p<0.05)
sum(output$fdr<0.05)

bc <- bacon(effectsizes=output$beta,standarderrors=output$se)
bc
estimates(bc)
inflation(bc)
bias(bc)
beta_bacon=bacon::es(bc)
se_bacon=bacon::se(bc)
p_bacon=pval(bc)
fdr_bacon=p.adjust(p_bacon, method = "fdr")
out=cbind(beta_bacon,se_bacon,p_bacon,fdr_bacon)
colnames(out)=c('beta_bacon','se_bacon','p_bacon','fdr_bacon')

###### QQ plot ######
library(qqman)
lambda=inflation(bc)
lambda=round(lambda,2)
qqman::qq(output$p, main = paste0("lambda=",round(lambda,4)),cex.main=1, pch = 18, col = "blue4", cex = 0.8, las = 1) # xlim = c(0, 6),ylim = c(0, 12),
# bacon adjusted
bc <- bacon(effectsizes=beta_bacon,standarderrors=se_bacon)
lambda=inflation(bc)
lambda=round(lambda,2)
qqman::qq(p_bacon, main = paste0("lambda=",round(lambda,4)),cex.main=1, pch = 18, col = "blue4", cex = 0.8, las = 1) # xlim = c(0, 6),ylim = c(0, 12),
