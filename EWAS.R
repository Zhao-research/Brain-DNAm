SVs=read.delim("SVs_AD.txt")
dt=read.delim("dataset_phenotype.txt")
dnam=read.delim("M_QN_BMIQ.txt")

dt=dt[dt[,'ID_o']%in%rownames(SVs),]
dim(dt)
SVs=SVs[match(dt[,'ID_o'],rownames(SVs)),]
dt=cbind(dt,SVs)
celltype=read.delim("proportion.cell type.txt",header=T)
dt=merge(dt, celltype, by.x = "ID_o",by.y='id') 

dnam=dnam[,id%in%dt[,'ID_o']]
id=id[id%in%dt[,'ID_o']]
dim(dnam)
dt=dt[match(id,dt[,'ID_o']),]

output <- NULL
for (i in (1:nrow(dnam))){
  dt[,'cpg']=dnam[i,]
  m <- lm(cpg~outcome+msex+age_death+pmi+race+V1+V2+V3+prop.NeuN.pos+factor(Sample_Plate), data = dt)
  ctable <- coef(summary(m))  
  output <- rbind(output, ctable[rownames(ctable)=='outcome',]) 
  print(i)
}
output <- as.data.frame(output)
colnames(output) <- c("beta","se","t","p")
rownames(output) <- rownames(dnam)
output$fdr <- p.adjust(output$p, method = "fdr")
write.table(output,'output_AD.txt',row.names=F,col.names=T,sep="\t",quote = FALSE)
