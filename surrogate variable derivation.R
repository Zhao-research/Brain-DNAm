library(sva)
pheno=read.delim("Project_Samples_Table.txt")
dim(pheno) 
pheno[,'ID_o']=paste0(pheno[,"Sentrix.Barcode"],'_',pheno[,"Sample.Section"])
pheno[,'slide']=substring(pheno[,'BCD_Well'], 1, 1)
table(pheno[,'slide'])
pheno[,'plate_row']=substring(pheno[,'BCD_Well'], 2)
table(pheno[,'plate_row'])
celltype=read.delim("proportion.cell type.txt",header=T)
pheno=merge(pheno, celltype, by.x = "ID_o",by.y='id') 
dim(pheno)

data_cpg=read.delim("M_QN_BMIQ.txt")
pheno=pheno[pheno[,'ID_o']%in%colnames(data_cpg),]
dim(pheno) 
pheno=pheno[match(colnames(data_cpg),pheno[,'ID_o']),]

mod = model.matrix(~outcome+factor(msex)+age_death+pmi+race+prop.NeuN.pos, data = pheno)
print(dim(mod)) 
colnames(mod)
mod0 = model.matrix(~factor(msex)+age_death+pmi+race+prop.NeuN.pos, data = pheno)
print(dim(mod0)) 

n.sv=100
svobj = sva(data_cpg,mod,mod0,n.sv=n.sv) 
cpg_sv <- svobj$sv
dim(cpg_sv) 
rownames(cpg_sv)=colnames(data_cpg)
write.table(cpg_sv,'SVs_AD.txt',row.names=T,col.names=T,sep="\t",quote = FALSE)
