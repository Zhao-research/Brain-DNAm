a=read.csv('cell.type.csv')
id=read.csv('data_phenotype.csv')
snRNA=read.table('437 snRNA-seq IDs.csv',sep='\t', header = T)
id=id[id[,'projid']%in%snRNA[,1],]
a=a[a[,'individualID']%in%id[,'individualID'],]
id=id[match(a[,'individualID'],id[,'individualID']),]
a[,'projid']=id[,'projid']

### merge clinical data ###
library(dplyr)
library(tibble)
df=read.csv('cell-annotation.csv')
table(df[,'seurat.object'])
main.batch <- df %>% mutate(batch = gsub("-[A|B]$", "", batch)) %>% 
  count(individualID, batch) %>% 
  group_by(individualID) %>% 
  slice_max(order_by=n, n=1) %>% 
  column_to_rownames("individualID") %>% 
  dplyr::select(batch)

id=read.csv('data_phenotype.csv')
### 437 participants after QC
snRNA=read.table('437 snRNA-seq IDs.csv',sep='\t', header = T)
dim(snRNA)
id=id[id[,'projid']%in%snRNA[,1],]
id=id[match(rownames(main.batch),id[,'individualID']),]
id=cbind(id,main.batch,rownames(main.batch))

require(openxlsx)
RF=read.xlsx("dataset_pheno.xlsx")
RF[,'projid']=as.numeric(RF[,'projid'])
RF=RF[RF[,'projid']%in%id[,'projid'],] 
dim(RF) # 437 participants have clinical data & snRNA data
RF=merge(RF,id,by='projid')

### for each of 7 cell types 
cell.type=unique(a[,'cell'])
cell.type
# 'excitatory' 'Oligodendrocytes' 'OPCs' 'vascular.niche' 'astrocytes' 'inhibitory' 'microglia'
re.all=NULL
for (k in 1:length(cell.type)){
  expression=a[a[,'cell']==cell.type[k],]
  expression=expression[match(RF[,'projid'],expression[,'projid']),]
  dt=cbind(expression,RF)
  
  pair=read.csv('replicated genes.csv')
  length(unique(pair[,'gene']))
  pair=pair[pair[,'gene']%in%colnames(dt),]
  gene.list=unique(pair[,'gene'])
  dim(pair)
  for(i in 1:length(gene.list)){
    dt[,colnames(dt)==gene.list[i]]=qnorm((rank(dt[,colnames(dt)==gene.list[i]],na.last='keep')-0.5)/sum(!is.na(dt[,colnames(dt)==gene.list[i]])))
  }
  beta=se=p=rep(0,length(gene.list))
  for (i in 1:length(gene.list)){
    dt[,'gene']=dt[,colnames(dt)==gene.list[i]]
    f=lm(gene~AD+age_death+msex+pmi+batch,dt) 
    tab=coef(summary(f))
    beta[i]=tab[2,'Estimate']
    se[i]=tab[2,'Std. Error']
    p[i]=tab[2,4]
  }
  fdr=p.adjust(p, method = "fdr")
  sum(fdr<0.05)
  gene.list[fdr<0.05]
  re=cbind(cell.type[k],gene.list,beta,se,p)
  re.all=rbind(re.all,re)
  print(k)
}
fdr=p.adjust(as.numeric(re.all[,'p']), method = "fdr")
sum(fdr<0.05)
unique(re.all[fdr<0.05,2]) 
unique(re.all[as.numeric(re.all[,'p'])<0.05,2])

re.all=data.frame(re.all)
re.all[,'fdr']=fdr
cell.name=c('Oligodendrocytes','OPCs','vascular.niche','astrocytes','inhibitory','microglia','excitatory')
cell.alt=c('Oli','OPC','VN','Ast','Inh','Mic','Exc')
for (i in 1:length(cell.name)){
  re.all[re.all[,'V1']==cell.name[i],'V1']=cell.alt[i]
}

######## heatmap (each row is one gene, each column is one cell type)
re.all[,'sig']=ifelse(as.numeric(re.all[,'p'])<0.05,'*','')
re.all[as.numeric(re.all[,'fdr'])<0.05,'sig']='**'
min(as.numeric(re.all[,'beta']))
max(as.numeric(re.all[,'beta']))
re.all[,'beta']=as.numeric(re.all[,'beta'])

re.all=re.all[re.all[,2]%in%re.all[as.numeric(re.all[,'fdr'])<0.05,2],]
require("ggrepel")
ggplot(re.all, aes(x = V1, y = gene.list)) + geom_raster(aes(fill=beta)) + 
  scale_fill_gradient2(low = "blue", high = "red", mid='white', midpoint = 0,limit=c(-1,1),breaks=c(-1,0,1))+   # ,limit=c(-0.2,0.2),breaks=c(-0.2,0,0.2)
  labs(y=NULL, x=NULL, fill=expression(beta))+theme_bw()+theme_minimal() + 
  theme(axis.text.x = element_text(size = 11,face="bold"),legend.title=element_text(size=11,face="bold"),panel.grid.major=element_blank(),legend.text=element_text(size=11),
        plot.title=element_text(size=11,face="bold"),panel.grid.minor=element_blank(),axis.text.y = element_text(size = 11,face="bold"),legend.position="right")+
  labs(title="",x ="", y = "")+geom_text(aes(label=sig), color="black", size=7)
