re=read.delim("result/R1/bacon_AD.txt") # data contains 38 DMPs and regression coefficient

########################################### replication in BSHR
a=read.delim("/external_cohorts/BSHRI/BSHRI.txt",header=T)
a[,'sample_bcd_barcode.ch1']=gsub("\\-BCD.*","",a[,'sample_bcd_barcode.ch1'])
table(a[,'brain.region.ch1'],useNA='always') # IFG/STG
a[,'nia.r.ch1']=ifelse(a[,'nia.r.ch1']%in%c('high','intermediate'),1,0)
table(a[,'nia.r.ch1']) 

# IFG
load("/external_cohorts/BSHRI/M.value_IFG.RData")
cpg[1:3,1:4]
colnames(cpg)=sub('X','',colnames(cpg))
cpg.1=cpg
# STG
load("/external_cohorts/BSHRI/M.value_STG.RData")
cpg[1:3,1:4]
colnames(cpg)=sub('X','',colnames(cpg))
cpg=cbind(cpg.1,cpg)
a=cbind(a,cpg)

### weighted DNAm score
a[,'cpg']=0
for (i in 1:nrow(re)){ 
  a[,'cpg']=a[,'cpg']+as.numeric(a[,colnames(a)==re[i,'cpg']])*re[i,'beta']
}

var.con=c('braak.score.ch1','age.at.death.ch1','pmi.ch1','neun_pos.ch1','sva1.ch1','sva2.ch1','sva3.ch1')
for (i in 1:length(var.con)){
  a[,var.con[i]]=as.numeric(a[,var.con[i]])
}
a[,'ad']=a[,'nia.r.ch1'] 
library(lme4)
library(lmerTest)
m <- lmer(cpg~ad+age.at.death.ch1+factor(gender.ch1)+pmi.ch1+neun_pos.ch1+factor(brain.region.ch1)+factor(sample_bcd_barcode.ch1)+sva1.ch1+sva2.ch1+sva3.ch1+(1|donorid.ch1), data = a)
summary(m) 

############ plot ############
library(ggforce)
ggplot(a, aes(x=factor(brain.region.ch1),y=cpg,fill=factor(ad)))+#scale_x_discrete(breaks = c("No", "Yes"))+
  labs(y = "DNAm score",x = "") + geom_boxplot(outlier.shape = NA,width=0.5,notch=T,position = position_dodge(width = 0.7))+
  theme_bw()+scale_x_discrete(labels=c('inferior frontal\n gyrus (n=117)','superior temporal\n gyrus (n=127)'))+
  theme(plot.title = element_text(size=8,hjust = 0.5),axis.title = element_text(size=7),
        axis.text=element_text(size=6),legend.position = "none",panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_manual(name = '',values = c("limegreen",'lightcoral'))#+scale_y_continuous(limits=c(-1,-0.15))



############################ replication in ROSMAP EPIC v1 
a=read.delim("multi-region_data/data ROSMAP_PFC.txt")
a[,'outcome']=a[,'ad'] 
dim(a)

### weighted DNAm score
a[,'cpg']=0
for (i in 1:nrow(re)){ 
  a[,'cpg']=a[,'cpg']+as.numeric(a[,colnames(a)==re[i,'cpg']])*re[i,'beta']
}
a.pfc=a

a=read.delim("multi-region_data/data ROSMAP_ST.txt")
a[,'outcome']=a[,'ad'] 
dim(a)
### weighted DNAm score
a[,'cpg']=0
for (i in 1:nrow(re)){ 
  a[,'cpg']=a[,'cpg']+as.numeric(a[,colnames(a)==re[i,'cpg']])*re[i,'beta']
}
a.st=a
a=rbind(a.pfc,a.st)

m <- lmer(cpg~outcome+msex+age_death+race+pmi+V1+V2+V3+prop.NeuN.pos+region+(1|projid), data = a)
coef(summary(m)) 

############ plot ############
library(ggforce)
ggplot(a, aes(x=factor(region),y=cpg,fill=factor(ad)))+
  labs(y = "DNAm score",x = "") + geom_boxplot(outlier.shape = NA,width=0.5,notch=T,position = position_dodge(width = 0.7))+
  theme_bw()+#scale_x_discrete(labels=c('PFC (n=50)','striatum (n=50)'))+
  theme(plot.title = element_text(size=8,hjust = 0.5),axis.title = element_text(size=7),
        axis.text=element_text(size=6),legend.position = "none",panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_fill_manual(name = '',values = c("limegreen",'lightcoral'))

