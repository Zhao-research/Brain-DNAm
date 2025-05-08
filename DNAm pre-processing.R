library(minfi)
library(wateRmelon)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

baseDir='~/data'
targets=read.metharray.sheet(baseDir) 
dim(targets)
rgSet=read.metharray.exp(targets=targets,extended = TRUE) 

annotation(rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(rgSet)["annotation"] = "20a1.hg38"
annotation=getAnnotation(rgSet)

aa=rownames(rgSet)

###-------------------------------------------------------------------------------------------------------------------
### QC: Detecting samples with bisulfite conversion efficiency < 80%
###-------------------------------------------------------------------------------------------------------------------
#calculate bisulfite conversion statistic using bscon (wateRmelon package) and plots histogram
bs<-bscon(rgSet)
quantile(bs)
###Detecting bad samples 
GOOD.samples.Bisulifite<-as.character(colnames(rgSet[,which(bs>=80)]))										  
length(GOOD.samples.Bisulifite) 
# select a subset of samples
rgSet <- rgSet[,GOOD.samples.Bisulifite]
dim(rgSet) 
targets[,'ID_match']=paste0(targets[,"Slide"],'_',targets[,"Array"])
targets <- targets[targets[,'ID_match']%in%GOOD.samples.Bisulifite,]
dim(targets) 

###### 
mset.raw <- mapToGenome(rgSet)
methylSex <- getSex(mset.raw)
table(methylSex$predictedSex)
head(methylSex)

pheno=read.delim("Project_Samples_Table.txt")
dim(pheno) 
pheno[,'ID_o']=paste0(pheno[,"Sentrix.Barcode"],'_',pheno[,"Sample.Section"])
require(openxlsx)
a=read.xlsx('dataset_phenotype.xlsx', sheet = 1)
a[,'projid']=as.numeric(a[,'projid'])

### 28 control samples
sample_control=pheno[!pheno[,'Sample.ID']%in%a[,'projid'],'ID_o']
a=a[a[,'projid']%in%pheno[,'Sample.ID'],]
pheno=merge(pheno,a[,c('projid','msex')],by.x='Sample.ID',by.y='projid',all.x=T)
dim(pheno) 

# select a subset of samples who have phenotype data
pheno=pheno[pheno[,'ID_o']%in%c(rownames(methylSex),sample_control),]
sample_sub=rownames(methylSex)[rownames(methylSex)%in%pheno[,'ID_o']]
sample_sub=unique(sample_sub,sample_control)
rgSet <- rgSet[,sample_sub]
head(targets)
targets[,'ID_match']=paste0(targets[,"Slide"],'_',targets[,"Array"])
targets <- targets[targets[,'ID_match']%in%sample_sub,]
methylSex=methylSex[rownames(methylSex)%in%sample_sub,]
dim(targets) 
dim(rgSet)
dim(methylSex) 

# Get IDs of sex-check failing samples
pheno=pheno[pheno[,'ID_o']%in%sample_sub,]
dim(pheno) # 90*26
pheno=pheno[order(pheno[,'ID_o']),]
table(pheno[,'msex']) # 0/1
methylSex=methylSex[order(rownames(methylSex)),]
table(methylSex[,'predictedSex']) 
methylSex[,'predictedSex']=ifelse(methylSex[,'predictedSex']=='M',1,0)
# samples not pass sex prediction
sum(methylSex[,'predictedSex']!=pheno[,'msex']) 
GOOD.samples.sex=rownames(methylSex)[methylSex[,'predictedSex']==pheno[,'msex']]
# select a subset of samples
GOOD.samples.sex=unique(GOOD.samples.sex,sample_control)
rgSet <- rgSet[,GOOD.samples.sex]
dim(rgSet) 
targets <- targets[targets[,'ID_match']%in%GOOD.samples.sex,]
dim(targets)
mset.raw=mset.raw[,GOOD.samples.sex]
dim(mset.raw) 


###pfilter function within the wateRmelon package, with the following criteria used for exclusion: samples with a detection  P > 0.05 in more than 1% of probes, 
###probes with > three beadcount in 5% of samples and  probes having 1% of samples with a detection P value > 0.05).
###-------------------------------------------------------------------------------------------------------------------
pfltSet<-pfilter(rgSet, perc=1, perCount=5, pthresh=1) 

# Create MethylSet:
mset.pflt <- preprocessRaw(pfltSet)

# Identify outliers based on the wateRmelon outlyx function. No outliers identified:
outliers <- outlyx(mset.pflt, plot = TRUE)
#outliers
BAD.samples.outliers=rownames(outliers)[outliers[,'outliers']==TRUE] # outlier sample IDs
BAD.samples.outliers
if (length(BAD.samples.outliers)>0){
  mset.pflt=mset.pflt[,!colnames(mset.pflt)%in%BAD.samples.outliers]
  dim(mset.pflt) 
  targets <- targets[!targets[,'ID_match']%in%BAD.samples.outliers,]
  dim(targets)
}

###-------------------------------------------------------------------------------------------------------------------
### QC.c : detecting the samples with the mean intensity of methylated or unmethylated signals were three standard deviations above or below the mean
###-------------------------------------------------------------------------------------------------------------------
meth=getMeth(mset.pflt)
unmeth=getUnmeth(mset.pflt)

M.mean<-apply(meth, 2, mean)
length(M.mean) 
U.mean<-apply(unmeth, 2, mean)
length(U.mean) 

Mean.M<-mean(M.mean)
Mean.U<-mean(U.mean)
sd.M<-sd(M.mean)
sd.U<-sd(U.mean)
###Detecting the bad samples 
BAD.samples.Intensity<-as.character(targets[which(Mean.M +3*sd.M < M.mean |
                        M.mean < Mean.M-3*sd.M | Mean.U +3*sd.U < U.mean |
                        U.mean< Mean.U-3*sd.U),'Sample_ID'])
BAD.samples.Intensity 

# select a subset of samples
if (length(BAD.samples.Intensity)>0){
  mset.pflt <- mset.pflt[,!BAD.samples.Intensity]
  dim(mset.pflt) 
  targets <- targets[!targets[,'ID_match']%in%BAD.samples.Intensity,]
  dim(targets) 
}

### delete control samples
sum(sample_control%in%colnames(mset.pflt)) 
mset.pflt <- mset.pflt[,!colnames(mset.pflt)%in%sample_control]
dim(mset.pflt)
### delete control samples

# Normalize using the wateRmelon dasen method
print("normalization")
dim(mset.pflt)
sum(is.na(mset.pflt))
#sum(is.nan(mset.pflt))
#sum(is.infinite(mset.pflt))
mset.pflt <- dasen(mset.pflt)
dim(mset.pflt)
print("normalization")


# Map to genome (minfi):
mset.pflt <- mapToGenome(mset.pflt) 
dim(mset.pflt) 

# Identity low quality samples based on minfi QC:
minfiQC <- getQC(mset.pflt)
plotQC(minfiQC, badSampleCutoff = 10.5)


###-------------------------------------------------------------------------------------------------------------------
# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mset.pflt,snps=c('SBE','CpG'),maf=0.05)
dim(mSetSqFlt) 

# Filter out probes on sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% annotation$Name[annotation$chr %in% c("chrX", "chrY")])
table(keep) 
mSetSqFlt <- mSetSqFlt[keep,]

## Removing probes that have been demonstrated to map to multiple places in the genome.
###-------------------------------------------------------------------------------------------------------------------
###Removing the probes with evidence for cross-hybridizition
###-------------------------------------------------------------------------------------------------------------------
xReactiveProbes <- read.csv(file="48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
cpg_name=sub("_.*", "", featureNames(mSetSqFlt))
keep <- !(cpg_name %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt
###-------------------------------------------------------------------------------------------------------------------

dasenBetas <- getBeta(mSetSqFlt) 
dim(dasenBetas)
cpg=rownames(dasenBetas)
write.table(dasenBetas,'beta_QN.txt',row.names=T,col.names=T,sep="\t",quote = FALSE)


### BMIQ (inter-arrary/intra-sample normalization)
library(wateRmelon)
library(RPMM)
beta=read.delim("beta_QN.txt")
dim(beta) 
beta[1:3,1:4]
cpg_info=read.table('CpG_info.csv',sep=',',header=T)
cpg_info=cpg_info[cpg_info[,'IlmnID']%in%rownames(beta),]
dim(cpg_info)
table(cpg_info[,'Infinium_Design_Type'])
cpg_info[,'Infinium_Design_Type']=ifelse(cpg_info[,'Infinium_Design_Type']=='I',1,2)
cpg_info=cpg_info[match(rownames(beta),cpg_info[,'IlmnID']),]

type12=cpg_info[,'Infinium_Design_Type']
set.seed (1234)
betaQN_BMIQ <- apply(beta, 2,
                     function(x){
                       re <- BMIQ(x, design.v = type12, plots = FALSE)
                       return (re$nbeta)
                     }
)
dim(betaQN_BMIQ)
sum(is.na(betaQN_BMIQ))
colnames(betaQN_BMIQ)=colnames(beta)
rownames(betaQN_BMIQ)=rownames(beta)
write.table(betaQN_BMIQ,'beta_QN_BMIQ.txt',row.names=T,col.names=T,sep="\t",quote = FALSE)
betaQN_BMIQ=log2(betaQN_BMIQ/(1-betaQN_BMIQ))
write.table(betaQN_BMIQ,'M_QN_BMIQ.txt',row.names=T,col.names=T,sep="\t",quote = FALSE)

