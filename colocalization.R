
##### perform colocalization between CpG and gene

library(coloc)

cpg.gene <- read.table("rosmap_cpg_gene_pair_for_coloc.txt",header = TRUE)
folder <- c("As","En","Ex","In","Mi","Ol","Op","Pe")

best.causal.snp <- c()

for(i in 1:191){
  
  cpg <- cpg.gene[i,5]
  gene <- cpg.gene[i,6]
  
  if(is.na(gene)){
    
    next
    
  }else{
    
    if(!file.exists(paste0("./",cpg,"_mQTL_Summary_Statistics_coloc.txt"))){
      
      next
    }else if(file.size(paste0("./",cpg,"_mQTL_Summary_Statistics_coloc.txt"))==0){
      
      next
      
    }else{
      
      qtl <- read.table(paste0("./",cpg,"_mQTL_Summary_Statistics_coloc.txt"),header = TRUE)
      qtl.list <- list()
      qtl <- qtl[order(qtl[,2]),]
      qtl.list[[1]] <- qtl$beta
      qtl.list[[2]] <- (qtl$se)^2
      qtl.list[[3]] <- qtl$SNP
      qtl.list[[4]] <- qtl$pos
      qtl.list[[5]] <- "quant"
      qtl.list[[6]] <- 773
      qtl.list[[7]] <- qtl$MAF
      names(qtl.list) <- c("beta","varbeta","snp","position","type","N","MAF")
    }
    
    for(j in 1:8){
      
      if(!file.exists(paste0("./",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt"))){
        
        next
        
      }else if(file.size(paste0("./",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt")) == 58){
        
        next
        
      }else{
        
        file <- read.table(paste0("./",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt"),header = TRUE)
        file.delete1 <- which(duplicated(file[,3]) == TRUE)
        file.delete2 <- which((file[,6]<= 0)|(file[,6] >= 1))
        file.delete <- sort(unique(c(file.delete1,file.delete2)))
        
        if(length(file.delete) > 0){
          
          file <- file[-file.delete,]}
        file.list <- list()
        file <- file[order(file[,2]),]
        file.list[[1]] <- file$Beta
        file.list[[2]] <- (file$Se)^2
        file.list[[3]] <- file$SNP
        file.list[[4]] <- file$Pos
        file.list[[5]] <- "quant"
        file.list[[6]] <- 192
        file.list[[7]] <- file$MAF
        names(file.list) <- c("beta","varbeta","snp","position","type","N","MAF")
        
        possibleError <- tryCatch({res <- coloc.abf(dataset1=qtl.list,dataset2=file.list)}, error=function(e) e)
        
        
        if(inherits(possibleError, "error")){ 
          
          next
          
        }else if(res$summary[6] < 0.8){
          
          next 
          
        }else{
          
          cat(paste0(cpg,"_",gene,"_",folder[j]),"\n")
          best <- res$results[which((res$results$SNP.PP.H4 == max(res$results$SNP.PP.H4))==TRUE),c(1,2,12)]
          a1 <- cbind(cpg, gene,folder[j])
          a2 <- cbind(res$summary[1],res$summary[2],res$summary[3],res$summary[4],res$summary[5],res$summary[6])
          best.causal.snp <- rbind(best.causal.snp, cbind(a1, best, a2))
          
          write.table(subset(res$results,SNP.PP.H4>0.01), file = paste0("./eQTL_EWAS_Coloc_Sig_Res_",cpg,"_",gene,"_",folder[j],"_SNP_PP_H4_Greater_Than_1e-2.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
          o <- order(res$results$SNP.PP.H4,decreasing=TRUE)
          cs <- cumsum(res$results$SNP.PP.H4[o])
          w <- which(cs > 0.95)[1]
          write.table(res$results[o,][1:w,], file = paste0("./eQTL_EWAS_Coloc_Sig_Res_",cpg,"_",gene,"_",folder[j],"_Credible_Set.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
          
        } 
      }
    }
  }
}


colnames(best.causal.snp) <- c("CpG","Gene","CT","SNP","Pos","SNP.PP.H4", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

write.table(best.causal.snp, file = "./best_causal_snp_mQTL_eQTL_EWAS_Coloc.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


##### summarize top coloc results

eqtl.res <- read.table("./best_causal_snp_mQTL_eQTL_EWAS_Coloc.txt",header = TRUE)[,c(1,2,3,4,12)]
eqtl.res.summary <- c()

for(i in 1:dim(eqtl.res)[1]){
  
  cpg <- eqtl.res[i,1]
  gene <- eqtl.res[i,2]
  ct <- eqtl.res[i,3]
  sentinel.snp <- eqtl.res[i,4]
  
  file <- read.table(paste0("./",ct,"/",ct,"_",gene,"_eQTL_SumStats_GRCh38.txt"),header = TRUE)
  file.row.indx <- which(file[,3] == sentinel.snp)
  sentinel.snp.pos <- gsub("chr","",paste0(file[file.row.indx,1],":",file[file.row.indx,2]))
  file.minp.indx <- which(file[,7]==min(file[,7]))[1]
  file.minp <- file[file.minp.indx,c(1,2,3,7)]
  file.top.snp.pos <- gsub("chr","",paste0(file.minp[1],":",file.minp[2]))
  
  mc.ss <- read.table(paste0("./",cpg,"_mQTL_Summary_Statistics_coloc.txt"),header = TRUE)
  mc.ss.minp.indx <- which(mc.ss[,5]==min(mc.ss[,5]))[1]
  mc.ss.minp <- mc.ss[mc.ss.minp.indx,c(1,2,6,5)]
  mc.top.snp.pos <- paste0(mc.ss.minp[1],":",mc.ss.minp[2])
  
  eqtl.res.summary <- rbind(eqtl.res.summary,c(eqtl.res[i,1:4],sentinel.snp.pos,file.minp[3],file.top.snp.pos,file.minp[4],mc.ss.minp[3],mc.top.snp.pos,mc.ss.minp[4],eqtl.res[i,5]))
}

colnames(eqtl.res.summary) <- c("CpG","Gene","Cell Type","Sentinel Causal SNP","Pos","Top SNP (eQTL)","Pos (eQTL)", "Pvalue (eQTL)", "Top SNP (mQTL)", "Pos (mQTL)", "Pvalue (mQTL)", "PP4")

write.csv(eqtl.res.summary, file = "./CT_eQTL_Coloc_Top_Signal_Summary.csv", row.names = FALSE)






