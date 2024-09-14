library(tidyverse)
library(Repitools)
library(bsseq)
library(dmrseq)
library(data.table)
library(DMRcate)
library(parallel)
library(BiocParallel)
#library(DSS)
bs_filtered <- readRDS(snakemake@input[['rds']])
#loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_input_obj, type="Cov")==0) == 0)
#bs_filtered <-  bs_input_obj[loci.idx, ]
Covariate <- snakemake@params[['Covariate']]

if(Covariate == "Patient outcome"){ #needs change
  bs_filtered <- bs_filtered[,bs_filtered$Patient_outcome != "NA"]
  pData(bs_filtered)$Patient_outcome <- gsub("-", "_", pData(bs_filtered)$Patient_outcome)
  tissue <- factor(pData(bs_filtered)$Patient_outcome)
  design <- model.matrix(~tissue)
  colnames(design) <- gsub("tissue", "", colnames(design))
  colnames(design)[1] <- "lethal"
  rownames(design) <- colnames(bs_filtered)
  design
  #Methylation matrix design
  methdesign <- edgeR::modelMatrixMeth(design)
  methdesign
  cont.mat <- limma::makeContrasts(lethal_vs_nlethal = lethal - non_lethal,levels=methdesign)
  seq_annot <- sequencing.annotate(bs_filtered, methdesign, all.cov = TRUE,
                                   contrasts = TRUE, cont.matrix = cont.mat,
                                   coef = "lethal_vs_nlethal", fdr=0.05)  
}else if(Covariate =="Tissue type"){ #needs change
  pData(bs_filtered)$Tissue_type <- gsub(" ", "_", pData(bs_filtered)$Tissue_type)
  tissue <- factor(pData(bs_filtered)$Tissue_type)
  design <- model.matrix(~tissue)
  colnames(design) <- gsub("tissue", "", colnames(design))
  colnames(design)[1] <- "Tumour"
  rownames(design) <- colnames(bs_filtered)
  design
  #Methylation matrix design
  methdesign <- edgeR::modelMatrixMeth(design)
  methdesign
  cont.mat <- limma::makeContrasts(tumour_vs_normal = Tumour - Adjacent_normal,levels=methdesign)
  seq_annot <- sequencing.annotate(bs_filtered, methdesign, all.cov = TRUE,
                                   contrasts = TRUE, cont.matrix = cont.mat,
                                   coef = "tumour_vs_normal", fdr=0.05)  
}else if(Covariate == "condition"){
  #pData(bs_filtered)$condition <- gsub(" ", "_", pData(bs_filtered)$Tissue_type)
  
  #samples <- rownames(pData(bs_filtered))
  #sampnames <- sub("\\..*", "", colnames(meth))[-c(1:2)]
 # DSSres <- DMLtest(bs_filtered, group1=samples[grepl("^Condition1_", samples)], group2=samples[grepl("^Condition2_", samples)], smoothing=FALSE)
#  wgbsannot <- cpg.annotate("sequencing", DSSres)
#  wgbs.DMRs <- dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1)
  condition <- factor(pData(bs_filtered)$condition)
  design <- model.matrix(~0 +condition)
  #design <- as.data.frame(design)
  colnames(design) <- gsub("condition", "", colnames(design))
  #colnames(design)[1] <- "Intercept"
  #rownames(design) <- colnames(bs_filtered)
  #design$Condition1 <- ifelse(grepl("^Condition1_", rownames(design)), 1, 0)
  #design
  #design <- as.matrix(design)

  #Methylation matrix design
  methdesign <- edgeR::modelMatrixMeth(design)
  methdesign
  cont.mat <- limma::makeContrasts(Condition1_vs_Condition2 = Condition1 - Condition2,levels=methdesign)
  seq_annot <- sequencing.annotate(bs_filtered, methdesign, all.cov = TRUE,
                                   contrasts = TRUE, cont.matrix = cont.mat,
                                   coef = "Condition1_vs_Condition2", fdr=0.05)  
}


#seq_annot
dmrcate.res <- dmrcate(seq_annot, C=2, min.cpgs = 5, pcutoff = 0.05)
#dmrcate.res
resulting.ranges <- extractRanges(dmrcate.res, genome="hg19")
result_df <- annoGR2DF(resulting.ranges)

#result_df_filtered <- result_df[abs(result_df$meandiff) >= 0.1,]
#print(result_df)
write.table(result_df, snakemake@output[['out']], quote = FALSE, row.names = FALSE, sep = ";")








