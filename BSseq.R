
## ----libr, message=FALSE, warning=FALSE---------------------------------------
library(tidyverse)
library(bsseq)
library(dmrseq)
library(data.table)
library(parallel)
library(BiocParallel)

bs_filtered_smoothed <- readRDS(snakemake@input[['rds']])
#loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_input_obj, type="Cov")==0) == 0)
#bs_filtered <-  bs_input_obj[loci.idx, ]

#bs_filtered <- bs_filtered[seqnames(bs_filtered) %in% paste0("chr",c(1:22,"X","Y")),]
bs_filtered_smoothed <- BSmooth(bs_filtered_smoothed,BPPARAM = MulticoreParam(workers = 1, progressbar = TRUE), verbose = TRUE)


Covariate <- snakemake@params[['Covariate']]
BS.cov <- getCoverage(bs_filtered_smoothed)
##if Patient outcome
if(Covariate ==  "Patient outcome"){
  keepLoci.ex <- which(rowSums(BS.cov[, bs_filtered_smoothed$Patient_outcome == "lethal"] >= 2) >= 2 &
                       rowSums(BS.cov[, bs_filtered_smoothed$Patient_outcome == "non-lethal"] >= 2) >= 2)
  length(keepLoci.ex)
  bs_filtered_smoothed <- bs_filtered_smoothed[keepLoci.ex,]
  meta_info <- as.data.table(pData(bs_filtered_smoothed))
  gr1 <- meta_info[Patient_outcome=="lethal",Sample_title]
  gr2 <- meta_info[Patient_outcome=="non-lethal",Sample_title]
}else if(Covariate ==  "Tissue type"){
  keepLoci.ex <- which(rowSums(BS.cov[, bs_filtered_smoothed$Tissue_type == "Tumour"] >= 2) >= 2 &
                         rowSums(BS.cov[, bs_filtered_smoothed$Tissue_type == "Adjacent normal"] >= 2) >= 2)
  length(keepLoci.ex)
  bs_filtered_smoothed <- bs_filtered_smoothed[keepLoci.ex,]
  meta_info <- as.data.table(pData(bs_filtered_smoothed))
  gr1 <- meta_info[Tissue_type=="Tumour",Sample_title]
  gr2 <- meta_info[Tissue_type=="Adjacent normal",Sample_title]
}else if(Covariate == "condition"){
  keepLoci.ex <- which(rowSums(BS.cov[, bs_filtered_smoothed$condition == "Condition1"] >= 2) >= 2 &
                         rowSums(BS.cov[, bs_filtered_smoothed$condition == "Condition2"] >= 2) >= 2)
  length(keepLoci.ex)
  bs_filtered_smoothed <- bs_filtered_smoothed[keepLoci.ex,]
  meta_info <- as.data.table(pData(bs_filtered_smoothed))
  gr1 <- meta_info[condition=="Condition1",sample_names]
  gr2 <- meta_info[condition=="Condition2",sample_names]
}


bs_filtered_smoothed.tstat <- BSmooth.tstat(bs_filtered_smoothed, 
                                    group1 = gr1,
                                    group2 = gr2, 
                                    estimate.var = "same",#"group2", #less variability => control
                                    local.correct = TRUE,
                                    verbose = TRUE)


plot(bs_filtered_smoothed.tstat)
dmrs0 <- dmrFinder(bs_filtered_smoothed.tstat, cutoff = c(-5, 5))

dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(dmrs)
write.table(dmrs0, snakemake@output[['out']], quote = FALSE, row.names = FALSE)
write.table(dmrs, snakemake@output[['sig']], quote = FALSE, row.names = FALSE)


