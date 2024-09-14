#####Simulate dataset
library(tidyverse)
library(bsseq)
library(dmrseq)
library(data.table)

#load dataset
bs_filtered <- readRDS(snakemake@input[['rds']])
#loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_input_obj, type="Cov")==0) == 0)
#bs_filtered <-  bs_input_obj[loci.idx, ]
#bs_filtered <- bs_input_obj[seqnames(bs_input_obj) %in% paste0("chr",c(1:22,"X","Y")),]
Covariate <- snakemake@params[['Covariate']]
number_of_DMRs <- snakemake@params[['numDMRs']]
only_normal <-snakemake@params[['only_normal']]

if(only_normal == "TRUE"){
  if(Covariate =="Patient outcome"){
    bs <- bs_filtered[,bs_filtered$Patient_outcome == "non-lethal"]
  }else{
    bs <- bs_filtered[,bs_filtered$Tissue_type == "Adjacent normal"]
  }
}else if(only_normal == "FALSE"){
 bs <- bs_filtered[,bs_filtered$Sample_title != "514C"]
}


#controlgroup = non-lethal
#non_lethal_bs <- bs_filtered[,bs_filtered$Patient_outcome == "non-lethal"]

#if tumour <-> adjacent normal
#adj_normal <- bs_filtered[,bs_filtered$Tissue_type == "Adjacent normal"]


control_bs.sim <- simDMRs(bs=bs, 
                         num.dmrs=number_of_DMRs)

####from their github
#' # show the simulated DMRs GRanges object
#' show(BS.chr21.sim$gr.dmrs)
#' 
#' # show the updated BSseq object that includes the simulated DMRs
#' show(BS.chr21.sim$bs)
#' 
#' # examine effect sizes of the DMRs
#' head(BS.chr21.sim$delta)
#####################################

simulated_bs <- control_bs.sim$bs
simulated_bs$sample_names <- rownames(pData(simulated_bs)) 
simulated_meta <- as.data.frame(pData(simulated_bs))
simulated_meta$SA <-simulated_meta$sample_names
simulated_meta<- separate(simulated_meta,col = "SA",into = c("condition", "rep"), sep = "_")
pData(simulated_bs) <- simulated_meta

#need function -> input bs_seq output .bismark.cov auch part of workflow
write_simulated_dimmer <- function(bsobj){
  meth_data <- getCoverage(bsobj, type = "M", what = "perBase")
  unmeth_data <- getCoverage(bsobj, type = "Cov", what = "perBase") - meth_data
  cov_vals <- getCoverage(bsobj, type = "Cov", what = "perBase")
  perc_data <- (meth_data/cov_vals)* 100
  gr <- granges(bsobj)
  for (name in sampleNames(bsobj)) {
    print(name)
    cov_data <- data.frame(
      chr = as.character(seqnames(gr)),
      start = start(gr),
      end = start(gr),
      perc = perc_data[,name],
      methylated = meth_data[, name],
      unmethylated = unmeth_data[, name]
    )
  #View(cov_data)
    filename <- paste0(snakemake@output[['out_folder']],name,".cov")
    write.table(cov_data, file = filename, row.names=FALSE, sep="\t", quote = FALSE, col.names = FALSE)
  }
  sampleNames(bsobj)
  sa_table <- data.table("sample_names"=sampleNames(bsobj))
  sa_table <- sa_table[, sample :=paste0(snakemake@output[['out_folder']],sample_names,".cov")]
  sa_table$SA <- sa_table$sample_names
  sa_table <- separate(sa_table,col = "SA",into = c("condition", "rep"), sep = "_")
  sa_table$condition <- gsub("Condition", "", sa_table$condition)
  sa_table$rep <- gsub("Rep", "", sa_table$rep)
  write.csv(sa_table, file = paste0(snakemake@output[['out_folder']],"sample_annotation_simulated.csv"), row.names=FALSE, quote = FALSE)
}


write_simulated_dimmer(simulated_bs)
saveRDS(simulated_bs,snakemake@output[['out_obj']])#snakemake@output[['out_obj']]
saveRDS(control_bs.sim,snakemake@output[['dmr_obj']])#snakemake@output[['dmr_obj']]
simulated_DMRs_dt <- as.data.table(control_bs.sim$gr.dmrs)
simulated_DMRs_dt$L <- control_bs.sim$dmr.L
simulated_DMRs_dt$delta <- control_bs.sim$delta
simulated_DMRs_dt$mncov <- control_bs.sim$dmr.mncov

write.table(simulated_DMRs_dt,snakemake@output[['dmr_tab']], sep = "\t", quote = FALSE)
#,'/Users/leonoraka/Desktop/MA-project/Dataset/Simulated_DMR_dt.tsv'

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(simulated_bs, type="Cov")==0) == 0)
simulated_bs_filtered <-  simulated_bs[loci.idx, ]
saveRDS(simulated_bs_filtered,snakemake@output[['filtered_out_obj']])#snakemake@output[['filtered_out_obj']]
write.table(as.data.table(simulated_bs_filtered@rowRanges),snakemake@output[['all_CpGs']], sep = "\t", quote = FALSE)



