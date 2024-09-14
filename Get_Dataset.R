#Filter und BSmooth Dataset
library(bsseq)
library(XLConnect)
library(parallel)
library(BiocParallel)
library(tidyverse)
library(data.table)

##Prostate cancer
##BSseq object and Dimmer .cov files
print("Step 1: Set up object")
read_in_dataset <- read.table("/Users/leonoraka/Downloads/BigTable.tsv")
read_in_dataset<- read_in_dataset %>% rename_all(~stringr::str_replace_all(.,"^X",""))
read_in_dataset<- read_in_dataset %>% rename_all(~stringr::str_replace(.,"[.]C","_M"))
read_in_dataset<- read_in_dataset %>% rename_all(~stringr::str_replace_all(.,"[.]cov","_Cov"))

All_M <- select(read_in_dataset,ends_with("_M"))
All_M <- All_M %>% rename_all(~stringr::str_replace(.,"_M",""))

All_Cov <- select(read_in_dataset,ends_with("_Cov"))
All_Cov <- All_Cov %>% rename_all(~stringr::str_replace(.,"_Cov",""))

bsseq_whole <- BSseq(chr = read_in_dataset[,"chr"], pos = read_in_dataset[,"position"],
                     M = as.matrix(All_M), Cov = as.matrix(All_Cov), sampleNames = colnames(All_M))

print("Step 1 done")
#Metadata
print("Step 2: Add metadata & save whole object")
read_metadata <- readWorksheet(loadWorkbook("/Users/leonoraka/Desktop/GSE158927_series_matrix.xlsx"),sheet="Metadata")
rownames(read_metadata)<- read_metadata[,"Sample_title"]
pData(bsseq_whole)<- read_metadata

saveRDS(bsseq_whole, snakemake@output[["whole"]])

print("Step 2 done")
#Filter object
print("Step 3: filter object and save")
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bsseq_whole, type="Cov")==0) == 0)
bs_filtered <-  bsseq_whole[loci.idx, ]
bs_filtered <- bs_filtered[seqnames(bs_filtered) %in% paste0("chr",c(1:22,"X","Y")),]

saveRDS(bs_filtered, snakemake@output[["filtered"]])
print("Step 3 done")

#smooth
print("Step 4: smooth object and save")
bs_filtered_smoothed <- BSmooth(bs_filtered,BPPARAM = MulticoreParam(workers = 8, progressbar = TRUE), verbose = TRUE)
saveRDS(bs_filtered_smoothed, snakemake@output[["smoothed"]])
print("Step 4 done")

##Dimmer 
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
    filename <- paste0(snakemake@output[["dimmer_folder"]],name,".cov")
    write.table(cov_data, file = filename, row.names=FALSE, sep="\t", quote = FALSE, col.names = FALSE)
  }
  #sampleNames(bsobj)
}
print("Step 5: write dimmer files")
write_simulated_dimmer(bs_filtered)
print("Step 5 done")
#sample annotation files already set up

##Xie Data 
##BSseq object and Dimmer .cov files






















