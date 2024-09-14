#Pipeline-DMR-seq 
#DMR seq

library(dmrseq)
testCovariate <- snakemake@params[['Covariate']]
bs_filtered <- readRDS(snakemake@input[['rds']])
#loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_input_obj, type="Cov")==0) == 0)
#bs_filtered <-  bs_input_obj[loci.idx, ]
regions <- dmrseq(bs=bs_filtered,
                  cutoff = 0.05,
                  testCovariate=testCovariate)

sigRegions <- regions[regions$qval < 0.05,]


write.csv(as.data.frame(regions), 
          file=snakemake@output[['out']], quote=FALSE)

write.csv(as.data.frame(sigRegions), 
          file=snakemake@output[['sig']], quote = FALSE)





