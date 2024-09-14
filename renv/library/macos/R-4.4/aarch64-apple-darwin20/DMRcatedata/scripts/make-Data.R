#Make 450K filtering data
#crosshyb
require(readxl)
require(plyr)
require(GenomicFeatures)
require(Gviz)
require(rtracklayer)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

tmp <- tempfile(fileext = ".xlsx")
download.file(url = "http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx", destfile = tmp)
ch_450k <- c(read_xlsx(tmp, sheet = 1)$TargetID, read_xlsx(tmp, sheet=2)$TargetID, read_xlsx(tmp, sheet=3)$TargetID)
tab2 <- read.csv("https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv", stringsAsFactors = FALSE)
tab3 <- read.csv("https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM3_ESM.csv", stringsAsFactors = FALSE)
ch_EPIC <- unique(c(tab2$PROBE, tab3$PROBE))
crosshyb <- sort(unique(c(ch_450k, ch_EPIC)))
file.remove(tmp)

#snpsall
#EPIC
tab4 <- read.csv("https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv", stringsAsFactors = FALSE)
tab5 <- read.csv("https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv", stringsAsFactors = FALSE)
tab6 <- read.csv("https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM6_ESM.csv", stringsAsFactors = FALSE)

allprobes <- unique(c(tab4$PROBE, tab5$PROBE, tab6$PROBE))

getallprobeinf <- function(probe){
  df0 <- data.frame(snps=character(0), distances=integer(0), mafs=numeric(0))
  df <- df0
  if(probe %in% tab4$PROBE){
    t4entries <- tab4[tab4$PROBE==probe,]
    df4 <- df0
    for (i in 1:nrow(t4entries)){
      these.snps <- unlist(strsplit(t4entries$VARIANT_ID[i], ";"))
      these.distances <- rep(0, length(these.snps))
      these.mafs <- unlist(strsplit(t4entries$AF[i], ";"))
      df4 <- rbind(df4, data.frame(snps=as.character(these.snps), 
                                   distances=as.numeric(these.distances), 
                                   mafs=as.numeric(these.mafs)))
    }
  df <- rbind(df, df4)
  } 
  if(probe %in% tab5$PROBE){
    t5entries <- tab5[tab5$PROBE==probe,]
    df5 <- df0
    for (i in 1:nrow(t5entries)){
      these.snps <- unlist(strsplit(t5entries$VARIANT_ID[i], ";"))
      these.distances <- rep(-1, length(these.snps))
      these.mafs <- unlist(strsplit(t5entries$AF[i], ";"))
      df5 <- rbind(df5, data.frame(snps=as.character(these.snps), 
                                   distances=as.numeric(these.distances), 
                                   mafs=as.numeric(these.mafs)))
    }
    df <- rbind(df, df5)
  }
  if(probe %in% tab6$PROBE){
    t6entries <- tab6[tab6$PROBE==probe,]
    df6 <- df0
    for (i in 1:nrow(t6entries)){
      these.snps <- unlist(strsplit(t6entries$VARIANT_ID[i], ";"))
      probecoord <- IRanges(t6entries$MAPINFO[i], t6entries$MAPINFO[i])
      varcoord <- IRanges(t6entries$VARIANT_START[i], t6entries$VARIANT_END[i])
      these.distances <- rep(distance(probecoord, varcoord), length(these.snps))
      these.mafs <- unlist(strsplit(t6entries$AF[i], ";"))
      df6 <- rbind(df6, data.frame(snps=as.character(these.snps), 
                                   distances=as.numeric(these.distances), 
                                   mafs=as.numeric(these.mafs)))
    }
    df <- rbind(df, df6)
  }
  cbind(probe=rep(probe, nrow(df)), df)
}
probeinf <- lapply(allprobes, getallprobeinf)
epicsnps <- rbind.fill(probeinf)

#450K
tmp1 <- tempfile(fileext = ".xlsx")
download.file(url = "http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx", destfile = tmp1)
polyandsbe <- suppressWarnings(read_xlsx(tmp1, sheet = 1))
polyandsbe <- polyandsbe[!polyandsbe$PROBE %in% epicsnps$probe,]
withinprobe <- suppressWarnings(read_xlsx(tmp1, sheet = 2))
withinprobe <- withinprobe[!withinprobe$PROBE %in% epicsnps$probe,]

overcpg <- polyandsbe[grep("MAPINFO", polyandsbe$BASE_FROM_MAPINFO),]
overcpg <- data.frame(probe=overcpg$PROBE, snps=overcpg$SNP_ID, 
                      distances=0, mafs=overcpg$AF)
sbe <- polyandsbe[grep("SBE", polyandsbe$BASE_FROM_MAPINFO),]
sbe <- data.frame(probe=sbe$PROBE, snps=sbe$SNP_ID, distances=-1, mafs=sbe$AF)

withinprobe <- data.frame(probe=withinprobe$PROBE, snps=withinprobe$SNP_ID,
                          distances=distance(IRanges(withinprobe$MAPINFO, withinprobe$MAPINFO),
                                             IRanges(withinprobe$SNP_POS, withinprobe$SNP_POS)),
                          mafs=withinprobe$AF)

snpsall <- rbind(epicsnps, overcpg, sbe, withinprobe)

#XY.probes

locs450k <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
locsEPIC <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
XY.probes <- union(rownames(locs450k[locs450k$chr %in% c("chrX", "chrY"),]),
                  rownames(locsEPIC[locsEPIC$chr %in% c("chrX", "chrY"),]))


#Make tracks for hg19, hg38, mm10

#hg38

hg38 <- makeTxDbFromEnsembl(organism="Homo sapiens", release=102)
hg38 <- keepSeqlevels(hg38, as.character(c(1:22, "X", "Y", "MT")))
newStyle <- mapSeqlevels(seqlevels(hg38),"UCSC")
hg38 <- renameSeqlevels(hg38, newStyle)
hg38.grt <- GeneRegionTrack(hg38, genome="hg38", shape="box", fill = "lightblue", 
                            name = "Gene", showId = TRUE, geneSymbol = TRUE, 
                            transcriptAnnotation = "symbol", just.group="above")
gtf <- as.data.frame(rtracklayer::import("ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.chr.gtf.gz"))
m <- match(transcript(hg38.grt), gtf$transcript_id)
symbol(hg38.grt) <- gtf$gene_name[m]
hg38.generanges <- sapply(unique(gene(hg38.grt)), function (x) {these.ranges <- ranges(hg38.grt)[ranges(hg38.grt)$gene==x]
                                                                  GRanges(unique(seqnames(these.ranges)), 
                                                                          IRanges(sapply(unique(seqnames(these.ranges)), function (y) min(start(these.ranges[seqnames(these.ranges)==y]))), 
                                                                                  sapply(unique(seqnames(these.ranges)), function (y) max(end(these.ranges[seqnames(these.ranges)==y])))),
                                                                          strand=strand(these.ranges)[1], gene=x)})
hg38.generanges <- unlist(GRangesList(hg38.generanges))
m <- match(hg38.generanges$gene, gtf$gene_id)
hg38.generanges$symbol <- gtf$gene_name[m]

#hg19
hg19 <- makeTxDbFromEnsembl(organism="Homo sapiens", release=75)
hg19 <- keepSeqlevels(hg19, as.character(c(1:22, "X", "Y", "MT")))
newStyle <- mapSeqlevels(seqlevels(hg19),"UCSC")
hg19 <- renameSeqlevels(hg19, newStyle)
hg19.grt <- GeneRegionTrack(hg19, genome="hg19", shape="box", fill = "lightblue", 
                            name = "Gene", showId = TRUE, geneSymbol = TRUE, 
                            transcriptAnnotation = "symbol", just.group="above")
gtf <- as.data.frame(rtracklayer::import("ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"))
m <- match(transcript(hg19.grt), gtf$transcript_id)
symbol(hg19.grt) <- gtf$gene_name[m]
hg19.generanges <- sapply(unique(gene(hg19.grt)), function (x) {these.ranges <- ranges(hg19.grt)[ranges(hg19.grt)$gene==x]
                                                                  GRanges(unique(seqnames(these.ranges)), 
                                                                          IRanges(sapply(unique(seqnames(these.ranges)), function (y) min(start(these.ranges[seqnames(these.ranges)==y]))), 
                                                                                  sapply(unique(seqnames(these.ranges)), function (y) max(end(these.ranges[seqnames(these.ranges)==y])))),
                                                                          strand=strand(these.ranges)[1], gene=x)})
hg19.generanges <- unlist(GRangesList(hg19.generanges))
m <- match(hg19.generanges$gene, gtf$gene_id)
hg19.generanges$symbol <- gtf$gene_name[m]


#mm10
mm10 <- makeTxDbFromEnsembl(organism="Mus musculus", release=102)
mm10 <- keepSeqlevels(mm10, as.character(c(1:19, "X", "Y", "MT")))
newStyle <- mapSeqlevels(seqlevels(mm10),"UCSC")
mm10 <- renameSeqlevels(mm10, newStyle)
mm10.grt <- GeneRegionTrack(mm10, genome="mm10", shape="box", fill = "lightblue", 
                            name = "Gene", showId = TRUE, geneSymbol = TRUE, 
                            transcriptAnnotation = "symbol", just.group="above")
gtf <- as.data.frame(rtracklayer::import("ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.chr.gtf.gz"))
m <- match(transcript(mm10.grt), gtf$transcript_id)
symbol(mm10.grt) <- gtf$gene_name[m]
mm10.generanges <- sapply(unique(gene(mm10.grt)), function (x) {these.ranges <- ranges(mm10.grt)[ranges(mm10.grt)$gene==x]
                                                                  GRanges(unique(seqnames(these.ranges)), 
                                                                          IRanges(sapply(unique(seqnames(these.ranges)), function (y) min(start(these.ranges[seqnames(these.ranges)==y]))), 
                                                                                  sapply(unique(seqnames(these.ranges)), function (y) max(end(these.ranges[seqnames(these.ranges)==y])))),
                                                                          strand=strand(these.ranges)[1], gene=x)})
mm10.generanges <- unlist(GRangesList(mm10.generanges))
m <- match(mm10.generanges$gene, gtf$gene_id)
mm10.generanges$symbol <- gtf$gene_name[m]

###EPICv2 SNPs from manifest
#Additional File 4 from Peters et al. 2024
EPICv2manifest <- read.csv("~/analysis_tmp/EPICv2/AdditionalFile4.csv")
rownames(EPICv2manifest) <- EPICv2manifest$Probe_ID


getallprobeinf <- function(probe){
  probeentries <- EPICv2manifest[probe,]
  df <- data.frame(snps=unlist(strsplit(probeentries$SNP_ID, ";")),
                   distances=as.numeric(unlist(strsplit(probeentries$SNP_DISTANCE, ";"))),
                   mafs=as.numeric(unlist(strsplit(probeentries$SNP_MinorAlleleFrequency, ";"))))
  df$probe <- rep(probe, nrow(df))
  df <- df[,c(2:4, 1)]
  df
}
probeinf <- lapply(rownames(EPICv2manifest)[nchar(EPICv2manifest$SNP_ID) > 0], getallprobeinf)
epicv2snps <- rbind.fill(probeinf)

#EPICv2 data from Noguera-Castells et al. (2023)

system("wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222919/suppl/GSE222919%5Fprocessed%5Fdata.txt.gz")
esteller <- read.table(gzfile("GSE222919_processed_data.txt.gz"), sep = "\t", row.names = 1, header = T)
ALLbetas <- esteller[3:nrow(esteller), grep("BALL_0|TALL_0", esteller[2,])]
detPs <- esteller[3:nrow(esteller), grep("BALL_0|TALL_0", esteller[2,]) + 1]
colnames(ALLbetas) <- colnames(detPs) <- esteller[2, grep("BALL_0|TALL_0", esteller[2,])]
ALLbetas <- apply(ALLbetas, 2, as.numeric)
detPs <- apply(detPs, 2, as.numeric)
rownames(ALLbetas) <- rownames(detPs) <- rownames(esteller)[3:nrow(esteller)]

#Remove detPs > 0.05
rm <- apply(detPs, 1, function (x) any(x > 0.05))
table(rm)
# FALSE   TRUE 
#894902    926 

ALLbetas <- ALLbetas[!rm,]
ALLbetas <- data.matrix(ALLbetas)

#Offset
ALLbetas[ALLbetas==0] <- 0.001
ALLbetas[ALLbetas==1] <- 0.999

save(hg19.generanges, file="hg19.generanges.Rda")
save(hg38.generanges, file="hg38.generanges.Rda")
save(mm10.generanges, file="mm10.generanges.Rda")
save(mm10.grt, file="mm10.grt.Rda")
save(hg19.grt, file="hg19.grt.Rda")
save(hg38.grt, file="hg38.grt.Rda")
save(snpsall, file="snpsall.Rda")
save(crosshyb, file="crosshyb.Rda")
save(XY.probes, file="XY.probes.Rda")
save(epicv2snps, file="epicv2snps.Rda")
save(ALLbetas, file="ALLbetas.Rda")
