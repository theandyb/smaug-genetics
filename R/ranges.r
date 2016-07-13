library(devtools)
install_github('carjed/bedr')
library(bedr)
library(dplyr)

sites<-read.table("/net/bipolar/jedidiah/mutation/output/logmod_data/motifs/GC_TA/dp/GC_TA_GTCCTGT_dp.txt", header=F)

# Get inside/outside broad histone peak status for each site
histCol <- function(sites, mark){
  file <- paste0("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", mark, ".bed")
  feat_ranges <- bed_to_granges(file, header=F)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$V2), ranges=IRanges(start=sites$V1, end=sites$V1))
  return(as.integer(site_ranges %within% feat_ranges))
}

# Get CpG Island status of each site
cpgiCol <- function(sites){
  file <- "/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.bed"
  feat_ranges <- bed_to_granges(file, header=F)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$V2), ranges=IRanges(start=sites$V1, end=sites$V1))
  return(as.integer(site_ranges %within% feat_ranges))
}

# Get recombination rate at each site
rcrCol <- function(sites){
  file <- "/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed"
  feat_ranges <- bed_to_granges(file, header=T)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$V2), ranges=IRanges(start=sites$V1, end=sites$V1))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")

  indices[is.na(indices)]<-0
  ind_df<-data.frame(V1=sites$V1, V2=sites$V2, indices)

  feat_df<-as.data.frame(feat_ranges)
  feat_df$indices<-seq_along(1:nrow(feat_df))
  rates <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(V2, V1) %>%
    select(id)

  rates[is.na(rates)]<-0
  return(as.vector(rates))
}

# Loop to add histone marks to site data
hists<-c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")

dflist <- list()
for(i in 1:length(hists)){
  mark <- hists[i]
  hist <- histCol(sites, mark)
  colname <- make.names(mark)
  dflist[[i]] <- hist
}

df <- as.data.frame(do.call(cbind, dflist))
names(df) <- hists
sites<-cbind(sites, df)

# Add CpG Island status
sites$CpGI<-cpgiCol(sites)

# Add recombination rate
sites$RR<-rcrCol(sites)