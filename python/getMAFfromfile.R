#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The core functions and methodolody in this R code has been kindly shared by Alexander Nathaniel Gorelick, PhD from Massachusetts General Hospital (Harvard)
# This is R code to generate the MAF file based on the above 
# workflow. 
#
# - In R, load the new tumor and normal MAFs, format some 
#.  important fields, and merge the result
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inputs
args = commandArgs(trailingOnly=TRUE)
tumormaf <- args[1]
normalmaf <- args[2]
outfile <- args[3]

normed_filtered_tumor_maf_file <- tumormaf
normed_filtered_normal_maf_file <- normalmaf

library(data.table)
load_and_correct_null_maf <- function(file) {
  ## vcfs with 0 mutations cause problems
  ## here we load the maf, check if it has 0 mutation, and then re-load it with appropriate arguments
  maf <- fread(file)
  if( grepl('version',names(maf)[1]) & nrow(maf)==1 ) maf <- fread(file,skip=1)
  
  maf$ShortVariantID <- paste0(maf$Reference_Allele,maf$Start_Position,maf$Tumor_Seq_Allele2)
  maf$t_ref_count <- as.integer(maf$t_ref_count)
  maf$t_alt_count <- as.integer(maf$t_alt_count)
  maf$t_depth <- as.integer(maf$t_depth)
  maf$TumorVAF <- maf$t_alt_count / maf$t_depth
  maf[is.na(t_ref_count),t_ref_count:=0]
  maf[is.na(t_alt_count),t_alt_count:=0]
  maf[is.na(t_depth),t_depth:=0]
  maf[is.na(TumorVAF),TumorVAF:=0]
  
  fsplit <- function(x) {
    x[x=='.'] <- '0,0'
    x <- as.integer(strsplit(x,',')[[1]])
    x <- x[!is.na(x)]
    if(length(x)==2) {
      out <- as.list(x)
    } else {
      x <- rep(as.integer(NA),2)
      out <- as.list(x)
    }
    names(out) <- c('ref','alt')
    out
  }
  F1R2 <- rbindlist(lapply(maf$t_F1R2, fsplit))
  F2R1 <- rbindlist(lapply(maf$t_F2R1, fsplit))
  
  maf$n_ref_count <- as.integer(NA)
  maf$n_alt_count <- as.integer(NA)
  maf$n_depth <- as.integer(NA)
  maf$NormalVAF <- as.numeric(NA)
  maf
}

## load/format the tumor MAF
maf <- load_and_correct_null_maf(normed_filtered_tumor_maf_file)

## load/format the normal MAF
ref_maf <- load_and_correct_null_maf(normed_filtered_normal_maf_file)

## merge fields with normal allele read support (from the normal MAF) to the tumor MAF
ref_maf <- ref_maf[,c('ShortVariantID','t_alt_count','t_ref_count','t_depth','TumorVAF'),with=F]
names(ref_maf) <- c('ShortVariantID','n_alt_count','n_ref_count','n_depth','NormalVAF')
maf[,c('n_alt_count','n_ref_count','n_depth','NormalVAF'):=NULL]
maf <- merge(maf, ref_maf, by='ShortVariantID', all.x=T)
maf[is.na(n_alt_count),n_alt_count:=0]
maf[is.na(n_ref_count),n_ref_count:=0]
maf[is.na(n_depth),n_depth:=0]
maf[is.na(NormalVAF),NormalVAF:=NA]

## save the resulting maf data.table, this should be the complete MAF with variants called in the tumor, annotated with their read support in the normal bam.
write.table(maf, outfile, quote=FALSE, sep="\t")