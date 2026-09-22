#!/usr/bin/env Rscript

# Script to calculate median bin coverage per sample for all samples in a bincov
# matrix

# Adapted from GATK-SV's own src/WGD/bin/medianCoverage.R -- see
# bin/README.md's own note on this file for why and how. Same CLI,
# options, and output format as upstream; the internal computation was
# rewritten (data.table::fread + matrixStats column/row medians on a
# plain numeric matrix, instead of read.table + apply() over repeated
# as.data.frame() copies) because the original's memory use scaled far
# worse than the matrix's own size as this pipeline's cohorts grew --
# see bin/README.md for the full explanation.

# Load libraries
require(optparse)
require(data.table)
require(matrixStats)

# Define options
option_list <- list(
  make_option(c("-b","--binwise"), action="store_true", default=FALSE,
              help="compute medians of all samples per bin [default: median of all bins per sample]"),
  make_option(c("-m","--mad"), action="store_true", default=FALSE,
              help="compute median absolute deviation of all bins per sample [default: FALSE]"),
  make_option(c("-H","--header"), action="store_true", default=FALSE,
              help="input coverage matrix has header with sample IDs [default: FALSE]"))

# Get command-line arguments and options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2)
{cat("Incorrect number of required positional arguments\n\n")
  stop()}

# Read matrix. fread() in place of read.table(): same positional-column
# semantics (no rownames), but streams the file into a data.table rather
# than read.table's incremental per-column type-guessing/growing, which
# is both slower and far more memory-hungry on a matrix with hundreds of
# sample columns. comment.char="" (header mode) / "#" (no-header mode)
# preserved exactly as upstream, for the same reason upstream sets it:
# medianCoverage.R's own header-mode output starts with "#sample_id",
# which fread's default comment handling would otherwise swallow as a
# comment line the same way read.table's would.
if(opts$header==T){
  cov <- fread(args$args[1], header=T, sep="\t", check.names=F, comment.char="")
}else{
  cov <- fread(args$args[1], header=F, sep="\t", check.names=F, comment.char="#")
}

# Function to compute medians per sample
covPerSample <- function(cov,downsample=1000000,mad=F){
  # Downsample to 1M random rows if nrows > 1M (for computational efficiency)
  if(nrow(cov)>1000000){
    cov <- cov[sample(1:nrow(cov), downsample),]
  }
  # Single numeric matrix, built once, instead of upstream's repeated
  # as.data.frame(cov[,-c(1:3)]) at every one of the three median/mad
  # passes below -- each of those was a fresh full-matrix copy plus
  # apply()'s own per-row coercion overhead (apply() over a data.frame is
  # markedly worse than over a matrix, since a data.frame is internally a
  # list of columns that apply() has to rebind row-by-row). One matrix
  # here, reused; rowMedians/colMedians run in C over it directly.
  m <- as.matrix(cov[, -(1:3)])
  # Get medians with and without zero-cov bins
  bin_med <- rowMedians(m, na.rm=T)
  zerobins <- which(as.integer(bin_med) == 0)
  withzeros <- colMedians(m, na.rm=T)
  withoutzeros <- colMedians(m[-zerobins, , drop=FALSE], na.rm=T)
  #Get SDs with and without zero-cov bins (if optioned)
  if(mad==T){
    withzeros.mad <- colMads(m, na.rm=T)
    withoutzeros.mad <- colMads(m[-zerobins, , drop=FALSE], na.rm=T)
  }
  # Compile results df to return
  if(mad==T){
    res <- data.frame("ID"=paste("Sample",1:(ncol(cov)-3),sep=""),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros,
                      "MAD_withZeros"=withzeros.mad,
                      "MAD_withoutZeros"=withoutzeros.mad)
  }else{
    res <- data.frame("ID"=paste("Sample",1:(ncol(cov)-3),sep=""),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros)
  }
  # Replace sample IDs if input matrix has header
  if(opts$header==T){
    res$ID <- names(cov)[-(1:3)]
  }
  # Return output df
  return(res)
}

# Function to compute medians per bin
covPerBin <- function(cov,downsample=500,mad=F){
  # Downsample to 500 random samples if nsamples > 500 (for computational efficiency)
  if(ncol(cov)>503){
    cov <- cov[,sample(1:ncol(cov), downsample)]
  }
  # Same matrix-once, C-level-ops approach as covPerSample above, applied
  # row-wise (per bin, across samples) instead of column-wise.
  m <- as.matrix(cov[, -(1:3)])
  withzeros <- rowMedians(m, na.rm=T)
  withoutzeros <- vapply(seq_len(nrow(m)), function(i){
    vals <- m[i,]
    pos <- vals[vals>0]
    if(length(pos)>0) median(pos, na.rm=T) else NA_real_
  }, numeric(1))
  # Get standard deviations (if optioned)
  if(mad==T){
    withzeros.mad <- rowMads(m, na.rm=T)
    withoutzeros.mad <- vapply(seq_len(nrow(m)), function(i){
      vals <- m[i,]
      pos <- vals[vals>0]
      if(length(pos)>0) mad(pos, na.rm=T) else NA_real_
    }, numeric(1))
  }
  # compile results df to return
  if(mad==T){
    res <- data.frame("#chr"=cov[[1]],"start"=cov[[2]],"end"=cov[[3]],
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros,
                      "MAD_withZeros"=withzeros.mad,
                      "MAD_withoutZeros"=withoutzeros.mad,
                      check.names=FALSE)
  }else{
    res <- data.frame("#chr"=cov[[1]],"start"=cov[[2]],"end"=cov[[3]],
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros,
                      check.names=FALSE)
  }
  # Return output df
  return(res)
}

# Compute appropriate medians & write out
if(opts$binwise==TRUE){
  res <- covPerBin(cov,mad=opts$mad)
  names(res)[1] <- "#chr"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}else{
  res <- covPerSample(cov,mad=opts$mad)
  names(res)[1] <- "#sample_id"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}
