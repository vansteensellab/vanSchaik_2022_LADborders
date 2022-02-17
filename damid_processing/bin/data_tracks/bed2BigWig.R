#!/usr/bin/Rscript

# Bed2BigWig
## Tom van Schaik, 171030

# - Functionality:
# Take a bed-like file as input and use the 4th column to make a bigwig

# - Input
# (required)
# * Normalized file
# * Output directory
# * Chromosome sizes

# - Output
# * [basename].bw

# - Method
# Create a BigWig of the bed file

# - Version
# 1.0 - Initial version

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(normalized, output_dir, chrom_sizes, basename, bedGraphToBigWig,
                 cpm) {
  
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(rtracklayer))
  
  # Some setting-up
  if (is.null(basename)) {
    basename <- strsplit(basename(normalized), "\\.")[[1]][1]
  }
  
  dir.create(output_dir, showWarnings = FALSE)
  output_bw <- file.path(output_dir, paste0(basename, ".bw"))
  
  # Read the data
  normalized <- read.table(normalized,
                           sep = "\t", stringsAsFactors = FALSE, 
                           col.names = c("chr", "start", "end", "score"))
  
  # Note: the .txt file is 0-based (bed-like), while the GRanges is 1-based
  normalized$start <- normalized$start + 1
  
  # If overlapping bins, fix this
  if (normalized[1, 3] > normalized[2, 2]) {
    # So far, always 2 bases shift
    normalized <- FixGATCsequences(normalized)
  }
  
  # Remove NAs first
  normalized <- normalized[complete.cases(normalized), ]
  
  # Normalization (CPM)
  if (cpm) {
    normalized$score <- normalized$score / sum(normalized$score) * 1e6
  }

  # Make BigWig
  makeBigwig(as(normalized, "GRanges"),
             output_bw,
             chrom_sizes,
             bedGraphToBigWig)
  
}

##############################################################################
## Supporting functions ######################################################
##############################################################################

makeBigwig <- function(gr, file, chrom_sizes, bedGraphToBigWig) {
  # Use this stupid solution to create a BigWig file, as export.bw is crap
  # Assumes a GRanges object with a single column "score", and a file to
  # write to.
  
  # Use temp files
  x_tmp <- tempfile()
  x_tmp2 <- tempfile()
  
  # First create a bedgraph, sort this and then create a bigwig
  export.bedGraph(gr, con = x_tmp)
  system(paste0("sort -k1,1 -k2,2n ", x_tmp, " > ", x_tmp2))
  system(
    paste(
      bedGraphToBigWig,
      x_tmp2,
      chrom_sizes,
      file,
      sep = " "
    )
  )
  
  # Remove temp files
  file.remove(x_tmp)
  file.remove(x_tmp2)
}

FixGATCsequences <- function(df, bases = 2) {
  # The GATC fragments have the downside that they overlap, as they end with
  # GATC and also start with GATC. This function quickly fixes this issue.
  
  df.gr <- as(df, "GRanges")
  
  # Fix the start / end
  start(df.gr) <- start(df.gr) + bases
  end(df.gr) <- end(df.gr) - bases
  
  # Also trim the object to be within range
  df.gr <- trim(df.gr)
  
  as(df.gr, "data.frame")[, c("seqnames", "start", "end", "score")]
}

##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.0"

option_list <- list(
  make_option(c("-n", "--normalized"), type="character", default=NULL, 
              help="normalized file", 
              metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="output directory (required)", metavar="character"),
  make_option(c("-s", "--chrom_sizes"), type="character", default=NULL, 
              help="chromosome sizes (required)",
              metavar="character"),
  
  # Optional arguments
  make_option(c("-b", "--basename"), type="character", default=NULL, 
              help="basename (optional)",
              metavar="character"),
  make_option(c("-g", "--bedGraphToBigWig"), type="character", 
              default="/home/t.v.schaik/mydata/bin/bedGraphToBigWig", 
              help="bedGraphToBigWig (default: Tom's location)",
              metavar="character"),
  make_option(c("-c", "--cpm"), type="logical", default=FALSE, action="store_true",
              help="CPM normalization for counts (optional)",
              metavar="character")
)



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$normalized) | is.null(opt$output_dir) | is.null(opt$chrom_sizes)) {
  print_help(opt_parser)
  stop("Not all required arguments given\n", call.=FALSE)
} else {
  main(opt$normalized, opt$output_dir, opt$chrom_sizes, opt$basename, 
       opt$bedGraphToBigWig, opt$cpm)
}
