#!/usr/bin/Rscript

# HMM
## Tom van Schaik, 171030

# - Functionality:
# Define enriched domains with a HMM 

# - Input
# (required)
# * Normalized file
# * Output directory
# * centromeres location

# - Output
# * [basename]_HMM.bed
# * [basename]_AD.bed

# - Method
# Run HMM on the desired file

# - Version
# 1.0 - Initial version

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(normalized, output_dir, basename, centromeres) {
  
  suppressPackageStartupMessages(library(HMMt))
  suppressPackageStartupMessages(library(GenomicRanges))
  
  # Some setting-up
  if (is.null(basename)) {
    basename <- strsplit(basename(normalized), "\\.")[[1]][1]
  }
  
  dir.create(output_dir, showWarnings = FALSE)
  output_hmm <- file.path(output_dir, paste0(basename, "_HMM.txt.gz"))
  output_ad <- file.path(output_dir, paste0(basename, "_AD.bed.gz"))
  
  # Read-in the data
  normalized_data <- read.table(normalized, sep = "\t", stringsAsFactors = FALSE,
                                col.names = c("chr", "start", "stop", "score"))

  if (! is.null(centromeres)) {
    centromeres_data <- read.table(centromeres, sep = "\t", stringsAsFactors = FALSE)
    centromeres_data <- centromeres_data[, 2:4]
    names(centromeres_data) <- c("chr", "start", "stop")
  
    # The centromeres contain multiple regions per chromosome - get the range of this
    centromeres_data <- cbind(unique(centromeres_data$chr),
                              do.call(rbind, lapply(unique(centromeres_data$chr), 
                                      function(x) range(centromeres_data[centromeres_data$chr == x, c("start", "stop")]))))
    centromeres_data <- data.frame(centromeres_data)
    names(centromeres_data) <- c("chr", "start", "stop")
  }

  # Run the HMM and select "associated domains"
  hmm_calls <- HMM(normalized_data, na_solution = "keep")

  hmm_ranges <- getHMMRanges(as(hmm_calls, "GRanges"), score = "AD")
  
  if (! is.null(centromeres)) {
    hmm_ranges <- setdiff(hmm_ranges, as(centromeres_data, "GRanges"))
  }
  
  # Convert the unknown bins to NA for the HMM .txt file, but not the .bed file
  hmm_calls$model[is.na(normalized_data$score)] <- NA

  # Write data tables - gzip'd
  out <- gzfile(output_hmm, "w")
  write.table(hmm_calls, 
              file = out, 
              quote = F, sep = "\t", row.names = F, col.names = F)
  close(out)
  
  
  out <- gzfile(output_ad, "w")
  write.table(as(hmm_ranges, "data.frame")[, 1:3],
              file = out, 
              quote = F, sep = "\t", row.names = F, col.names = F)
  close(out)
}

##############################################################################
## Supporting functions ######################################################
##############################################################################

# HMM functions
HMM <- function(normalized, na_solution) {
  # Run a basic HMM for each data column of a damid data frame.
  if (!na_solution %in% c("NA", "keep", "-")) {
    stop("Unknown na_solution")
  }
  
  # Add index to easily convert to the original order
  normalized.tmp <- normalized
  normalized.tmp$idx <- 1:nrow(normalized)
  
  normalized.tmp <- normalized.tmp[order(normalized.tmp$chr, normalized.tmp$start), ]
  
  # HMM
  br <- bridge(normalized.tmp[, 1:4])
  
  # Flush the undesired output to the sink
  sink(file="/dev/null")
  fit <- BaumWelchT(x=br$x, series.length=br$series.length)
  sink()
  
  boundstate <- which.max(fit$mu)
  model <- 0 + (fit$ViterbiPath[br$nonvirtuals] == boundstate)
  model <- ifelse(model == 1, "AD", "iAD")
  
  df.hmm <- cbind(normalized.tmp[, c(1:3, 5)], model)
  
  # Finally, I want to make the assumption that bins with NA cannot be called 
  # LAD / ...
  if (na_solution == "NA") {
    df.hmm$model[is.na(normalized.tmp$score)] <- NA
  } else if (na_solution == "-") {
    df.hmm$model[is.na(normalized.tmp$score)] <- "-"
  }
  
  # Convert into the original order
  df.hmm <- df.hmm[order(df.hmm$idx), ]
  
  df.hmm[, c(1:3, 5)]
}
getHMMRanges <- function(gr, i = 1, score = 1) {
  # Get GRanges bins of a GRanges object with a binary score column
  
  # First, remove sequences with NA
  gr <- gr[! is.na(mcols(gr)[, i])]
  
  gr <- gr[mcols(gr)[, i] == score]
  gr <- reduce(gr)

  gr
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
  
  # Optional arguments
  make_option(c("-b", "--basename"), type="character", default=NULL, 
              help="basename (optional)",
              metavar="character"),
  make_option(c("-c", "--centromeres"), type="character", default=NULL, 
              help="centromeres to mask from bed file (optional)",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$normalized) | is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Not all required arguments given\n", call.=FALSE)
} else {
  main(opt$normalized, opt$output_dir, opt$basename, opt$centromeres)
}
