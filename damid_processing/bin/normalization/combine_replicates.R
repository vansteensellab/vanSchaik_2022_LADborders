#!/usr/bin/Rscript

# Combine replicates (version 1.0)
## Tom van Schaik, 180116

# - Functionality:
# Simply, combine multiple replicates. This is simply the mean log2 signal.
# This script does not deserve more comments or a logger

# - Input
# (required)
# * comma-seperated replicates
# * output directory
#
# (optional)
# * basename - basename of the output files (default: basename from file)

# - Output (within output_dir)
# * [BASENAME].norm.txt - bed-like output format

# - Method
# Simply, mean log2

# - Version
# 1.0 - First version

# - To do
# 

suppressPackageStartupMessages(library("optparse"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(replicates, output_dir, basename=NULL) {
    # Simply, read in all the replicates and make a mean log2 table
    
    # Setting-up
    replicates <- strsplit(replicates, ",")[[1]]
    
    if (is.null(basename)) {
        basename <- strsplit(basename(replicates[1]), "\\.")[[1]][1]
    }
    
    # Read in the normalized data
    norm.tmp <- read.table(replicates[1], sep = "\t")
    
    if (length(replicates) > 1) {
        for (i in 2:length(replicates)) {
            norm.tmp <- cbind(norm.tmp, read.table(replicates[i], sep = "\t")[, 4])
        }
    }
    
    # Create the mean log2 score, where NA is ignored
    norm <- norm.tmp[, 1:3]
    if (length(replicates) > 1) {
        norm[, 4] <- rowMeans(norm.tmp[, 4:ncol(norm.tmp)], na.rm = T)
    }
    
    # Write the table
    dir.create(output_dir, showWarnings = FALSE)
    output_file <- file.path(output_dir,
                             paste0(basename, "-combined.norm.txt.gz"))
	output_file <- gzfile(output_file, "w")
	write.table(norm, file = output_file, quote = FALSE, sep = "\t", col.names = FALSE,
    			row.names = F)
	close(output_file)
    
}

##############################################################################
## Supporting functions ######################################################
##############################################################################



##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.0"

# Get input
option_list <- list(
    make_option(c("-r", "--replicates"), type="character", default=NULL, 
                help="replicates, comma-separated (required)", metavar="character"),
    make_option(c("-o", "--output_dir"), type="character", default=NULL,
                help="output directory (required)", metavar="character"),
                
    # Optional arguments            
    make_option(c("-n", "--basename"), type="character", default=NULL, 
                help="basename (default: basename first replicate)", metavar="character")     
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$replicates) | is.null(opt$output_dir)) {
  	print_help(opt_parser)
  	stop("Not all required arguments present.\n", call.=FALSE)
} else {
	main(opt$replicates, opt$output_dir, opt$basename)
}
