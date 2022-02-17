#!/opt/R/4.0.5/bin/Rscript

# Pipeline basic statistics
## Tom van Schaik, 171030

# - Functionality:
# Merge all statistics documents into one table + figures

# - Input
# (required)
# * Basenames
# * Output directory
#

# - Output
# * pipeline.statistics.txt

# - Method
# Combine all the individual statistics files into one document

# - Version
# 1.0 - Initial version
# 1.1 - Moved from reshape2 to tidyverse

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(basenames, output_dir, input_directory) {
  
  library(ggplot2)
  #library(reshape2)
  library(RColorBrewer)
  library(tidyverse)
  
  # Some setting-up
  dir.create(output_dir, showWarnings = FALSE)
  output <- file.path(output_dir, "pipeline.statistics.txt")
  
  # Create list of basenames
  basenames <- strsplit(basenames, ",")[[1]]

  stopifnot(! any(duplicated(basenames)))
  
  # For each basename, get all statistics
  statistics.df <- data.frame()
  for (basename in basenames) {
    statistics <- read.table(file.path(input_directory, paste0(basename, ".statistics.txt")),
                             header = T, sep = "\t", stringsAsFactors = FALSE)
    
    statistics.df <- rbind(statistics.df, statistics)
  }
  
  # Write output table
  write.table(statistics.df, output, sep = "\t", row.names = F, quote = F)
  
  # Plot this as well
  #x <- melt(statistics.df, id.vars = "basename")
  x <- as_tibble(statistics.df) %>%
    gather(key, value, -basename)
  x$basename <- factor(x$basename, levels = statistics.df$basename)
  
  plt <- ggplot(x, aes(x = basename, y = value, fill = key)) + 
    geom_bar(stat = "identity", position = "dodge") +
    ylab("Count") +
    xlab("Sample") +
    scale_fill_grey() +
    theme_bw()
  
  plot_name <- file.path(output_dir, "pipeline.statistics.pdf")
  
  pdf(plot_name, height = 4.5, width = 3 + nrow(statistics.df))
  plot(plt)
  dev.off()
  
}

##############################################################################
## Supporting functions ######################################################
##############################################################################

# None

##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.0"

option_list <- list(
  make_option(c("-b", "--basenames"), type="character", default=NULL, 
              help="basename list, comma separated (required)", 
              metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="output directory (required)", metavar="character"),
  
  # Optional arguments
  make_option(c("-d", "--input_directory"), type="character", default="results/statistics", 
              help=paste0("directory where statistics can be found ", 
                          "(default: 'results/statistics')"),
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$basenames) | is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Not all required arguments given\n", call.=FALSE)
} else {
  main(opt$basenames, opt$output_dir, opt$input_directory)
}
