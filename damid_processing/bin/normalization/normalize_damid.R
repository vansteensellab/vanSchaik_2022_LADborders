#!/usr/bin/Rscript

# Normalize DamID (version 1.2)
## Tom van Schaik, 160831

# - Functionality:
# Normalize a DamID counts table with a given method. Several options are 
# possible, depending on available data sets and preference. More details
# are listed below.

# - Input
# (required)
# * counts table
# * method (see below)
# * output directory
#
# (optional)
# * basename - basename for output files (default: input basename)
# * dam-only - dam-only file for various methods (default: NULL)
# * pseudo - pseudocount to prevent unallowed division (default: 1)
# * scaling-factor - scaling factor for the dam-substr method (default: 1)
# * stdout - report to stdout (default: false)
# * quiet - be quiet (default: false)

# - Output (within output_dir)
# * [BASENAME].norm.txt - bed-like output format
# * [BASENAME].norm.statistics.txt - accompanying statistics

# - Method
# Read in counts table and perform a normalization method. Possibilities:
# - reads:          normalize to reads / M
# - dam-log2:       normalize to reads / M, determine log2-ratio over dam-only
#                   with pseudocounts
# - dam-substr:     normalize to reads / M, substract target - dam-only * scaling
#                   factor, divided by the root of the sum of both (+ pseudo)
# - OE:             observed over expected (not implemented yet)

# - Version
# 1.2

# - History
# * 1.01 - added norm-substr method
# * 1.1 - no reads == NA
# * 1.2 - <10 read == NA

# - To do
# * Add OE support
# * Create some figures
# * Better statistics file
# * Test normalization methods

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tools"))

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(counts, method, output_dir, basename, dam_only=NULL, pseudo=1,
                 scaling=1, stdout=FALSE, zip=FALSE, min_reads=10) {
	# Main function
    	 
	# Some sanity checks
    logger(paste0("----- normalize-damid version: ", script.version, " -----"))
    logger("Setting up...")
	stopifnot(file.exists(counts))
    
    stopifnot(method %in% c("reads", "dam-log2", "dam-substr", "OE"))
    
	if (method %in% c("dam-log2", "dam-substr")) {
		stopifnot(!is.null(dam_only))
        stopifnot(file.exists(dam_only))
	}
    
    # Setting-up
    if (is.null(basename)) {
        # Basename = remove directory + everything after the first "." [dot]
        basename <- sub("\\..*", "", basename(counts))
    }

    dir.create(output_dir, showWarnings = FALSE)
    
    if (! stdout) {
        out <- file.path(output_dir, paste0(basename, ".norm.txt"))
		if (zip) {
			out <- paste0(out, ".gz")
		}
    } else {
        out <- ""
    }
    out_stat <- file.path(output_dir, paste0(basename, ".norm.statistics.txt"))
    
    logger(paste0("\tcount file: ", counts))
    logger(paste0("\tmethod: ", method))
    logger(paste0("\toutput file: ", out, " (Note: [nothing] = stdout)"))
    logger(paste0("\tstatistics file: ", out_stat))
    
    # Read counts table
    logger("\n----------")
    logger(paste0("Reading counts ", counts, "..."))
    counts <- read_counts(counts)
    logger(paste0(nrow(counts), " rows read"))
    
    # Normalize, depending on method
    logger(paste0("Normalizing with method ", method))
    if (method == "reads") {
        norm <- norm.reads(counts)
    } else if (method == "dam-log2") {
        norm <- norm.dam_log2(counts, dam_only, pseudo, min_reads)
    } else if (method == "dam-substr") {
        x <- norm.dam_substr(counts, dam_only, pseudo, scaling, output_dir, basename)
        norm <- x[[1]]
        scaling <- x[[2]]
    } else if (method == "OE") {
        norm <- norm.OE(counts, OE=NULL)
    } else {
        stop("Method not detected. Please look into this.")
    }
		
    # Write output table
    logger("\n----------")
    logger("Write output")
    
    # Prevent scientific notation
    options(scipen = 999)
	if (endsWith(out, ".gz")) {
		out <- gzfile(out, "w")
    	write.table(norm, file = out, quote = FALSE, sep = "\t", col.names = FALSE,
        			row.names = F)
		close(out)
	} else { 
    	write.table(norm, file = out, quote = FALSE, sep = "\t", col.names = FALSE,
        			row.names = F)
	}
	            
    # Create statistics
    logger("\n----------")
    logger(paste0("Total counts: ", sum(counts$counts)))
    logger(paste0("Mean normalization score: ", round(mean(na.omit(norm$score)), 4)))
    
    stat <- data.frame(n = c("reads", "norm_mean", "norm_range"),
                       c = c(nrow(counts), 
                             round(mean(na.omit(norm$score)), 4),
                             paste0(round(min(na.omit(norm$score)), 4), 
                                    ":", 
                                    round(max(na.omit(norm$score)), 4))))
                                    
    if (!is.null(scaling)) {
        stat <- rbind(stat, data.frame(n="scaling_factor", 
                                       c=as.character(round(scaling, 3))))
    }
    
    write.table(stat, out_stat, quote = FALSE, sep = "\t", col.names = FALSE,
                row.names = FALSE)
    # Add versions
    write(paste0("## Versions ##"), out_stat, append = TRUE)
    write(paste0("# R: ", version["version.string"][[1]]), out_stat, append = TRUE)
    write(paste0("# normalize-damid: ", script.version), out_stat, append = TRUE)
    
    
    # Report session info
    logger("\n----------")
    logger("Session info:")
    logger(capture.output(sessionInfo()))
    logger("")            
	
}

##############################################################################
## Normalization methods #####################################################
##############################################################################

# Here, we can list all normalization methods, which all take the same input
# format and return the same output format. (Both bed-like data frames). In
# between, you can do whatever you fancy.

norm.reads <- function(counts) {
    # This is a target-only reads / M normalization
    
    # Normalize for library size (in reads / M)
    counts$counts <- counts$counts / sum(counts$counts) * 1e6
    
    # Create norm table
    norm <- counts[, 1:3]
    norm$score <- counts$counts
    
    norm
}

norm.dam_log2 <- function(counts, dam_only, pseudo, min_reads) {
    # This is the dam-only log2 normalization
    
    # Read dam_only table
    logger(paste0("Reading counts ", dam_only, "..."))
    dam_only <- read_counts(dam_only)
    
    idx.too_few_reads <- which(counts$counts + dam_only$counts < min_reads)
    
    # Normalize for library size (in reads / M)
    counts$counts <- counts$counts / sum(counts$counts) * 1e6
    dam_only$counts <- dam_only$counts / sum(dam_only$counts) * 1e6
        
    # Add pseudocount for later divisions
    counts$counts <- counts$counts + pseudo
    dam_only$counts <- dam_only$counts + pseudo
    
    # Create normalized table
    norm <- counts[, 1:3]
    norm$score <- log2(counts$counts / dam_only$counts)
    
    # No reads = NA
    norm[idx.too_few_reads, "score"] <- NA
    
    norm
}

norm.dam_substr <- function(counts, dam_only, pseudo, scaling_factor=NULL,
                            output_dir, basename) {
    # This is the dam-only substration + total signal division-normalization
    #   Note that I haven't decided how valid this method is, and how I could 
    #   choose a reasonable scaling factor.
    
    # Read dam_only table
    logger(paste0("Reading counts ", dam_only, "..."))
    dam_only <- read_counts(dam_only)
    
    # Normalize for library size (in reads / M)
    counts$counts <- counts$counts / sum(counts$counts) * 1e6
    dam_only$counts <- dam_only$counts / sum(dam_only$counts) * 1e6
    
    # No scaling - testing
    norm.noscaling <- (counts$counts / (dam_only$counts+pseudo))
    idx.inf <- which(is.infinite(norm.noscaling))
    norm.noscaling[idx.inf] <- counts$counts[idx.inf]
    
    # If not given, determine scaling factor at 10% of the data
    if (is.null(scaling_factor)) {
        scaling_factor <- quantile(na.omit(norm.noscaling), 0.1)[[1]]
    } 
        
    # Normalize using scaling factor
    norm.counts <- counts$counts - dam_only$counts * scaling_factor
    norm.counts <- (norm.counts / (dam_only$counts+pseudo))
    
    if (pseudo == 0) {
        # Inf -> dam-only == 0, replace with original value
        idx.inf <- which(is.infinite(norm.counts))
        norm.counts[idx.inf] <- counts$counts[idx.inf]
        # NaN -> both == 0, replace with 0
        norm.counts[is.na(norm.counts)] <- 0
        # when target == 0, score == 0, keep that
    }

    norm <- counts[, 1:3]
    norm$score <- norm.counts
  
    # Also create a plot
    d.noscaling = density(na.omit(norm.noscaling), from=-2, to=6)
    d.scaling = density(na.omit(norm$score), from=-2, to=6)
    pdf(file.path(output_dir, paste0(basename, ".norm.density.pdf")), width=4.5, height=4)
    plot(d.noscaling, xlab="score", lwd=2, main="density and scaling", 
         ylim=c(0,max(d.noscaling$y, d.scaling$y)))
    abline(v=0, lty=2)
    abline(v=scaling_factor, lty=2, col="red")
    lines(d.scaling, lwd=2, col="red")
    dev.off()
    
    return(list(norm, scaling_factor))
    
    }

norm.OE <- function(counts, OE) {
    stop("Not implemented yet")
}

##############################################################################
## Supporting functions ######################################################
##############################################################################

logger <- function(s, quiet=opt$quiet, log_file=opt$log) {
    # Create logging messages. (Only printed when quiet == False)
    if (! quiet) {
        if (is.null(log_file)) {
            cat(paste0(s, "\n"), file = stderr())
        } else {
            cat(paste0(s, "\n"), file = log_file, append = TRUE)
        }
    }
}

read_counts <- function(file) {
    # Read counts table
    
    counts <- read.table(file, 
                         sep = "\t")
                         
    # I assume that columns 1-3 are location, and 4 is the count
    counts <- counts[, 1:4]
    names(counts) <- c("chr", "start", "stop", "counts")
    for (i in 2:4) {
        counts[, i] <- as.numeric(counts[, i])
    }
    									
    counts
}

##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.2"

option_list <- list(
    make_option(c("-c", "--counts"), type="character", default=NULL, 
                help="counts file (required)", metavar="character"),
    make_option(c("-m", "--method"), type="character", default=NULL, 
                help="normalization method (options: ...) (required)", metavar="character"),
    make_option(c("-o", "--output_dir"), type="character", default=NULL, 
                help="output directory (required)", metavar="character"),
                
    # Optional arguments                       
    make_option(c("-b", "--basename"), type="character", default=NULL,
                help="basename (default: basename count file)", metavar="character"),
    make_option(c("-d", "--dam_only"), type="character", default=NULL, 
                help="dam-only counts file (default: NULL)", metavar="character"),
    make_option(c("-p", "--pseudo"), type="numeric", default=1, 
                help="pseudocounts (default: 1)", metavar="character"),
    make_option(c("-s", "--scaling"), type="numeric", default=NULL, 
                help="scaling factor for dam-substr (default: 10% of data)", metavar="character"),
    make_option(c("-S", "--stdout"), type="logical", default=FALSE, action="store_true",
                help="report to stdout", metavar="character"),
	make_option(c("-Z", "--zip"), type="character", default=FALSE, action="store_true",
                help="Write gzipped output files", metavar="character"),
    make_option(c("-r", "--min_reads"), type="numeric", default=10, 
                help="Minimum number of reads in a bin (default: 10)", metavar="character"),
    make_option(c("-Q", "--quiet"), type="logical", default=FALSE, action="store_true", 
                help="quiet", metavar="character"),
    make_option(c("-L", "--log"), type="character", default=NULL, 
                help="Write debug to log file instead of stderr", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$counts) | is.null(opt$method) | is.null(opt$output_dir)) {
  	print_help(opt_parser)
  	stop("Not all required arguments given\n", call.=FALSE)
} else {
	main(opt$counts, opt$method, opt$output_dir, opt$basename, opt$dam_only,
         opt$pseudo, opt$scaling, opt$stdout, opt$zip, opt$min_reads)
}
