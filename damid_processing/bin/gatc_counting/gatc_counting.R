#!/usr/bin/Rscript

# GATC counting (version 1.1)
## Tom van Schaik, 160828

# - Functionality:
# Count reads to GATC fragments (or any kind of gff-file in theory). Using
# GRanges objects, overlaps are determined (specifically at the ends / starts
# of fragments), and a bed-like file is returned. Several options exist to 
# alter the behavious of the counting.

# - Input
# (required)
# * bam file
# * GATC gff file - with the fragments
# * output directory
#
# (optional)
# * basename - basename of the output files (default: basename from file)
# * maxgap - maximum gap allowed between the start/end and the read (default: 0)
# * shift - shift reads (default 0)
# * onehit - only allow each GATC once (for single cell?)
# * bin - bin reads (not working yet, how to bin?)
# * genome_size - file with genome sizes, as returned from fetchChromSizes of UCSC
#   (required when binning)
# * mapqual - filter based on a mapping quality
# * multimap - filter out multimapping reads
# * dup - filter out duplicated reads
# * quiet - be quiet please (default: False)
# * ...

# - Output (within output_dir)
# * [BASENAME].counts.txt - bed-like output format
# * [BASENAME].counts.statistics.txt - accompanying statistics (i.e. reads counted)

# - Method
# Create GRanges from the fragment GATC. Read alignments, and filter those as given.
# Determine overlap on the fragment ends, taking strand information into account.
# (Bin if required.) Return bed-like format.

# - Version
# 1.0 - First version
# 1.1 - Binning now doesn't require a separate run

# - To do
# * Chunking of reads, instead of all in one go (saves memory)
# * Paired end reads (readGAlignmentPairs), but then: what is unique?!
#   And when do they overlap fragment ends? Both, or just one?
# * More statistics to be returned - which ones?
# * Create some figures - which figures?!
# * Allow for non-end-overlapping counting for more general use

suppressPackageStartupMessages(library("optparse"))
options(warn=1)

##############################################################################
## Main function #############################################################
##############################################################################

main <- function(bam, gff, output_dir, basename=NULL, maxgap=0, shift=0,
                 onehit=FALSE, bin = NULL, genome_size = NULL, mapqual=NULL, 
                 multimap=FALSE, dup=FALSE, stdout=FALSE, zip=FALSE, folder=FALSE) {
    # Main function with the flow of the module. See comments for a short 
    # description of each step. 
    
    # Setting-up
    logger(paste0("----- gatc-counting version: ", script.version, " -----"))
    logger("Setting-up")
    
    # First, separate bam files
    bam <- strsplit(bam, ",")[[1]]
    
    if (length(bam) > 1 & is.null(basename)) 
        stop(sprintf("basename is required if input files > 1"))
    
    for (b in bam) {
        if (!file.exists(b))
            stop(sprintf("sample file %s does not exist, aborting.", b))
    }
    
    if (!file.exists(gff))
        stop(sprintf("sample file %s does not exist, aborting.", gff))
    
    if (is.null(basename)) {
        basename <- sub(".bam", "", basename(bam))
    }
    
    dir.create(output_dir, showWarnings = FALSE)
	
    if (! stdout) {
        if (folder) {
            out <- file.path(output_dir, "bin-gatc", paste0(basename, "-gatc.counts.txt"))
        } else {
            out <- file.path(output_dir, paste0(basename, "-gatc.counts.txt"))
        }
        if (zip) {
            out <- paste0(out, ".gz")
        }
    } else {
        out <- ""
    }
    out_stat <- file.path(output_dir, paste0(basename, ".counts.statistics.txt"))
    
    if (! bin == FALSE) {
        if (is.null(genome_size) | ! file.exists(genome_size)) {
            stop("Please provide a valid genome_size file.")
        }
    }
    
    logger(paste0("\tbam file: ", bam))
    logger(paste0("\tgff file: ", gff))
    logger(paste0("\tOutput file: ", out, " (Note: [nothing] = stdout)"))
    
    # Load dependencies (only now, because it's quite slow)
    logger("Loading dependencies")
    
    suppressWarnings(suppressPackageStartupMessages(library("GenomicRanges")))
    suppressWarnings(suppressPackageStartupMessages(library("GenomicAlignments")))
    suppressWarnings(suppressPackageStartupMessages(library("rtracklayer")))
    
    # Read data
    logger("\n----------")
    logger(paste0("Reading gff ", gff, "..."))
    gff.ranges <- read_gff(gff)
    logger(paste0(length(gff.ranges), " features read"))
    
    reads <- bam_reader(bam, mapqual, multimap, dup)
    logger(paste0(length(reads), " alignments read"))
    
    # Shift reads
    if (shift != 0) {
        reads <- shift(reads, shift*ifelse(strand(reads) == '+', 1, -1))
    }
    
    # Determine overlaps
    logger("\n----------")
    logger("Overlapping features and alignments")
    fwd_index <- strand(reads) == "+"
    fwd <- countOverlaps(gff.ranges, reads[fwd_index],
                         type = "start", maxgap = maxgap)
    rev <- countOverlaps(gff.ranges, reads[!fwd_index],
                         type = "end", maxgap = maxgap)
    total <- fwd + rev
    
    # One hit-GATCs
    if (onehit) {
        logger("Reducing hits to 1 / fragment")
        total[total > 1] <- 1
    }
    
    mcols(gff.ranges)[, "count"] <- total
    
    # For now, keep only the "normal" chromosomes
    normal_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
    
    if (! is.null(normal_chromosomes)) {
        
        logger(paste0("Selecting for chromosomes: ",
                      normal_chromosomes))
        
        # Do notify how many counts are NOT in the normal chromosomes
        n_total <- sum(mcols(gff.ranges)[, "count"])
        n_removed <- sum(mcols(gff.ranges)[! seqnames(gff.ranges) %in% normal_chromosomes, "count"])
        
        logger(paste("Of",
                     n_total,
                     "reads total,",
                     n_removed,
                     "were removed as they belonged to the wrong chromosome."))
        
        # Mouse data
        seqlevels(gff.ranges, pruning.mode="coarse") <- c(paste0("chr", 1:19), "chrX", "chrY")
        
    }
    
    # Count reads in bins and write the output immediately
    if (! bin == FALSE) {

    # Split the different bin widths requested 
      bin <- as.integer(strsplit(bin, ",")[[1]])
      if (is.unsorted(bin)) {
        stop("Please provide bins in increasing size")
      }
      
      # If bins are not present, use the GATC fragments. Otherwise the last used bins
      bins <- NULL
      for (w in bin) {
        
        logger("\n----------")
        logger(paste0("Counting reads in ", w, "kb bins"))
        if (folder) {
            file_name <- file.path(output_dir, 
                                   paste0("bin-", w, "kb"),
                                   paste0(basename, "-", w, "kb.counts.txt"))
        } else {
            file_name <- file.path(output_dir, 
                                   paste0(basename, "-", w, "kb.counts.txt"))
        }
        if (zip) {
          file_name <- paste0(file_name, ".gz")
        }
        w <- w*1000
        
        # Create and count the bins
        if (is.null(bins)) {
          bins <- get_bins(gff.ranges, genome_size, w)
          bins <- count_bins(bins, gff.ranges)
          write_counts(bins, file_name)
        } else {
          bins.new <- get_bins(gff.ranges, genome_size, w)
          bins.new <- count_bins(bins.new, bins)
          bins <- bins.new
          write_counts(bins, file_name)
        }
      }
    } 
    
    # Write output
    logger("\n----------")
    logger("Writing GATC output")   
    
    write_counts(gff.ranges, out)
    
    # Create statistics
    logger("\n----------")
    logger(paste0("Reads used: ", length(reads)))
    logger(paste0("Reads counted: ", sum(total)))
    
    stat <- data.frame(n = c("reads", "counted"),
                       c = c(length(reads), sum(total)))
    write.table(stat, out_stat, quote = FALSE, sep = "\t", col.names = FALSE,
                row.names = FALSE)
    
    # Add versions
    write(paste0("## Versions ##"), out_stat, append = TRUE)
    write(paste0("# R: ", version["version.string"][[1]]), out_stat, append = TRUE)
    write(paste0("# gatc-counting: ", script.version), out_stat, append = TRUE)
    
    # Report session info
    logger("\n----------")
    logger("Session info:")
    logger(capture.output(sessionInfo()))
    logger("")
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

read_gff <- function(gff) {
    # Import gff.file into GRanges
    gff <- import.gff(gff)
    
    #gff <- as(gff, "GRanges")
    isCircular(gff)[TRUE] <- FALSE
    gff
}

bam_reader <- function(bam, mapqual, multimap, dup) {
    # Read alignments from a bam-file, applying given filters on it.
    param <- get_param(mapqual, dup)
    
    reads <- GAlignments()
    for (b in bam) {
        logger(paste0("Reading bam file ", b, "..."))
        reads <- c(reads, readGAlignments(b, param = param))
    }
    
    if (multimap) {
        # Bowtie-specific! - the XS-flag is only present when multiple 
        # alignments are present
        reads <- reads[is.na(mcols(reads)$XS)]
    }
    
    # Convert to GRanges
    reads <- as(reads, "GRanges")
    reads
}

get_param <- function(mapqual, dup) {
    # Get the parameter options for reading alignments
    flags <- scanBamFlag(isUnmappedQuery = FALSE,
                         isDuplicate = dup)
    param <- ScanBamParam(flag = flags, 
                          tag = "XS",
                          mapqFilter = mapqual)
    param
}

get_bins <- function(gff, genome_size, bin) {
    # Create bins of width "bin", using a genome size table (file)
    
    # First, read the genome size and select chromosomes
    seqnames <- levels(seqnames(gff))
    g <- read.table(genome_size, sep = "\t")
    names(g) <- c("chr", "size")
    g <- g[g$chr %in% seqnames, ]
    
    # For each chromosome, create bins
    start <- unlist(sapply(g$chr, function(x) seq(from=1, by=bin, to=g[g$chr == x, 2])))
    end <- unlist(sapply(g$chr, function(x) {
        s <- seq(from=bin, by=bin, to=ceiling(g[g$chr == x, 2]/bin)*bin)
        s[length(s)] <- g[g$chr == x, 2]
        s
    }))
    
    # Assert that everything is okay
    if (length(start) != length(end)) {
        stop("Oops, something went wrong with the binning")
    }
    
    # Create GRanges
    bins <- GRanges(seqnames = rep(g$chr, ceiling(g$size / bin)),
                    ranges = IRanges(start, end))
    bins <- sortSeqlevels(bins)
    bins <- sort(bins)
    bins
}

count_bins <- function(bins, gff.ranges) {
    # Collapse fragment counts into bins. For this, use the middle of the bin
    
    gff.middle <- start(gff.ranges) + floor((end(gff.ranges)-start(gff.ranges))/2)
    gff.middle <- GRanges(seqnames(gff.ranges),
                          IRanges(gff.middle, gff.middle))
                          
    ovl <- findOverlaps(gff.middle, bins, type='within')
    # stopifnot(all.equal(seq_along(gff.middle), queryHits(ovl)))
    if (all.equal(seq_along(gff.middle), queryHits(ovl)) != TRUE) {
        logger("Note: not all fragments counted into bins!")
    }
    
    # Note: put the unique subjectHits in order (chromosome order)    
    mcols(bins)$count <- 0 
    mcols(bins)$count[sort(unique(subjectHits(ovl)))] <- tapply(
        queryHits(ovl),
        subjectHits(ovl),
        function(x) sum(mcols(gff.ranges)$count[x]))
    
    bins
}

write_counts <- function(gff.ranges, output) {
    d <- as(gff.ranges, "data.frame")
    d <- d[, c("seqnames", "start", "end", "count")]
    
    # Output format is bed-like, which is 0-based. R is 1-based. 
    # Correct for this.
    d$start <- d$start - 1
        
    options(scipen = 999) 
    
    if (endsWith(output, ".gz")) {
        output <- gzfile(output, "w")
        write.table(d, file = output, quote = FALSE, sep = "\t", col.names = FALSE,
                    row.names = F)
        close(output)
    } else { 
        write.table(d, file = output, quote = FALSE, sep = "\t", col.names = FALSE,
                    row.names = F)
    }
}

##############################################################################
## User input ################################################################
##############################################################################

script.version <- "1.1"

# Get input
option_list <- list(
    make_option(c("-b", "--bam"), type="character", default=NULL, 
                help="bam file(s), comma-separated (required)", metavar="character"),
    make_option(c("-f", "--gff"), type="character", default=NULL, 
                help="GATC gff file (required)", metavar="character"),
    make_option(c("-o", "--output_dir"), type="character", default=NULL,
                help="output directory (required)", metavar="character"),
                
    # Optional arguments            
    make_option(c("-n", "--basename"), type="character", default=NULL, 
                help="basename (default: basename bam file)", metavar="character"),
    make_option(c("-a", "--maxgap"), type="integer", default=0, 
                help="maxgap for overlap (see IRanges, findOverlaps) (default: 0)", 
                metavar="character"),
    make_option(c("-s", "--shift"), type="integer", default=0, 
                help="shift for overlap (see IRanges, findOverlaps) (default: 0)", 
                metavar="character"),
    make_option(c("-l", "--onehit"), type="logical", default=FALSE, action="store_true",
                help="no duplicated GATC fragments (thus 0 or 1) (default: False)", 
                metavar="character"),
    make_option(c("-B", "--bin"), type="character", default=FALSE, 
                help="binning width(s) in kb, comma separated (default: GATC fragments only)", 
				metavar="character"),
    make_option(c("-G", "--genome_size"), type="character", default="", 
                help="genome sizes, required for binning (see 'hg19.chrom.sizes' in \
                script directory) (default: NULL)", metavar="character"),
    make_option(c("-q", "--mapqual"), type="integer", default=0, 
                help="mapping quality threshold (default: 0)", metavar="character"),
    make_option(c("-m", "--multimap"), type="logical", default=FALSE, action="store_true", 
                help="no multimappers (bowtie XS-tag present) (default: False)", 
                metavar="character"),
    make_option(c("-d", "--dup"), type="logical", default=NA, action="store_false",
                help="no duplicated reads (default: NA)", metavar="character"),
    make_option(c("-S", "--stdout"), type="logical", default=FALSE, action="store_true",
                help="output to stdout (still creates the statistics file) (default: False)", 
                metavar="character"),
    make_option(c("-Z", "--zip"), type="character", default=FALSE, action="store_true",
                help="Write gzipped output files", metavar="character"),
    make_option(c("-F", "--folder"), type="character", default=FALSE, action="store_true",
                help="Create folders for each (bin-) size", metavar="character"),
    make_option(c("-Q", "--quiet"), type="logical", default=FALSE, action="store_true",
                help="Quiet", metavar="character"),
    make_option(c("-L", "--log"), type="character", default=NULL, 
                help="Write debug to log file instead of stderr", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$bam) | is.null(opt$gff) | is.null(opt$output_dir)) {
  	print_help(opt_parser)
  	stop("Not all required arguments present.\n", call.=FALSE)
} else {
	main(opt$bam, opt$gff, opt$output_dir, opt$basename, opt$maxgap, opt$shift, 
         opt$onehit, opt$bin, opt$genome_size, opt$mapqual, opt$multimap, opt$dup,
         opt$stdout, opt$zip, opt$folder)
}
