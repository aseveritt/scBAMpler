#Author: Amanda Everitt
#in R

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-b", "--bam_file"), type="character",  default=NULL, action="store",
                     help="Path to bamfile", 
                     dest="bam_file")
parser <- add_option(parser, c("-o", "--outdir"), type="character",  default=NULL, action="store",
                     help="Path to output diretory", 
                     dest="out_dir")
parser <- add_option(parser, c("-l", "--peak_length"), type="integer",  default=500, action="store",
                     help="Length to standardize peaks to", 
                     dest="peaklen")
parser <- add_option(parser, c("--cores"), type="integer",  default=1, action="store",
                     help="number of cores for mclapply", 
                     dest="ncores")
parser <- add_option(parser, c("--txdb"), type = "character", default = "TxDb.Hsapiens.UCSC.hg38.knownGene", action="store",
                     help = "TxDb package name to load [default %default]",
                     dest="txdb_package")
parser <- add_option(parser, c("--exclusion_file"), type = "character", default = NULL, action="store",
                     help = "If there are regions to exclude from peak set, provide the full path",
                     dest="exclusion_list")
parser <- add_option(parser, c("-s", "--summit_file"), type="character",  default=NULL, action="store",
                     help="Path to summit file if you're not running macs", 
                     dest="summit_file")
parser <- add_option(parser, c("--union_files"), type="character", default=NULL, action="store",
                     help="Comma-separated list of standardized BED files to merge into a union peak set",
                     dest="union_files")
parser <- add_option(parser, c("--union_outfile"), type="character", default="union_peaks.bed", action="store",
                     help="Output filename for union peak set [default %default]",
                     dest="union_outfile")

opt <- parse_args(parser)

if (is.null(opt$out_dir)) { stop("ERROR: Missing --outdir argument.") }
if (is.null(opt$union_files) && is.null(opt$summit_file) && is.null(opt$bam_file)) {
    stop("ERROR: Provide --bam_file (to run macs3), or --summit_file (to standardize an existing macs3 summit file).")
}

################################################################################################

#resolve the helper path relative to this script, so call_peaks.R can be run from any directory
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
source(file.path(script_dir, "call_peak_functions.R"))


if (!is.null(opt$union_files)) {
    infiles <- strsplit(opt$union_files, ",")[[1]]
    make_union(infiles  = infiles,
               outdir   = opt$out_dir,
               outfile  = opt$union_outfile,
               ncores   = opt$ncores)
} else {
    if (is.null(opt$summit_file)) { 
        prefix = gsub(".bam", "", basename(opt$bam_file)) 
        call_macs(opt$bam_file, opt$out_dir, prefix) 
        opt$summit_file = file.path(opt$out_dir, paste0(prefix, "_summits.bed"))
    } 
    standardize_summits(summit_file    = opt$summit_file, 
                        out_dir        = opt$out_dir, 
                        exclusion_list = opt$exclusion_list, 
                        peaklen        = opt$peaklen, 
                        txdb           = opt$txdb_package, 
                        ncores         = opt$ncores)
}


