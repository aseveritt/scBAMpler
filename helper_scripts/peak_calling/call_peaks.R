#Author: Amanda Everitt
#in R

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-b", "--bam"), type="character",  default=NULL, action="store",
                     help="Path to bamfile", 
                     dest="bam_file")
parser <- add_option(parser, c("-o", "--outdir"), type="character",  default=NULL, action="store",
                     help="Path to output diretory", 
                     dest="out_dir")
parser <- add_option(parser, c("-len", "--peak_len"), type="float",  default=500, action="store",
                     help="Length to standardize peaks to", 
                     dest="peaklen")
parser <- add_option(parser, c("--cores"), type="integer",  default=1, action="store",
                     help="number of cores for mclapply", 
                     dest="ncores")
parser <- add_option(parser, c("--txdb"), type = "character", default = "TxDb.Hsapiens.UCSC.hg38.knownGene", action="store",
                     help = "TxDb package name to load [default %default]",
                     dest="txdb_package")
parser <- add_option(parser, c("--exclusion_list"), type = "character", default = NULL, action="store",
                     help = "If there are regions to exclude from peak set, provide the full path",
                     dest="exclusion_list")
parser <- add_option(parser, c("-s", "--summitfile"), type="character",  default=NULL, action="store",
                     help="Path to summit file if you're not running macs", 
                     dest="summit_file")

opt <- parse_args(parser)

if (is.null(opt$out_dir)) { stop("ERROR: Missing --out_dir argument.") }

if (!opt$skipmacs){ if (is.null(opt$bam_file)) { stop("ERROR: Missing --bam argument.") } }
if (opt$skipmacs){ if (is.null(opt$summit_file)){ stop("ERROR: Must add summit_file argument if macs isnt being run.") }}

################################################################################################


source("helper_scripts/peak_calling/call_peak_functions.R")
if (!is.null(summit_file)) { 
    prefix = gsub(".bam", "", basename(opt$bam_file)) 
    
    call_macs(opt$bam_file, opt$out_dir, prefix) 
    opt$summit_file = paste0(opt$out_dir, prefix, "_summits.bed")
} 

standardize_summits(summit_file    = opt$summit_file, 
                    out_dir        = opt$out_dir, 
                    exclusion_list = opt$exclusion_list, 
                    peaklen        = opt$peaklen, 
                    txdb           = opt$txdb, 
                    ncores         = opt$ncores)
    
#########################################################################

