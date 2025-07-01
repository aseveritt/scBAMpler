#Author: Amanda Everitt
#in R

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures) 
  library(dplyr)
  library(rtracklayer)
  library(stringr)
  library(data.table)
  library(parallel)
  library(txdb_pkg, character.only = TRUE)
})

`%ni%` <- Negate(`%in%`)


################################################################################################

call_macs <- function(bam_file, out_dir, prefix){
  
  template = c('macs3 callpeak',
              '-g 2.9e+09',
              '--name %s',
              '--treatment %s',
              '--outdir %s',
              '--call-summits',
              '--keep-dup all',
              '--nomodel',
              '--nolambda',
              '--shift -100',
              '--extsize 200', 
              '-q 0.05')
  macs_call <- sprintf(paste(template, collapse = " "),
                       prefix, bam_file, out_dir)
  
  output <- system(macs_call, intern = TRUE)
  
}

################################################################################################
#### Internal functions

cleanHitsObject <- function(gr, hits_df){
  
  #remove self and mirrored hits. 
  hits_df_cleaned = hits_df %>% 
    filter(queryHits != subjectHits) %>%
    filter(queryHits < subjectHits) 
  
  #build a table that shows how many overlaps each hit has so we can filter later. 
  freq_table = data.frame(table(c(hits_df_cleaned$queryHits, hits_df_cleaned$subjectHits)))
  freq_table$Var1 <- as.numeric(as.character(freq_table$Var1))
  
  #add useful info columns
  hits_df_cleaned = hits_df_cleaned %>%
    left_join(freq_table, by = c("subjectHits" = "Var1")) %>% 
    left_join(freq_table, by = c("queryHits" = "Var1"), suffix = c(".subject", ".query")) %>%
    mutate(
      queryScore = mcols(gr)$score[queryHits],
      subjectScore = mcols(gr)$score[subjectHits],
      maxHit = ifelse(queryScore >= subjectScore, queryHits, subjectHits),
      minHit = ifelse(queryScore >= subjectScore, subjectHits, queryHits),
      maxScore = ifelse(queryScore >= subjectScore, queryScore, subjectScore),
      maxFreq = pmax(Freq.subject, Freq.query)
    ) 
  
  return(hits_df_cleaned)
}



getuniquePeaks <- function(gr, hits_df_cleaned){
  allidx = 1:length(gr)
  overlapped_somewhere <- unique(c(hits_df_cleaned$subjectHits, hits_df_cleaned$queryHits))
  unique_peaks = allidx[allidx %ni% overlapped_somewhere]
  return(unique_peaks)
}



iterativeOverlap <- function(hits_df_cleaned, quiet=FALSE){
  
  # Filter to only relevant rows and rank
  dt <- as.data.table(hits_df_cleaned)
  dt <- dt[maxFreq > 1][, rank := frank(-maxScore, ties.method = "random")]
  setorder(dt, rank)
  if (! quiet) { message(sprintf("%s peaks overlap > 2 regions", n_distinct(c(dt$queryHits, dt$subjectHits))))}
  
  max_greaterthan2 <- integer(0)
  seen_hits <- integer(0)
  all_maxes = dt$maxHit 

  for (i in seq_along(all_maxes)) {
    current_hit <- all_maxes[i]
    
    if (current_hit %in% seen_hits) next 
    #if the entry was already in a comparison we don't need it.   
    
    max_greaterthan2 = c(max_greaterthan2, current_hit) 
    #if we haven't seen it, and we're going in rank order, than by default its the max of this comparision and should be retained. 
    
    toremove <- unique(dt[queryHits == current_hit | subjectHits == current_hit, c(queryHits, subjectHits)]) 
    #get ids of all associated peaks including the max 
    
    seen_hits <- c(seen_hits, toremove) 
    #remove them all
    
  } 
  
  return(max_greaterthan2)
}


wrapper_for_per_chromosome <- function(gr, maxgap){
  
  hits_df <- as.data.table(findOverlaps(gr, gr, maxgap = maxgap, ignore.strand = TRUE, select="all"))
  hits_df_cleaned = cleanHitsObject(gr, hits_df)
  
  #First, lets find all peaks that didn't overlap anywhere
  unique_peaks <- getuniquePeaks(gr, hits_df_cleaned)
  
  #Next, lets quickly remove all the overlaps that are just with one other to save time. 
  max_groupsof2 = hits_df_cleaned %>% dplyr::filter(maxFreq == 1) %>% pull(maxHit)

  #Finally, lets do the iterative overlap for all remaining entries
  max_greaterthan2 <- iterativeOverlap(hits_df_cleaned, quiet=TRUE)
  
  #note: remember that the hits_df values are relative to this chunk. 
  #so here we want to return the idx in the full ranges object. 
  return(list(gr[unique_peaks]$idx, gr[max_groupsof2]$idx, gr[max_greaterthan2]$idx))
}



checkfinalhits <- function(gr){
  hits_df <- as.data.frame(findOverlaps(gr, gr, ignore.strand = TRUE, select="all"))
  non_self_hits <- hits_df[hits_df$queryHits != hits_df$subjectHits, ]
  if (nrow(non_self_hits) > 0){print(paste("ERROR, there were overlapping regions at the final step -- you should evaluate how and why that happened before proceeding."))}
}



################################################################################################



standardize_summits <- function(summit_file, out_dir, exclusion_list, peaklen, txdb, ncores){

    if (!requireNamespace(txdb, quietly = TRUE)) {
        stop(sprintf("Package '%s' not installed.", txdb))
    }
    
  ########################
  # Load Summits File
  summits <- read.table(summit_file, header = FALSE, stringsAsFactors = FALSE)
  prefix = gsub("summits.bed", "", summit_file)
  
  gr <- GRanges(
    seqnames = summits$V1,
    ranges = IRanges(start = summits$V2, end = summits$V3),
    strand = "*",
    score = summits$V5,
    name = gsub(prefix, "", summits$V4)
  )
  gr$idx = 1:length(gr)
    
  ########################
  data_chunks <- as.list(split(gr, seqnames(gr)))
  results <- mclapply(seq_along(data_chunks), 
                    function(i){wrapper_for_per_chromosome(data_chunks[[i]], peaklen)},
                    mc.cores = ncores)
  
  unique_peaks = unlist(lapply(results, `[[`, 1))
  max_groupsof2 = unlist(lapply(results, `[[`, 2))
  max_greaterthan2 = unlist(lapply(results, `[[`, 3))
  union_peaks = c(unique_peaks, max_groupsof2, max_greaterthan2)
  
  message(sprintf("%s total summits to begin with", length(gr)))
  message(sprintf("%s peaks overlaped no other peak", length(unique_peaks)))
  message(sprintf("%s peaks overlaped one other peak -- removing %s peaks", length(max_groupsof2), length(max_groupsof2)))
  message(sprintf("the remaining peaks were summarized into %s peaks", length(max_greaterthan2)))
  message(sprintf("Ultimately, there are %s peaks left", length(union_peaks)))
  
  gr_filt = gr[union_peaks]
  checkfinalhits(gr_filt)
  
  ########################
  # Extend Summits
  gr_filt_resized <- resize(gr_filt, width = peaklen, fix = "center")
  
  ########################
  # Make sure chromosome boundaries aren't invalidated
  chrom_lengths <- seqlengths(txdb)
  out_of_bounds <- gr_filt_resized[end(gr_filt_resized) > chrom_lengths[as.character(seqnames(gr_filt_resized))]]
  if (length(out_of_bounds) > 0){print(paste("ERROR, regions out of bounds", out_of_bounds))}
  
  
  ########################
  # Remove regions overlapping blacklist regions
  if (!is.null(exclusion_list)) { 
      bl <- read.table(, sep="\t", header = FALSE, stringsAsFactors = FALSE)
      bl_gr <- GRanges(
        seqnames = bl$V1,
        ranges = IRanges(start = bl$V2, end = bl$V3),
        strand = "*"
      )
  
      toremove <- data.frame(findOverlaps(gr_filt_resized, bl_gr, ignore.strand = TRUE, minoverlap = 150))
      gr_filt_resized = gr_filt_resized[-unique(toremove$queryHits)] #remove bl regions
  }
  
  ########################
  # Output them
  if (unique(width(gr_filt_resized)) != peaklen){print(paste("ERROR, regions not equal length", unique(width(gr_filt_resized))))}
  gr_filt_resized = sort(gr_filt_resized) #YOU HAVE TO SORT ITTT
  
  cl = gsub("_summits.bed", "", basename(summit_file))
  export(gr_filt_resized, paste0(out_dir, cl, "_filt_", peaklen, "bp.exclusion.bed"), format = "bed")
  
}





################################################################################################


make_union <- function(infiles, outdir, outfile, ncores=opt$ncores){
  
  ########################
  #Read in and concatenate all peaks
  allpeaks = data.frame()
  
  for (peak_file in infiles){
    f1 <- read.table(paste0(outdir, peak_file), header = FALSE, stringsAsFactors = FALSE)
    f1$norm_score = (f1$V5 / sum(f1$V5))*1e6 #normalize scores: I'm sure there are better ways of doing this, but this will work for our purposes
    f1$cell_line = str_split(peak_file, "_")[[1]][1]
    allpeaks <- rbind(allpeaks, f1)
  }
  
  allpeaks <- allpeaks %>% distinct(V1, V2, V3, .keep_all = TRUE) # filter any identical ones out
  
  gr <- GRanges(
    seqnames = allpeaks$V1,
    ranges = IRanges(start = allpeaks$V2, end = allpeaks$V3),
    strand = "*",
    score = allpeaks$norm_score,
    origin_score = allpeaks$V5,
    origin = allpeaks$cell_line
  )
  gr$idx = 1:length(gr)
    
  ########################
  data_chunks <- as.list(split(gr, seqnames(gr)))
  results <- mclapply(seq_along(data_chunks), 
                      function(i){wrapper_for_per_chromosome(data_chunks[[i]], -1)},
                      mc.cores = ncores)
    
  unique_peaks = unlist(lapply(results, `[[`, 1))
  max_groupsof2 = unlist(lapply(results, `[[`, 2))
  max_greaterthan2 = unlist(lapply(results, `[[`, 3))
  union_peaks = c(unique_peaks, max_groupsof2, max_greaterthan2)
  
  message(sprintf("%s total summits to begin with", length(gr)))
  message(sprintf("%s peaks overlaped no other peak", length(unique_peaks)))
  message(sprintf("%s peaks overlaped one other peak -- removing %s peaks", length(max_groupsof2), length(max_groupsof2)))
  message(sprintf("the remaining peaks were summarized into %s peaks", length(max_greaterthan2)))
  message(sprintf("Ultimately, there are %s peaks left", length(union_peaks)))
  
  final_gr = gr[union_peaks]
  checkfinalhits(final_gr)
  final_gr = sort(final_gr) 
  export(final_gr, paste0(outdir, outfile), format = "bed")
  
}





########################################################################

quickunittest <- function(){
  
  gr <- GRanges(seqnames = Rle(rep("chr1", 8)),
                ranges = IRanges(start = c(100, 2000, 2500, 4000, 6000, 6200, 6600, 7000),
                                 end = c(100, 2000, 2500, 4000, 6000, 6200, 6600, 7000)+1),
                score = c(10, 15, 20, 5, 30, 26, 10, 25))
  
  
  ########################
  #Get all overlapping regions
  hits_df <- as.data.frame(findOverlaps(gr, gr, maxgap = 500, ignore.strand = TRUE, select="all"))
  hits_df_cleaned = cleanHitsObject(gr, hits_df)
  
  #First, lets find all peaks that didn't overlap anywhere
  unique_peaks <- getuniquePeaks(gr, hits_df_cleaned)
  if (all(unique_peaks != c(1,4))){message("ERROR IN UNIQUE PEAKS")}
  
  
  #Next, lets quickly remove all the overlaps that are just with one other to save time. 
  max_groupsof2 = hits_df_cleaned %>% dplyr::filter(maxFreq == 1) %>% pull(maxHit)
  if (all(max_groupsof2 != c(3))){message("ERROR IN GROUPS OF 2")}
  
  ########################
  #Finally, lets do the iterative overlap for all remaining entries
  max_greaterthan2 <- iterativeOverlap(hits_df_cleaned)
  if (all(max_greaterthan2 != c(5,8))){message("ERROR IN ITERATIVE OVERLAP")}
  
  ########################
  union_peaks = c(unique_peaks, max_groupsof2, max_greaterthan2)
  if (all(sort(union_peaks) !=  c(1,3,4,5,8))){message("ERROR IN UNION PEAKS")}
  
  message("unit test done\n\n")
  
}

#quickunittest()