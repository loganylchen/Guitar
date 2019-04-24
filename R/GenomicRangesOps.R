GRangesListmapToTranscripts <- function(x, mapFilterTranscript,transcripts, ignore.strand=FALSE)
{ 
  if(class(x) == "CompressedGRangesList")
  {  
    names(x) <- 1:length(x)
    xWidthes <- sum(width(x))
    names(xWidthes) <- names(x)
     x_unlisted <- unlist(x)
  }else{
    names(x) <- 1:length(x)
    xWidthes <- width(x)
    names(xWidthes) <- names(x)
     x_unlisted <- x
  }
  
  tx_coord <- mapToTranscripts(x_unlisted, transcripts,ignore.strand=FALSE)
  
  xHit_txHit_joint <- paste(names(tx_coord), tx_coord$transcriptsHits, sep='-')
  tx_coord_grouped <- split(tx_coord, xHit_txHit_joint)
  mapping_reduced <- reduce(tx_coord_grouped)
  
  # some sites map to tx has mutiple regions after reduce because of isoform of tx
  mapping_reduced_width <- width(mapping_reduced)
  mapping_region_nums <- lapply(mapping_reduced_width, function(x) length(x))
  index_of_continous <- which(mapping_region_nums == 1)
  mapping_filter <- mapping_reduced[index_of_continous]
  tx_coord_filtered <- unlist(mapping_filter)
  
  xHit_txHit <- strsplit(names(tx_coord_filtered), '-')
  xHits <- as.numeric(lapply(xHit_txHit, `[`, 1))
  txHits <- as.numeric(lapply(xHit_txHit, `[`, 2))
  mcols(tx_coord_filtered) <- data.frame(xHits = xHits, txHits= txHits)
  
  # remove hits whoes length smaller than sites because of isoform
  if(mapFilterTranscript) {
    tx_coord_filtered_width <- width(tx_coord_filtered)
    tx_coord_filtered <- tx_coord_filtered[tx_coord_filtered_width == xWidthes[tx_coord_filtered$xHits]]
  }
  return(tx_coord_filtered)
}
