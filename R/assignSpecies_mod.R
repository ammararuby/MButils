#' @title Return long-form species assignments with DADA2's assignSpecies
#'
#' @description This function is a modification of DADA2's \code{assignSpecies},
#'   which uses exact matching against a reference fasta to identify the
#'   genus-species binomial classification of the input sequences. Where the
#'   original function run with \code{allowMultiple = TRUE} returns a
#'   concatenated string of all exactly matched species, results are returned
#'   here as distinct rows, one per match.
#'
#' @usage assignSpecies_mod(seqs, refFasta, tryRC = FALSE, n = 2000)
#'
#' @param seqs A character vector of the sequences to be assigned
#' @param refFasta The path to the reference fasta file, or an R connection. Can
#'   be compressed. This reference fasta file should be formatted so that the ID
#'   lines correspond to the genus-species of the associated sequence:
#'
#'   >SeqID genus species ACGAATGTGAAGTAA......
#' @param tryRC (Optional). Default FALSE. If TRUE, the reverse-complement of
#'   each sequences will also be tested for exact matching to the reference
#'   sequences.
#' @param n (Optional). Default 2000. The number of sequences to perform
#'   assignment on at one time. This controls the peak memory requirement so
#'   that large numbers of sequences are supported.
#'
#' @import dada2
#' @import Biostrings
#' @import ShortRead
#'
#' @return A two-column data frame matching ASVs to their exact match in the
#'   reference, with multiple matches indicated by the presence of more than one
#'   row per ASV.
#'
#' @export
#'
#'

assignSpecies_mod <- function(seqs, refFasta, tryRC=FALSE, n=2000) {
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  # Read in the reference fasta
  refsr <- ShortRead::readFasta(refFasta)
  ids <- as(ShortRead::id(refsr), "character")
  # Crude format check
  if(!length(unlist(strsplit(ids[[1]], "\\s"))) >= 3) {
    if(length(unlist(gregexpr(";", ids[[1]]))) >= 3) {
      stop("Incorrect reference file format for assignSpecies (this looks like a file formatted for assignTaxonomy).")
    } else {
      stop("Incorrect reference file format for assignSpecies.")
    }
  }
  ex <- '>\\d+\\s'
  species <- sapply(unlist(ids), sub, pattern = ex, replacement = '')
  names(species) <- NULL

  # Identify the exact hits
  hits <- vector("list", length(seqs))
  lens <- nchar(seqs)
  for(len in unique(lens)) { # Requires all same length sequences
    i.len <- which(lens==len) # Find indexes of ASVs equal to that length
    n.len <- length(i.len) # Number of ASVs equal to that length

    # This throttles the number of sequences that can be analyzed at any one time
    j.lo <- 1
    j.hi <- min(n,n.len)

    while(j.lo <= n.len) {
      i.loop <- i.len[j.lo:n.len] # Creates a list of indexes to loop over
      seqdict <- Biostrings::PDict(seqs[i.loop])
      vhit <- (Biostrings::vcountPDict(seqdict, ShortRead::sread(refsr))>0)
      if(tryRC) vhit <- vhit | (Biostrings::vcountPDict(seqdict, Biostrings::reverseComplement(ShortRead::sread(refsr)))>0)
      hits[i.loop] <- lapply(seq(nrow(vhit)), function(x) vhit[x,]) # A logical vector of which reference sequences were "hit"
      j.lo <- j.lo + n; j.hi <- min(j.hi+n, n.len)
      rm(seqdict)
      gc()
    }
  }
  # Get genus species return strings
  rval <- cbind(unlist(sapply(hits, mapHits_mod, refs=species)))
  colnames(rval) <- c("Species")
  rownames(rval) <- seqs
  # Reformat rval dataframe
  rval <-
    rval %>%
    as.data.frame() %>%
    rownames_to_column(var = 'asv') %>%
    separate_rows(Species, sep = '\\/') # "Lengthen" data-- split /-separated column entries out to multiple rows

  # TODO: Last common ancestor??
  gc()
  rval
}

# Helper function for assignSpecies_mod
mapHits_mod <- function(x, refs, sep="/") {
  hits <- refs[x]
  if(length(unique(hits))<=Inf) {
    rval <- do.call(paste, c(as.list(sort(unique(hits))), sep=sep))
  } else { rval <- NA_character_ }
  if(length(rval)==0) rval <- NA_character_
  rval
}
