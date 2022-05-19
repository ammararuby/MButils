#' @title Extract the finest taxonomic level of ASVs from a phyloseq taxonomy
#'
#' @description Identifies the lowest taxonomic name available for each ASV in a
#'   phyloseq taxonomy table and aggregates them in a new column.
#'
#' @param taxtab The taxonomy table of the phyloseq object to be renamed
#'
#' @import phyloseq
#'
#' @return An updated taxonomy table, with an added column, 'name', that
#'   contains the name of the lowest phylogenetic level to which that ASV is
#'   identified.
#' @export
#'
#'

lowest_level <- function(taxtab){
     # Update taxa names from ASV sequence to identified taxon at the most
     # precise phylogenetic level possible

     # This gets the right-most, non-NA value
     lowest.index <- max.col(!is.na(taxtab), 'last')
     taxtab$name <- taxtab[cbind(seq_along(lowest.index),
                                         lowest.index)]

     taxtab
}
