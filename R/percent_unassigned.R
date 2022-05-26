#' @title Find the percentage of reads without a taxon assignment
#'
#' @description Identifies sequence variants with no taxon assignment (*i.e.*,
#'   those for which every column of the taxonomy table is NA).  Then,
#'   determines what percentage of overall reads in the dataset correspond to
#'   unidentified reads.  This is helpful for determining if the current
#'   reference database for taxonomic assignment is a good fit to the data.
#'
#' @param ps The phyloseq object to operate on.
#' @param by_sample (Default: FALSE) Whether or not to return the calculations
#'   for the whole object, or broken out by individual samples in the dataset.
#'
#' @import phyloseq
#'
#' @return The percentage of unmapped reads, either for the overall dataset or
#'   on a per-sample basis.
#'
#' @export
#'

percent_unassigned <- function(ps,
                               by_sample = FALSE){

  # Get ASVs without assignment (all entries NA)
  taxtab <- data.frame(ps@tax_table)
  unassigned <- apply(taxtab, 1, function(x){all(is.na(x))})

  # Get proportion of reads from these ASVs in ASV table
    asvtab <- data.frame(ps@otu_table)

  if(by_sample == FALSE){
    sum(asvtab[, unassigned])/sum(asvtab)}
  else {
    # If only one ASV, need to calculate differently
    if (is.null(dim(unassigned))){
      asvtab[, unassigned]/rowSums(asvtab)
      # If more, vectorize
    } else {
      rowSums(asvtab[, unassigned])/rowSums(asvtab)
    }
  }

}
