#' Example ASV assignments from assignSpecies_mod
#'
#' \describe{
#' \item{asv}{Amplicon sequence variant (ASV) identified}
#' \item{index}{The NCBI accession number of the sequence's exact match in the
#' reference database}
#' \item{label}{The taxonomic name associated with the NCBI accession number}
#' }
#'
"assignment"

#' Example taxonomy data structure
#'
"taxonomy"

#' Example sequence table
#'
#' Rows are samples (named here by their plate and well location), and columns
#' are *trnL* amplicon sequence variants (ASVs). Values indicate counts of an
#' ASV detected in a given sample inferred by DADA2.
#'
#' @format ## `seqtab` A data frame with 52 rows and 110 columns:
#' \describe{
#'   \item{ATCACGTTTTCCGAAAACAAACAACGGTTCAGAAAGCGAAAATCAAAAAG}{ASV1}
#'   \item{ATCCCTTTTTTGAAAAACAAGTGGTTCTCAAACTAGAACCCAAAGGAAAAG}{ASV2}
#'   ...
#' }
"seqtab"
