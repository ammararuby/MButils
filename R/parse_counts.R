#' @title Parse read count QC files generated at pipeline steps
#'
#' @description This function parses text files of sample read counts of
#'   generated during the metabarcoding analysis pipeline for QC purposes.
#'
#' @usage parse_counts(f)
#'
#' @param f The path to a text file that alternates fastq.gz filenames with
#'   their read counts.
#'
#' @importFrom rlang :=
#'
#' @return A two-column data frame matching samples to their read counts.
#'
#' @export
#'

parse_counts <- function(f){
    counts.text <- utils::read.table(f)

  # Odd rows are sample names, evens are read counts
  n <- dim(counts.text)[1]

  # Get file basename
  name <- gsub('.txt', '', basename(f))

  counts.raw <-
    # Pull sample name lines
    data.frame(sample = counts.text$V1[seq(from = 1, to = n, by = 2)]) %>%
    # Remove Illumina suffix
    mutate(sample = gsub(pattern = '_.*', '', sample),
           # Pull count lines
           !!name := as.numeric(counts.text$V1[seq(from = 2, to = n,
                                                   by = 2)]))

  counts.raw
}
