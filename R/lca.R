#' @title Find the last common ancestor of assignSpecies taxonomic assignments
#'
#' @description The output of the \code{assignSpecies_mod} function is a
#'   two-column dataframe of ASVs and their taxonomic assignments based on exact
#'   matching to reference sequences: if multiple sequences are matched, the ASV
#'   will be repeated across rows, one per match.  This function collapses each
#'   ASV to a single row labeled with the last common ancestor of all its
#'   matched taxa (*e.g.*, an ASV matching to an ASV matching to both wheat
#'   [Triticum aestivum] and rye [Secale cereale] was relabeled as Poaceae, the
#'   family shared by both genera).
#'
#' @param assignment The output of \code{assignSpecies_mod} (*i.e.* matches of
#'   ASVs to sequences in a metabarcoding reference database).
#'
#' @return A taxonomy table with taxonomic assignment per ASV, compatible for
#'   use in a phyloseq object.
#' @export
#'

lca <- function(assignment){
  # Group assignments by ASV
  lca <-
    assignment %>%
    group_by(asv) %>%
    summarize_all(n_distinct,
                  # Want to keep NAs in case one match is specified to
                  # a different level than another
                  na.rm = FALSE) %>%
    column_to_rownames(var = 'asv')

  # Now, relabel all those with >1 name at a particular level as NA
  # As a placeholder, keep only the first species, knowing it will be overwritten
  assignment.lca <-
    assignment %>%
    group_by(asv) %>%
    summarize_all(first) %>%
    column_to_rownames(var = 'asv')

  # Relabel
  assignment.lca[lca > 1] = NA

  assignment.lca
}
