#' @title Join ASV assignments to taxonomy from a taxmap object
#'
#' @description The output of the \code{assignSpecies_mod} function is a
#'   two-column dataframe of ASVs and their taxonomic assignments based on exact
#'   matching to reference sequences.  This joins the named taxonomic assignment
#'   to a full taxonomy queried from NCBI.
#'
#' @param taxmap An \code{R6Class} object of class \code{Taxmap} produced by
#'   the \code{metacoder::lookup_tax_data} function.
#'
#' @rawNamespace import(dplyr, except = id)
#'
#' @return A join of each ASV assignment to its full taxonomy, done from finer
#'   taxonomic levels upward.
#' @export
#'

asv_to_taxonomy = function(taxmap){
  #### Get taxonomy ####
  # Need this in two forms: one using Latin names, and one using taxon IDs
  # Named taxonomy
  taxonomy <- metacoder::taxonomy_table(taxmap,
                                        use_ranks = c('superkingdom',
                                                      'kingdom',
                                                      'phylum',
                                                      'order',
                                                      'family',
                                                      'genus',
                                                      'species',
                                                      'subspecies',
                                                      'varietas',
                                                      'forma'),
                                        add_id_col = TRUE)

  # ID taxonomy
  taxonomy.ids <- metacoder::taxonomy_table(taxmap,
                                            value = 'taxon_ids',
                                            use_ranks = c('superkingdom',
                                                          'kingdom',
                                                          'phylum',
                                                          'order',
                                                          'family',
                                                          'genus',
                                                          'species',
                                                          'subspecies',
                                                          'varietas',
                                                          'forma'))

  #### Join to taxmap ####

  taxtab <-
    full_join(taxmap$data$query_data,
              taxonomy)

  # There will be missing entries here if they are nested (e.g. Cucumis
  # sativus var. hardwickii will have an entry and C. sativus will not).

  # Fill these in using the ID-based taxonomy

  taxtab.add <-
    taxtab %>%
    filter(is.na(superkingdom) & taxon_id != 'unknown') %>%
    select(taxon_id:label)

  missing <- unique(taxtab.add$taxon_id) # Duplicates if multiple accessions

  # Find the row of each missing taxon in the ID-based taxonomy table
  # Helper function to get first index (all matches should be identical):
  first_index <- function(m){
    # m is a missing taxon ID
    which(taxonomy.ids == m, arr.ind = TRUE)[1, ]
  }

  rows <-
    lapply(missing, first_index) %>%
    bind_rows() %>%
    mutate(taxon_id = missing) %>%
    mutate(pad_na = ncol(taxonomy.ids) - col)

  # Now, extract each row, pad with NAs, and build out taxonomy entries to add
  # TODO: Vectorize this, do inside mutate?
  taxonomy.add <- data.frame()
  # Remove taxon ID to synchronize indices now that join has been done
  taxonomy <- taxonomy[, 2:ncol(taxonomy)]
  for (i in seq(nrow(rows))){
    entry <- taxonomy[rows$row[i], 1:rows$col[i]] # Extract
    entry <-
      entry %>%
      as.character() %>%
      c(rep(NA, rows$pad_na[i])) # Pad

    names(entry) <- names(taxonomy) # Name

    taxonomy.add <- bind_rows(taxonomy.add, entry) # Join
  }

  # Join back to ASVs by lowest-level label
  taxtab.add <-
    taxonomy.add %>%
    MButils::lowest_level() %>%
    left_join(taxtab.add, ., by = c('label' = 'name'))

  #### Make full taxonomy table ####

  taxtab <-
    taxtab %>%
    filter(!(label %in% taxtab.add$label)) %>%  # Remove entries we looked up
    bind_rows(taxtab.add) %>% # Rejoin
    select(asv, superkingdom:forma)

  taxtab
}
