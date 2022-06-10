#' @title Join ASV assignments to taxonomy
#'
#' @description The output of the \code{assignSpecies_mod} function is a
#'   two-column dataframe of ASVs and their taxonomic assignments based on exact
#'   matching to reference sequences.  This joins the named taxonomic assignment
#'   to a full taxonomy queried from NCBI, assuming that all assignments occur
#'   to at least the species level, and potentially to subspecies levels
#'   (currently hard-coded for subspecies, variety, and form).
#'
#' @param assignment The output of \code{assignSpecies_mod}
#' @param taxonomy An \code{R6Class} object of class \code{Taxonomy} produced by
#'   the \code{metacoder::taxonomy} function.
#'
#' @rawNamespace import(dplyr, except = id)
#'
#' @return A join of each ASV assignment to its full taxonomy, done from finer
#'   taxonomic levels upward.
#' @export
#'

asv_to_taxonomy = function(assignment, taxonomy){
  # If present, remove '(nom. inval.)' notation, which interferes with joins
  taxonomy <-
    mutate(taxonomy,
           across(.fns = ~gsub(pattern = ' \\(nom\\. inval\\.\\)$',
                               replacement = '',
                               .x)))

  # The joins need to be done at sub-species level designations, which in rank
  # order are:
  # - subspecies
  # - varietas
  # - forma

  # TODO Do join based on taxid, rather than level by level?

  #### Separate taxonomy table into these components ####
  # Form
  taxonomy.f <-
    filter(taxonomy, !is.na(forma))
  # Variety
  taxonomy.var <-
    filter(taxonomy, !is.na(varietas) & is.na(forma))
  # Subspecies
  taxonomy.ssp <-
    filter(taxonomy, !is.na(subspecies) & is.na(varietas) & is.na(forma))
  # Species
  # For simplicity, just make a list of all species-level entries; this will
  # ultimately be used in a right join so okay if there are some that won't be
  # incorporated because they're assigned to a finer level
  taxonomy.spp <-
    taxonomy %>%
    select(superkingdom:species) %>%
    distinct()

  #### Join each to name, accession, and ASV sequence from assignments ####
  # Form
  taxtab.f <-
    assignment %>%
    left_join(taxonomy.f, by = c('label' = 'forma')) %>%
    filter(!is.na(species)) %>%
    select(asv,
           superkingdom,
           kingdom,
           phylum,
           order,
           family,
           genus,
           species,
           subspecies,
           varietas,
           forma = label)
  # Varietas
  taxtab.var <-
    assignment %>%
    left_join(taxonomy.var, by = c('label' = 'varietas')) %>%
    filter(!is.na(species)) %>%
    select(asv,
           superkingdom,
           kingdom,
           phylum,
           order,
           family,
           genus,
           species,
           subspecies,
           varietas = label,
           forma)
  # Subspecies
  taxtab.ssp <-
    assignment %>%
    left_join(taxonomy.ssp, by = c('label' = 'subspecies')) %>%
    filter(!is.na(species)) %>%
    select(asv,
           superkingdom,
           kingdom,
           phylum,
           order,
           family,
           genus,
           species,
           subspecies = label,
           varietas,
           forma)
  # Species
  taxtab.spp <-
    assignment %>%
    filter(!(asv %in% c(taxtab.f$asv,
                        taxtab.var$asv,
                        taxtab.ssp$asv))) %>%
    left_join(taxonomy.spp, by = c('label' = 'species')) %>%
    select(asv,
           superkingdom,
           kingdom,
           phylum,
           order,
           family,
           genus,
           species = label) %>%
    # Fill missing columns with NA
    mutate(subspecies = NA,
           varietas = NA,
           forma = NA) %>%
    distinct()

  #### Join all together ####
  taxtab <-
    bind_rows(taxtab.spp,
              taxtab.ssp,
              taxtab.var,
              taxtab.f)

  taxtab
}
