#' @title Update the taxonomic assignment of a phyloseq object
#'
#' @description Replaces the \code{tax_table} slot of a phyloseq object with
#'   updated assignments based on a revised version of a reference database.
#'
#' @param ps Phyloseq object to be updated
#' @param ref Path to the reference database (a FASTA file) to use for
#'   assignment
#' @param use_function One of either \code{'species'} or \code{'taxonomy'},
#'   indicating if reference database is use_functionted for species
#'   (\code{assignSpecies}) or taxonomy (\code{assignTaxonomy}) assignment.
#'   Defaults to \code{'species'} with the MButils \code{assignSpecies_mod}
#'   function.
#'
#' @import phyloseq
#'
#' @return An updated phyloseq object.
#' @export
#'

update_assignment <- function(ps, ref, use_function = 'species'){
  # Get ASV table from object
  asvtab <- ps@otu_table@.Data

  # Do taxonomic assignment
  if (use_function == 'species'){
    # Make assignment
    cat('Matching ASVs to reference species...\n')
    taxtab.species <- MButils::assignSpecies_mod(asvtab,
                                                 refFasta = ref,
                                                 tryRC = TRUE)

    # Get full taxonomy of results
    # Separate ID from species name for querying
    taxtab.species <-
      taxtab.species %>%
      separate(col = Species,
               into = c('index', 'label'),
               sep = '\\s',
               extra = 'merge')

    taxmap <- metacoder::lookup_tax_data(taxtab.species,
                                         type = 'seq_id',
                                         column = 'index')

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

    # Join the two together
    taxtab <-
      full_join(taxmap$data$query_data,
                taxonomy)

    # There will be some missing entries here if they are nested (e.g. Cucumis
    # sativus var. hardwickii will have an entry and C. sativus will not). Fill
    # these in

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

    # Add back to full set of results
    taxtab <-
      taxtab %>%
      filter(!(label %in% taxtab.add$label)) %>%  # Remove entries we looked up
      bind_rows(taxtab.add) # Rejoin

    # Collapse multiple identifications to their last common ancestor
    cat('Calculating last common answer of matched species...\n')
    taxtab <- MButils::lca(taxtab)
  }

  if (use_function == 'taxonomy'){
    taxtab <-
      dada2::assignTaxonomy(asvtab,
                            taxLevels = c(
                              'kingdom',
                              'phylum',
                              'class',
                              'order',
                              'family',
                              'genus',
                              'species',
                              'subspecies'
                            ),
                            refFasta = ref,
                            tryRC = TRUE)
  }

  # Replace in phyloseq object
  tax_table(ps) <- as.matrix(taxtab)

  # Return updated object
  ps

}
