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
                                                        'forma'))

    # Join the two together
    taxtab <- MButils::asv_to_taxonomy(taxtab.species,
                                       taxonomy)

    # Collapse multiple identifications to their last common ancestor
    cat('Calculating last common answer of matched species...\n')
    taxtab <- MButils::lca(taxtab)
  }

  if (use_function == 'taxonomy'){
    taxtab <-
      assignTaxonomy(asvtab,
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
