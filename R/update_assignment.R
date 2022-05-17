#' @title Update the taxonomic assignment of a phyloseq object
#'
#' @description Replaces the \code{tax_table} slot of a phyloseq object with
#'   updated assignments based on a revised version of a reference database.
#'
#' @param ps Phyloseq object to be updated
#' @param ref Path to the reference database (a FASTA file) to use for
#'   assignment
#' @param format One of either \code{'species'} or \code{'taxonomy'}, indicating
#'   if reference database is formatted for species (\code{assignSpecies}) or
#'   taxonomy (\code{assignTaxonomy}) assignment. Defaults to \code{'species'}
#'   with the MButils \code{assignSpecies_mod} function.
#'
#' @import phyloseq
#'
#' @return An updated phyloseq object.
#' @export
#'

update_assignment <- function(ps, ref, format = 'species'){
  # Get ASV table from object
  asvtab <- ps@otu_table@.Data

  # Do taxonomic assignment
  if (format == 'Species'){
    # Make assignment
    taxtab.species <- assignSpecies_mod(seqtab.merged,
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
    taxtab <- MButils::lca(taxtab)
  }

  if (format == 'taxonomy'){

  }

  # Replace in phyloseq object
  tax_table(ps) <- taxtab

}
