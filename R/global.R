.onLoad <- function(libname, pkgname)
{
  # Make data set names global to avoid CHECK notes
  # assignSpecies_mod
  utils::globalVariables("Species")
  # asv_to_taxonomy
  utils::globalVariables("superkingdom")
  utils::globalVariables("kingdom")
  utils::globalVariables("phylum")
  utils::globalVariables("order")
  utils::globalVariables("family")
  utils::globalVariables("genus")
  utils::globalVariables("species")
  utils::globalVariables("subspecies")
  utils::globalVariables("varietas")
  utils::globalVariables("forma")
  utils::globalVariables("label")
  utils::globalVariables("taxon_id")
  utils::globalVariables("index")
  utils::globalVariables(".")
  # biplot
  utils::globalVariables("X1")
  utils::globalVariables("X2")
  utils::globalVariables("PC1")
  utils::globalVariables("PC2")
  utils::globalVariables("slope")
  utils::globalVariables("ang")
  utils::globalVariables("variable")
  # lca
  utils::globalVariables("asv")
  # lookup_tax_data
  utils::globalVariables("correct_taxon_names")
  # merge_into_phyloseq
  utils::globalVariables("row_names")

  invisible()
}
