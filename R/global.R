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
  # lca
  utils::globalVariables("asv")
  # merge_into_phyloseq
  utils::globalVariables("row_names")

  invisible()
}
