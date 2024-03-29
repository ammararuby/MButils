% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lca.R
\name{lca}
\alias{lca}
\title{Find the last common ancestor of assignSpecies taxonomic assignments}
\usage{
lca(assignment)
}
\arguments{
\item{assignment}{The output of \code{assignSpecies_mod} (*i.e.* matches of
ASVs to sequences in a metabarcoding reference database).}
}
\value{
A taxonomy table with taxonomic assignment per ASV, compatible for
  use in a phyloseq object.
}
\description{
The output of the \code{assignSpecies_mod} function is a
  two-column dataframe of ASVs and their taxonomic assignments based on exact
  matching to reference sequences: if multiple sequences are matched, the ASV
  will be repeated across rows, one per match.  This function collapses each
  ASV to a single row labeled with the last common ancestor of all its
  matched taxa (*e.g.*, an ASV matching to an ASV matching to both wheat
  [Triticum aestivum] and rye [Secale cereale] was relabeled as Poaceae, the
  family shared by both genera).
}
