% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_assignment.R
\name{update_assignment}
\alias{update_assignment}
\title{Update the taxonomic assignment of a phyloseq object}
\usage{
update_assignment(ps, ref, use_function = "species")
}
\arguments{
\item{ps}{Phyloseq object to be updated}

\item{ref}{Path to the reference database (a FASTA file) to use for
assignment}

\item{use_function}{One of either \code{'species'} or \code{'taxonomy'},
indicating if reference database is use_functionted for species
(\code{assignSpecies}) or taxonomy (\code{assignTaxonomy}) assignment.
Defaults to \code{'species'} with the MButils \code{assignSpecies_mod}
function.}
}
\value{
An updated phyloseq object.
}
\description{
Replaces the \code{tax_table} slot of a phyloseq object with
  updated assignments based on a revised version of a reference database.
}
