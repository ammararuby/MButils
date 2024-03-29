% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylosequtils.R
\name{merge_into_phyloseq}
\alias{merge_into_phyloseq}
\title{Add New Sample Data into Phyloseq Object}
\usage{
merge_into_phyloseq(dataframe, ps_object, variable)
}
\arguments{
\item{dataframe}{The dataframe containing new data to be incorporated.}

\item{ps_object}{The phyloseq object to be updated with new data.}

\item{variable}{The variable (column in data) as a string by which to merge
new data into phyloseq object. Variable name and values must be identical
between the dataframe and the phyloseq object}
}
\value{
An updated phyloseq object with the new data incorporated into sample
  data.
}
\description{
Allows incorporation of new sample data (a dataframe) to be
  merged into an existing phyloseq object.
}
