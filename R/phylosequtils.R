#' @title Add New Sample Data into Phyloseq Object
#'
#' @description Allows incorporation of new sample data (a dataframe) to be merged into an existing phyloseq object.
#'
#' @param dataframe The dataframe containing new data to be incorporated.
#' @param ps_object The phyloseq object to be updated with new data.
#' @param variable The variable (column in data) by which to merge new data into phyloseq object. Variable name and values must be identical between the dataframe and the phyloseq object
#'
#' @import dplyr
#' @import tibble
#' @import phyloseq
#'
#' @return An updated phyloseq object with the new data incorporated into sample data.
#' @export
#'
#'

merge_into_phyloseq = function(dataframe, ps_object, variable){

  sam= as.data.frame(sample_data(ps_object)) %>%
    rownames_to_column("row_names") %>%
    remove_rownames() %>%
    as_data_frame() %>%
    select(row_names, everything())

  complete_data = left_join(sam, dataframe, by= variable) %>%
    column_to_rownames("row_names")

  complete_data_ps = sample_data(complete_data)

  sample_data(ps_object) = complete_data_ps

  return(ps_object)
}
